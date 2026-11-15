#include "pipeline.h"

#include <cmath>
#include <set>

namespace rawfuse {
namespace {

struct Observation {
    size_t index;
    double offset, rateSignal, excess, sigma;
};

std::optional<Observation> observe(const RawFrame& frame, size_t i,
                                   const std::array<float, 2>& model, const Options& options) {
    if (frame.rejected[i] & (saturated | bloom | badRect | outsideValidity)) return std::nullopt;
    const auto& layout = frame.layout;
    const int y = int(i / layout.width), x = int(i % layout.width);
    const double full = frame.white - frame.black(y, x);
    std::vector<double> neighbours;
    neighbours.reserve(8);
    // Use the same Bayer phase and compare black-subtracted sensor counts.
    for (int dy = -2; dy <= 2; dy += 2)
        for (int dx = -2; dx <= 2; dx += 2) {
            const int row = y + dy, col = x + dx;
            if ((!dy && !dx) || row < 0 || col < 0 || row >= layout.height || col >= layout.width)
                continue;
            if (frame.rejected[layout.index(row, col)] & (saturated | bloom | excluded | outsideValidity)) continue;
            const double black = frame.black(row, col);
            const double signal = frame.raw(row, col) - black;
            if (signal <= options.high * (frame.white - black)) neighbours.push_back(signal);
        }
    if (neighbours.size() < 3) return std::nullopt;

    const double local = median(neighbours);
    for (double& v : neighbours) v = std::abs(v - local);
    const double texture = 1.4826 * median(neighbours);
    const auto noise = frame.noise[layout.colour(y, x)];
    const double net = frame.raw(y, x) - frame.black(y, x);
    const double variance = full * (noise.shot * std::max(0.0, net)) + noise.read * full * full;
    const double sigma = std::max(1.0, std::hypot(std::sqrt(1.25 * variance), texture));
    return Observation{i, model[0], model[1] * frame.metadata.seconds, net - local, sigma};
}

struct NormalEquations {
    double xx = 0, yy = 0, xy = 0, xz = 0, yz = 0;

    void add(const Observation& obs, double weight) {
        const double x = obs.offset, y = obs.rateSignal, z = obs.excess;
        xx += weight * x * x;
        yy += weight * y * y;
        xy += weight * x * y;
        xz += weight * x * z;
        yz += weight * y * z;
    }

    std::array<double, 2> solve(std::string& method) const {
        const double determinant = xx * yy - xy * xy;
        if (xx == 0 && yy == 0) return {0, 0};
        if (xx == 0 || yy == 0) {
            method = xx == 0 ? "rate_only" : "offset_only";
            return {xx ? std::max(0.0, xz / xx) : 0, yy ? std::max(0.0, yz / yy) : 0};
        }
        if (determinant <= 1e-4 * xx * yy) {
            const double norm = xx + 2 * xy + yy;
            if (norm <= 1e-4 * (xx + yy)) return {0, 0};
            method = "shared_scale";
            const double scale = std::max(0.0, (xz + yz) / norm);
            return {scale, scale};
        }
        method = "separate_scales";
        const double a = (xz * yy - yz * xy) / determinant;
        const double b = (yz * xx - xz * xy) / determinant;
        if (a < 0) return {0, std::max(0.0, yz / yy)};
        if (b < 0) return {std::max(0.0, xz / xx), 0};
        return {a, b};
    }
};

} // namespace

DarkFit fitDarkModel(std::vector<RawFrame>& frames, const DarkModel& model, const Options& o) {
    DarkFit fit;
    std::vector<std::vector<Observation>> byFrame(frames.size());
    std::set<size_t> informativeSites;
    for (size_t k = 0; k < frames.size(); ++k)
        for (const auto& entry : model) {
            const auto obs = observe(frames[k], entry.first, entry.second, o);
            if (!obs) continue;
            byFrame[k].push_back(*obs);
            ++fit.observations;
            if (std::max(std::abs(obs->offset), std::abs(obs->rateSignal)) >= obs->sigma)
                informativeSites.insert(obs->index);
        }

    fit.sites = informativeSites.size();
    for (auto& frame : frames) frame.hotRejected = model.size();
    if (fit.sites < 8) return fit;

    // Normalize leverage by read noise, local texture and predicted defect amplitude.
    // Huber weights keep scene edges and failed defect predictions from driving the fit.
    std::array<double, 2> scale{1, 1};
    auto predict = [&](const Observation& obs) { return scale[0] * obs.offset + scale[1] * obs.rateSignal; };
    for (int pass = 0; pass < 12; ++pass) {
        NormalEquations equations;
        for (const auto& samples : byFrame)
            for (const auto& obs : samples) {
                const double modelMagnitude = std::abs(obs.offset) + std::abs(obs.rateSignal);
                const double spread = std::hypot(obs.sigma, 0.15 * modelMagnitude);
                const double residual = std::abs(obs.excess - predict(obs));
                const double robust = residual > 1.5 * spread ? 1.5 * spread / residual : 1;
                equations.add(obs, robust / (spread * spread));
            }
        scale = equations.solve(fit.method);
    }
    if (fit.method == "unavailable") return fit;
    if (fit.method != "rate_only") fit.offsetScale = scale[0];
    if (fit.method != "offset_only") fit.rateScale = scale[1];

    std::vector<double> relativeResiduals, relativeNoise;
    for (const auto& samples : byFrame)
        for (const auto& obs : samples) {
            const double prediction = predict(obs);
            if (std::abs(prediction) < 3 * obs.sigma) continue;
            relativeResiduals.push_back((obs.excess - prediction) / std::abs(prediction));
            relativeNoise.push_back(obs.sigma / std::abs(prediction));
        }
    if (!relativeResiduals.empty()) {
        const double centre = median(relativeResiduals);
        for (double& v : relativeResiduals) v = std::abs(v - centre);
        const double scatter = std::hypot(centre, 1.4826 * median(relativeResiduals));
        const double noise = median(relativeNoise);
        fit.relativeScatter = std::sqrt(std::max(0.0, scatter * scatter - noise * noise));
    }
    fit.uncertainty = std::max(o.darkModelError, fit.relativeScatter.value_or(0));

    // A small, repeatable mismatch may be invisible in every individual exposure.
    // Pool its significance across the sweep before admitting that model site.
    struct Evidence { double mismatch = 0, strength = 0; };
    std::map<size_t, Evidence> evidence;
    for (const auto& samples : byFrame)
        for (const auto& obs : samples) {
            const double prediction = predict(obs);
            auto& sum = evidence[obs.index];
            sum.mismatch += prediction * (obs.excess - prediction) / (obs.sigma * obs.sigma);
            sum.strength += (prediction / obs.sigma) * (prediction / obs.sigma);
        }
    const double relativeTolerance = std::min(0.5, 3 * fit.uncertainty);
    for (auto it = evidence.begin(); it != evidence.end();) {
        const auto [mismatch, strength] = it->second;
        if (std::abs(mismatch) > std::max(3 * std::sqrt(strength), relativeTolerance * strength))
            ++it;
        else it = evidence.erase(it);
    }
    for (size_t k = 0; k < frames.size(); ++k) {
        auto& frame = frames[k];
        frame.darkUncertainty = fit.uncertainty;
        for (const auto& obs : byFrame[k]) {
            const double prediction = predict(obs);
            const double tolerance = std::max(4 * obs.sigma, relativeTolerance * std::abs(prediction));
            if (evidence.count(obs.index) || std::abs(obs.excess - prediction) > tolerance) continue;
            frame.dark.emplace_back(obs.index, float(prediction));
            frame.rejected[obs.index] = (frame.rejected[obs.index] & ~excluded) | modeled;
            ++frame.hotKept;
            --frame.hotMissing;
            --frame.hotRejected;
        }
    }
    return fit;
}
}

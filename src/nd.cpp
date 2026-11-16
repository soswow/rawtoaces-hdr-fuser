#include "pipeline.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace rawfuse {
namespace {

constexpr int blockSize = 64;
constexpr size_t minimumSamples = 2000;
// The reference merge needs SNR ≥ 10 at a calibration site. An attenuated sample
// is selected by the signal the reference predicts for it (SNR ≥ 5), not by its
// own noisy reading, which would favour samples that happened to read high.
// Thousands of such samples still give a median to about 0.01 EV. Drift samples
// need SNR ≥ 5, because thousands of them still reveal a large scene change.
constexpr double referenceVariance = 0.01, sampleVariance = 0.04, driftVariance = 0.04;
struct Samples {
    std::array<std::vector<float>, 3> ev;
    std::vector<float> drift;
};

std::string seconds(double value) {
    std::ostringstream text;
    if (value < 0.5) text << "1/" << std::setprecision(4) << 1 / value << " s";
    else text << std::setprecision(3) << value << " s";
    return text.str();
}

std::string overlapError(const std::vector<RawFrame>& reference, const std::vector<RawFrame>& filtered,
                         const NDResult& estimate) {
    std::ostringstream message;
    message << "insufficient overlap with the calibration reference (R/G/B " << estimate.samples[0] << '/'
            << estimate.samples[1] << '/' << estimate.samples[2] << " usable sites; need " << minimumSamples
            << " per colour). A usable site is unclipped with SNR of at least 10 in the reference merge, "
               "and an attenuated frame records it unclipped with a predicted SNR of at least 5";
    double strongest = -1;
    for (int c = 0; c < 3; ++c)
        if (estimate.samples[c] >= 50) strongest = std::max(strongest, estimate.measured[c]);
    if (strongest >= 0) {
        // The longest attenuated frame records what this unattenuated exposure would.
        const double equivalent = filtered.back().effective * std::exp2(-strongest);
        const double gap = std::log2(reference.front().effective / equivalent);
        message << ". At the measured " << std::fixed << std::setprecision(2) << strongest
                << " EV, the longest attenuated frame equals an unattenuated " << seconds(equivalent)
                << ", " << std::fabs(gap) << " EV " << (gap > 0 ? "darker than" : "brighter than")
                << " the shortest reference frame (" << seconds(reference.front().effective) << ")";
        if (gap > 0) message << "; the brackets do not overlap, so extend the attenuated bracket or "
                                "shorten the reference bracket by at least that much";
        else message << "; the brackets overlap, but too little of the scene falls in the shared range";
    }
    message << ". Widen the calibration_mask or set a measured ev for this group";
    return message.str();
}

double floatMedian(const std::vector<float>& input) {
    std::vector<double> values(input.begin(), input.end());
    return median(values);
}

double relativeVariance(const RawFrame& frame, int colour, double fraction, double extra) {
    const auto noise = frame.noise[colour];
    return (noise.shot * fraction + noise.read + extra) / (fraction * fraction);
}

} // namespace

NDResult calibrateND(const MergeResult& merged, const std::vector<RawFrame>& reference,
                     const std::vector<RawFrame>& filtered, const Mask& calibration, const Mask& drift,
                     const std::optional<RGB>& manualEV, const Options& o, Warnings& warnings) {
    const auto& layout = reference.front().layout;
    const int across = (layout.width + blockSize - 1) / blockSize;
    const int down = (layout.height + blockSize - 1) / blockSize;
    const bool checkDrift = !drift.empty();
    // Odd stride samples every Bayer phase while keeping the overlap pool bounded on large images.
    int stride = 1;
    constexpr size_t sampleBudget = 4'000'000;
    while (double(layout.size()) * filtered.size() / (double(stride) * stride) > sampleBudget &&
           stride < 7) stride += 2;
    std::vector<Samples> blocks(size_t(across) * down);
    // Without `predicted` attenuation, samples are selected by their own reading for a first estimate.
    auto gather = [&](int stride, const std::array<bool, 3>& include, bool collectDrift,
                      const std::optional<RGB>& predicted) {
        parallelRows(down, o.threads, [&](int by) {
            for (int bx = 0; bx < across; ++bx) {
                auto& block = blocks[size_t(by) * across + bx];
                for (int c = 0; c < 3; ++c) if (include[c]) block.ev[c].clear();
                if (collectDrift) block.drift.clear();
                for (int y = by * blockSize; y < std::min((by + 1) * blockSize, layout.height); y += stride) {
                    for (int x = bx * blockSize; x < std::min((bx + 1) * blockSize, layout.width); x += stride) {
                        const size_t i = layout.index(y, x);
                        const int colour = layout.colour(y, x);
                        const bool covered = calibration.empty() || calibration[i];
                        const bool unfiltered = collectDrift && drift[i];
                        if (!(covered && include[colour]) && !unfiltered) continue;
                        const auto& measurement = merged.sites[i];
                        if (measurement.flags & (lowerBound | missing | moving | adjacentClipFallback) ||
                            measurement.signal <= 0)
                            continue;
                        const double reliable = measurement.variance /
                                                (double(measurement.signal) * measurement.signal);
                        if (covered && include[colour] && reliable < referenceVariance) {
                            for (const auto& frame : filtered) {
                                const auto s = sample(frame, y, x, i, o);
                                if (s.state != SampleState::usable || s.fraction <= 0) continue;
                                const double expected = predicted ?
                                    frame.effective * measurement.signal * std::exp2(-(*predicted)[colour]) : s.fraction;
                                if (expected < o.low || expected > 0.8 * o.high ||
                                    relativeVariance(frame, colour, expected, s.extra) >= sampleVariance)
                                    continue;
                                constexpr size_t maxBlockSamples = 2048;
                                auto& values = block.ev[colour];
                                if (values.size() < maxBlockSamples)
                                    values.push_back(float(std::log2(frame.effective * measurement.signal / s.fraction)));
                            }
                        }
                        if (unfiltered && reliable < driftVariance) {
                            // Drift sites lie outside the group's validity, so read them directly.
                            for (const auto& frame : filtered) {
                                if (frame.rejected[i] & excluded) continue;
                                const float black = frame.black(y, x);
                                const double full = frame.white - black;
                                const double raw = frame.raw(y, x);
                                if (raw >= frame.white - o.clipMargin * full) continue;
                                const double fraction = (raw - black) / full;
                                if (fraction < o.low || fraction > o.high ||
                                    relativeVariance(frame, colour, fraction, 0) >= driftVariance)
                                    continue;
                                block.drift.push_back(float(std::log2(fraction / (frame.effective * measurement.signal))));
                            }
                        }
                    }
                }
            }
        });
    };
    gather(stride, {true, true, true}, checkDrift, std::nullopt);

    NDResult result;
    result.measured.fill(std::numeric_limits<double>::quiet_NaN());
    result.mad = result.measured;
    std::vector<double> changes;
    for (const auto& b : blocks) changes.insert(changes.end(), b.drift.begin(), b.drift.end());
    if (checkDrift && stride > 1 && changes.size() < 100) {
        gather(1, {}, true, std::nullopt);
        changes.clear();
        for (const auto& b : blocks) changes.insert(changes.end(), b.drift.begin(), b.drift.end());
    }
    result.sceneDriftSamples = changes.size();
    std::string driftWarning;
    if (checkDrift && changes.size() < 100)
        warnings.push_back("too few unfiltered samples (" + std::to_string(changes.size()) +
                           ") to measure scene drift");
    if (changes.size() >= 100) {
        result.sceneDrift = median(changes);
        if (std::abs(*result.sceneDrift) > 0.2) {
            std::ostringstream message;
            message << "scene drift of " << std::fixed << std::setprecision(2) << *result.sceneDrift
                    << " EV exceeds 0.20 EV in the unfiltered area";
            driftWarning = message.str();
            if (manualEV && !o.allowSceneDrift)
                throw std::runtime_error(driftWarning +
                    "; capture matching brackets or pass --allow-scene-drift to keep the manual ev");
            warnings.push_back(driftWarning);
        }
    }
    auto estimate = [&](int c) {
        result.samples[c] = 0;
        std::vector<double> values;
        size_t count = 0;
        for (const auto& b : blocks) count += b.ev[c].size();
        values.reserve(count);
        for (const auto& b : blocks) values.insert(values.end(), b.ev[c].begin(), b.ev[c].end());
        if (!values.empty()) {
            const double centre = median(values);
            std::vector<double> deviations(values.size());
            std::transform(values.begin(), values.end(), deviations.begin(),
                           [centre](double value) { return std::abs(value - centre); });
            const double cutoff = std::max(0.1, 4 * median(deviations));
            values.erase(std::remove_if(values.begin(), values.end(), [&](double value) {
                return std::abs(value - centre) > cutoff;
            }), values.end());
            result.samples[c] = values.size();
            if (!values.empty()) result.measured[c] = median(values);
            if (values.size() >= minimumSamples) {
                deviations.resize(values.size());
                std::transform(values.begin(), values.end(), deviations.begin(),
                               [&](double value) { return std::abs(value - result.measured[c]); });
                result.mad[c] = median(deviations);
            }
        }
    };
    // A first estimate from self-selected samples, then the final one from samples
    // selected by their predicted signal.
    RGB rough{};
    for (int c = 0; c < 3; ++c) estimate(c);
    {
        std::vector<double> known;
        for (int c = 0; c < 3; ++c) if (result.samples[c] >= 20) known.push_back(result.measured[c]);
        const double fallback = known.empty() ? std::numeric_limits<double>::quiet_NaN() : median(known);
        for (int c = 0; c < 3; ++c)
            rough[c] = manualEV ? (*manualEV)[c] : result.samples[c] >= 20 ? result.measured[c] : fallback;
    }
    if (std::all_of(rough.begin(), rough.end(), [](double ev) { return std::isfinite(ev); })) {
        gather(stride, {true, true, true}, false, rough);
        std::array<bool, 3> retry{};
        for (int c = 0; c < 3; ++c) {
            estimate(c);
            retry[c] = result.samples[c] < minimumSamples;
        }
        if (stride > 1 && std::any_of(retry.begin(), retry.end(), [](bool value) { return value; })) {
            gather(1, retry, false, rough);
            for (int c = 0; c < 3; ++c) if (retry[c]) estimate(c);
        }
    }
    const bool enough = std::all_of(result.samples.begin(), result.samples.end(),
        [](size_t count) { return count >= minimumSamples; });

    if (manualEV) {
        result.ev = *manualEV;
        for (int c = 0; c < 3; ++c)
            if (result.samples[c] >= minimumSamples && std::abs(result.measured[c] - (*manualEV)[c]) > 0.2) {
                warnings.push_back("measured attenuation differs from the supplied ev by more than 0.2 EV");
                break;
            }
    } else {
        if (!enough) {
            const auto error = overlapError(reference, filtered, result);
            throw std::runtime_error(driftWarning.empty() ? error : error + ". " + driftWarning);
        }
        result.ev = result.measured;
        for (const double ev : result.ev)
            if (ev < -0.2 || ev > 64)
                throw std::runtime_error("measured attenuation is outside the expected range; check masks and scene overlap");
    }
    for (int c = 0; c < 3; ++c)
        if (result.samples[c] < minimumSamples) result.measured[c] = std::numeric_limits<double>::quiet_NaN();

    size_t flagged = 0;
    for (size_t i = 0; i < blocks.size(); ++i) {
        const auto& b = blocks[i];
        const size_t total = b.ev[0].size() + b.ev[1].size() + b.ev[2].size();
        if (total < 50) continue;
        BlockCheck block{int(i / across) * blockSize, int(i % across) * blockSize, {}, {}, false};
        block.ev.fill(std::numeric_limits<double>::quiet_NaN());
        for (int c = 0; c < 3; ++c) {
            block.samples[c] = b.ev[c].size();
            if (b.ev[c].empty()) continue;
            block.ev[c] = floatMedian(b.ev[c]);
            if (block.samples[c] >= 12 && std::abs(block.ev[c] - result.ev[c]) > 0.25) block.flagged = true;
        }
        flagged += block.flagged;
        result.blocks.push_back(block);
    }
    if (flagged)
        warnings.push_back(std::to_string(flagged) + " calibration blocks differ from the group's ev by more "
                           "than 0.25 EV (glare, filter edge or scene change); see blocks_flagged");
    return result;
}
}

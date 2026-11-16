#include "pipeline.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace rawfuse {
namespace {

constexpr size_t minimumRatioSamples = 2000;
constexpr double maximumCorrectionEv = 0.25;

// Odd strides visit every Bayer phase while bounding the work per frame pair.
int strideFor(const Layout& layout, double budget) {
    int stride = 1;
    while (double(layout.size()) / (double(stride) * stride) > budget && stride < 15) stride += 2;
    return stride;
}

// Calls visit(xa, xb) for sites that both frames record as plain usable samples.
template<class Visit>
void forUsablePairs(const RawFrame& a, const RawFrame& b, const Options& o, int stride, Visit visit) {
    const auto& layout = a.layout;
    for (int y = 0; y < layout.height; y += stride)
        for (int x = 0; x < layout.width; x += stride) {
            const size_t i = layout.index(y, x);
            const auto sa = sample(a, y, x, i, o);
            if (sa.state != SampleState::usable || sa.extra != 0) continue;
            const auto sb = sample(b, y, x, i, o);
            if (sb.state != SampleState::usable || sb.extra != 0) continue;
            visit(double(sa.fraction), double(sb.fraction));
        }
}

} // namespace

double robustVariance(std::vector<double> values) {
    auto copy = values;
    const double centre = median(copy);
    for (auto& v : copy) v = std::abs(v - centre);
    double sigma = 1.4826 * median(copy);
    if (!(sigma > 0)) {
        double sum = 0;
        for (double v : values) sum += (v - centre) * (v - centre);
        sigma = std::sqrt(sum / double(values.size()));
    }
    double variance = sigma * sigma;
    for (int pass = 0; pass < 3 && variance > 0; ++pass) {
        const double limit = 16 * variance;
        double sum = 0, mean = 0;
        size_t count = 0;
        for (double v : values)
            if ((v - centre) * (v - centre) <= limit) { mean += v; ++count; }
        if (count < 2) break;
        mean /= double(count);
        for (double v : values)
            if ((v - centre) * (v - centre) <= limit) sum += (v - mean) * (v - mean);
        variance = sum / double(count - 1);
    }
    return variance;
}

namespace {

double fullScale(const RawFrame& frame) {
    double black = 0;
    for (float value : frame.blackPattern) black += value;
    return frame.white - black / double(frame.blackPattern.size());
}

double opticalBlackRead(const RawFrame& frame) {
    if (!(frame.opticalBlackSigma > 0)) return 0;
    const double sigma = frame.opticalBlackSigma / fullScale(frame);
    return sigma * sigma;
}

struct NoiseRow {
    double iso, x, r, variance, readA, readB;
    size_t samples;
};

constexpr int tileSize = 8;
constexpr size_t minimumTileSamples = 48, minimumTiles = 50, minimumCellTiles = 20;
constexpr double percentile = 0.1;

double lowPercentile(std::vector<double> values) {
    auto it = values.begin() + std::ptrdiff_t(percentile * double(values.size() - 1));
    std::nth_element(values.begin(), it, values.end());
    return *it;
}

// The 10th percentile of n-sample variances of Gaussian noise, relative to the
// true variance (Wilson–Hilferty approximation to the chi-squared quantile).
double percentileBias(double samples) {
    const double k = samples - 1, z = -1.2816;
    return std::pow(1 - 2 / (9 * k) + z * std::sqrt(2 / (9 * k)), 3);
}

void fitLeastSquares(const std::vector<NoiseRow>& rows, NoiseFit& fit) {
    if (rows.size() < 4) return;
    // Relative weights: each point has a similar fractional uncertainty.
    double aa = 0, ab = 0, bb = 0, av = 0, bv = 0;
    for (const auto& row : rows) {
        const double w = 1 / (row.variance * row.variance);
        const double a = row.x * (1 + row.r), b = 1 + row.r * row.r;
        aa += w * a * a; ab += w * a * b; bb += w * b * b;
        av += w * a * row.variance; bv += w * b * row.variance;
    }
    const double det = aa * bb - ab * ab;
    double shot = det > 0 ? (av * bb - bv * ab) / det : 0;
    double read = det > 0 ? (aa * bv - ab * av) / det : 0;
    if (!(read > 0)) {
        read = 0;
        shot = av / aa;
    }
    fit.shot = shot;
    fit.read = std::max(read, shot * 1e-4);
    fit.bins = rows.size();
}

std::string format(double value) {
    std::ostringstream text;
    text << std::setprecision(3) << value;
    return text.str();
}

std::vector<ExposureStep> measureSteps(const std::vector<RawFrame>&, const std::string& group,
                                       const Options&, Warnings&);

} // namespace

std::vector<std::vector<ExposureStep>> measureExposures(const FrameGroups& groups, const Options& o,
                                                        Warnings& warnings) {
    std::vector<std::vector<ExposureStep>> result(groups.size());
    if (o.exposures != "measured") return result;
    for (size_t g = 0; g < groups.size(); ++g)
        result[g] = measureSteps(*groups[g], o.groups[g].name, o, warnings);
    // Shutter error is a property of the camera at each setting, so a step that
    // cannot be measured takes the correction measured for the same pair of
    // settings in another group, typically a plain bracket for an ND bracket.
    for (size_t g = 0; g < groups.size(); ++g) {
        size_t unmeasured = 0;
        for (auto& step : result[g]) {
            if (step.applied) continue;
            for (size_t d = 0; d < groups.size() && step.borrowedFrom.empty(); ++d)
                for (const auto& donor : result[d])
                    if (d != g && donor.applied && donor.borrowedFrom.empty() &&
                        donor.fromSeconds == step.fromSeconds && donor.toSeconds == step.toSeconds &&
                        std::abs(donor.exifRatio / step.exifRatio - 1) < 1e-6) {
                        step.ratio = step.exifRatio * donor.measuredRatio / donor.exifRatio;
                        step.borrowedFrom = o.groups[d].name;
                        break;
                    }
            unmeasured += step.borrowedFrom.empty();
        }
        auto& frames = *groups[g];
        std::vector<double> chain(frames.size(), 1);
        for (size_t k = 1; k < frames.size(); ++k) chain[k] = chain[k - 1] * result[g][k - 1].ratio;
        // Keep the EXIF geometric mean so the absolute scale follows the metadata.
        double logExif = 0, logChain = 0;
        for (size_t k = 0; k < frames.size(); ++k) {
            logExif += std::log(frames[k].exifEffective);
            logChain += std::log(chain[k]);
        }
        const double anchor = std::exp((logExif - logChain) / double(frames.size()));
        for (size_t k = 0; k < frames.size(); ++k) frames[k].effective = chain[k] * anchor;
        if (unmeasured)
            warnings.push_back(o.groups[g].name + ": " + std::to_string(unmeasured) + " of " +
                               std::to_string(result[g].size()) + " exposure steps could not be measured; "
                               "using EXIF ratios for them");
    }
    return result;
}

namespace {

std::vector<ExposureStep> measureSteps(const std::vector<RawFrame>& frames, const std::string& group,
                                       const Options& o, Warnings& warnings) {
    std::vector<ExposureStep> steps;
    for (size_t k = 0; k + 1 < frames.size(); ++k) {
        const auto& a = frames[k];
        const auto& b = frames[k + 1];
        ExposureStep step;
        step.from = a.inputIndex + 1;
        step.to = b.inputIndex + 1;
        step.exifRatio = b.effective / a.effective;
        step.fromSeconds = a.metadata.seconds;
        step.toSeconds = b.metadata.seconds;
        step.ratio = step.exifRatio;
        // Both readings well above the noise and below the nonlinear top of the range.
        const double minimum = std::max(0.03, o.low);
        std::vector<double> ratios;
        auto collect = [&](int stride) {
            ratios.clear();
            forUsablePairs(a, b, o, stride, [&](double xa, double xb) {
                if (xa >= minimum) ratios.push_back(xb / xa);
            });
        };
        const int stride = strideFor(a.layout, 1'500'000);
        collect(stride);
        // Short frames of a bright scene can overlap only near highlights.
        if (ratios.size() < minimumRatioSamples && stride > 1) collect(1);
        step.samples = ratios.size();
        if (ratios.size() >= minimumRatioSamples) {
            step.measuredRatio = median(ratios);
            const double correction = std::log2(step.measuredRatio / step.exifRatio);
            if (step.measuredRatio > 0 && std::abs(correction) <= maximumCorrectionEv) {
                step.ratio = step.measuredRatio;
                step.applied = true;
            } else {
                std::ostringstream message;
                message << group << ": measured exposure ratio of " << a.path.filename().string() << " to "
                        << b.path.filename().string() << " differs from EXIF by " << std::fixed
                        << std::setprecision(2) << correction << " EV, more than a shutter error; "
                        << "the frames likely disagree (motion, an occluder or flare). Not using it";
                warnings.push_back(message.str());
            }
        }
        steps.push_back(step);
    }
    return steps;
}

} // namespace

std::vector<NoiseFit> assignNoise(const FrameGroups& groups, const Options& o, Warnings& warnings) {
    std::vector<RawFrame*> pending;
    for (auto* group : groups)
        for (auto& frame : *group) {
            const auto& profile = o.noiseOverride.empty() ? frame.metadata.noise : o.noiseOverride;
            if (profile.empty()) {
                pending.push_back(&frame);
                continue;
            }
            if (profile.size() != 2 && profile.size() != 6)
                throw std::runtime_error(frame.path.string() + ": expected one or three NoiseProfile pairs");
            for (size_t c = 0; c < 3; ++c) {
                const size_t index = profile.size() == 6 ? 2 * c : 0;
                frame.noise[c] = {profile[index], profile[index + 1]};
                if (frame.noise[c].shot < 0 || frame.noise[c].read < 0 ||
                    (frame.noise[c].shot == 0 && frame.noise[c].read == 0))
                    throw std::runtime_error("noise coefficients must be nonnegative and include a positive term");
            }
            frame.noiseSource = o.noiseOverride.empty() ? "dng" : "custom";
        }
    if (pending.empty()) return {};

    // Temporal differences of neighbouring frames cancel the scene: with
    // r = E_b / E_a, var(x_b − r·x_a) = S·x_b·(1 + r) + O_b + r²·O_a.
    // Motion, misregistration and texture only add variance, so the estimates
    // use a low percentile of small-tile variances, corrected for its sampling bias.
    std::vector<NoiseRow> tiles;
    std::vector<size_t> pairOf;
    size_t pairs = 0;
    const double lowest = 1e-3;
    for (auto* group : groups) {
        for (size_t k = 0; k + 1 < group->size(); ++k) {
            const auto& a = (*group)[k];
            const auto& b = (*group)[k + 1];
            if (a.metadata.iso != b.metadata.iso || !a.noiseSource.empty() || !b.noiseSource.empty())
                continue;
            const double r = b.effective / a.effective;
            const auto& layout = a.layout;
            const int stride = strideFor(layout, 1'500'000);
            const int span = tileSize * stride;
            std::vector<double> difference;
            for (int top = 0; top < layout.height; top += span)
                for (int left = 0; left < layout.width; left += span) {
                    difference.clear();
                    double level = 0;
                    for (int y = top; y < std::min(top + span, layout.height); y += stride)
                        for (int x = left; x < std::min(left + span, layout.width); x += stride) {
                            const size_t i = layout.index(y, x);
                            const auto sa = sample(a, y, x, i, o), sb = sample(b, y, x, i, o);
                            if (sa.state != SampleState::usable || sb.state != SampleState::usable ||
                                sa.extra != 0 || sb.extra != 0)
                                continue;
                            difference.push_back(sb.fraction - r * sa.fraction);
                            level += sb.fraction;
                        }
                    if (difference.size() < minimumTileSamples) continue;
                    level /= double(difference.size());
                    if (level < lowest) continue;
                    tiles.push_back({a.metadata.iso, level, r, robustVariance(difference),
                                     opticalBlackRead(a), opticalBlackRead(b), difference.size()});
                    pairOf.push_back(pairs);
                }
            ++pairs;
        }
    }

    std::vector<NoiseFit> fits;
    std::vector<double> isos;
    for (const auto& tile : tiles)
        if (std::find(isos.begin(), isos.end(), tile.iso) == isos.end()) isos.push_back(tile.iso);
    for (const double iso : isos) {
        std::vector<size_t> subset;
        for (size_t t = 0; t < tiles.size(); ++t) if (tiles[t].iso == iso) subset.push_back(t);
        NoiseFit fit{iso};
        const bool known = std::all_of(subset.begin(), subset.end(), [&](size_t t) {
            return tiles[t].readA > 0 && tiles[t].readB > 0;
        });
        auto readOf = [&](const NoiseRow& tile) { return tile.readB + tile.r * tile.r * tile.readA; };
        if (known) {
            // Read noise from optical black; S from shot-dominated tiles.
            auto shotOf = [&](const NoiseRow& tile) {
                return (tile.variance - readOf(tile)) / (tile.x * (1 + tile.r));
            };
            auto estimate = [&](auto keep) {
                std::vector<double> values;
                double samples = 0;
                for (size_t t : subset)
                    if (keep(tiles[t])) {
                        values.push_back(shotOf(tiles[t]));
                        samples += double(tiles[t].samples);
                    }
                if (values.size() < minimumTiles) return std::make_pair(0.0, size_t(0));
                return std::make_pair(lowPercentile(values) / percentileBias(samples / double(values.size())),
                                      values.size());
            };
            auto first = estimate([](const NoiseRow& tile) { return tile.x >= 0.02; });
            if (first.second) {
                auto refined = estimate([&](const NoiseRow& tile) {
                    return first.first * tile.x * (1 + tile.r) >= 2 * readOf(tile);
                });
                if (!refined.second) refined = first;
                fit.shot = refined.first;
                fit.bins = refined.second;
                fit.opticalBlack = true;
                double read = 0;
                for (size_t t : subset) read += (tiles[t].readA + tiles[t].readB) / 2;
                fit.read = read / double(subset.size());
            }
        } else {
            // Per pair and level, the low-percentile variance becomes one point of a
            // weighted least-squares fit of S and O.
            constexpr int bins = 24;
            const double step = std::log(o.high / lowest) / bins;
            std::map<std::pair<size_t, int>, std::vector<size_t>> cells;
            for (size_t t : subset)
                cells[{pairOf[t], std::min(bins - 1, int(std::log(tiles[t].x / lowest) / step))}].push_back(t);
            std::vector<NoiseRow> rows;
            for (const auto& cell : cells) {
                if (cell.second.size() < minimumCellTiles) continue;
                std::vector<double> variances, levels;
                double samples = 0;
                for (size_t t : cell.second) {
                    variances.push_back(tiles[t].variance);
                    levels.push_back(tiles[t].x);
                    samples += double(tiles[t].samples);
                }
                const auto& any = tiles[cell.second.front()];
                rows.push_back({iso, median(levels), any.r,
                                lowPercentile(variances) / percentileBias(samples / double(variances.size())),
                                0, 0, cell.second.size()});
            }
            fitLeastSquares(rows, fit);
        }
        if (std::isfinite(fit.shot) && fit.shot >= 1e-8 && fit.shot <= 1e-2 &&
            std::isfinite(fit.read) && fit.read <= 1e-3)
            fits.push_back(fit);
    }

    size_t generic = 0;
    const RawFrame* firstGeneric = nullptr;
    for (auto* frame : pending) {
        const double iso = frame->metadata.iso;
        const NoiseFit* nearest = nullptr;
        for (const auto& fit : fits)
            if (!nearest || std::abs(std::log(fit.iso / iso)) < std::abs(std::log(nearest->iso / iso)))
                nearest = &fit;
        const double ob = opticalBlackRead(*frame);
        Noise noise;
        if (nearest) {
            const double scale = iso / nearest->iso;
            noise = {nearest->shot * scale, ob > 0 ? ob : nearest->read * scale};
            frame->noiseSource = "measured";
        } else {
            // Typical of current sensors at base ISO, within a few times.
            noise = {3e-5 * iso / 100, ob > 0 ? ob : 2e-8 * iso / 100};
            frame->noiseSource = "generic";
            if (!generic++) firstGeneric = frame;
        }
        frame->noise.fill(noise);
    }
    if (generic) {
        const auto& noise = firstGeneric->noise[0];
        warnings.push_back("no NoiseProfile and too little bracket overlap to measure noise for " +
                           std::to_string(generic) + " frames; using a generic profile (S=" + format(noise.shot) +
                           ", O=" + format(noise.read) + " at the first such frame's ISO)");
    }
    return fits;
}
}

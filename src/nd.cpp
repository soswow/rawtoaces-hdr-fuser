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
struct Samples {
    std::array<std::vector<float>, 3> ev;
    std::vector<float> drift;
};

std::string overlapError(const std::vector<RawFrame>& plain, const std::vector<RawFrame>& nd,
                         const NDResult& estimate, const Options& o) {
    const double shortestPlain = plain.front().effective;
    const double longestND = nd.back().effective;
    const double overlapEv = std::log2(longestND * o.high / (shortestPlain * o.ndMinSignal));
    std::ostringstream message;
    message << "insufficient ND overlap (R/G/B " << estimate.samples[0] << '/' << estimate.samples[1]
            << '/' << estimate.samples[2] << " inliers; need " << minimumSamples << " per colour). "
            << "Longest ND exposure is " << nd.back().metadata.exposure
            << " s, shortest plain exposure is " << plain.front().metadata.exposure
            << " s. The exposure windows can overlap up to " << std::fixed << std::setprecision(2)
            << overlapEv << " EV of attenuation. Capture a longer ND exposure with usable scene overlap "
            << (o.groups.empty() ? "or supply --nd-ev EV" : "or supply a measured ev for this group");
    double measured = -1;
    for (double ev : estimate.measured) if (std::isfinite(ev)) measured = std::max(measured, ev);
    if (measured >= 0) {
        const double required = shortestPlain * std::exp2(measured) * o.ndMinSignal / o.high;
        message << ". At the strongest measured channel (" << measured << " EV), overlap requires at least "
                << std::defaultfloat << std::setprecision(3)
                << required * nd.back().metadata.seconds / longestND << " s at this ND ISO/aperture";
    }
    return message.str();
}

double floatMedian(const std::vector<float>& input) {
    std::vector<double> values(input.begin(), input.end());
    return median(values);
}

} // namespace

NDResult calibrateND(const MergeResult& plainMerge, const std::vector<RawFrame>& plain,
                     const std::vector<RawFrame>& nd, const Mask& coverage,
                     const std::optional<RGB>& manualEV, bool checkUncovered,
                     const Options& o, Warnings& warnings) {
    const auto& layout = plain.front().layout;
    const int across = (layout.width + blockSize - 1) / blockSize;
    const int down = (layout.height + blockSize - 1) / blockSize;
    auto behindFilter = [&](size_t i) { return coverage.empty() || coverage[i]; };
    // Odd stride samples every Bayer phase while keeping the overlap pool bounded on large images.
    int stride = 1;
    constexpr size_t sampleBudget = 4'000'000;
    while (double(layout.size()) * nd.size() / (double(stride) * stride) > sampleBudget &&
           stride < 7) stride += 2;
    std::vector<Samples> blocks(size_t(across) * down);
    auto gather = [&](int stride, const std::array<bool, 3>& include, bool collectDrift) {
        parallelRows(down, o.threads, [&](int by) {
            for (int bx = 0; bx < across; ++bx) {
                auto& block = blocks[size_t(by) * across + bx];
                for (int c = 0; c < 3; ++c) if (include[c]) block.ev[c].clear();
                if (collectDrift) block.drift.clear();
                for (int y = by * blockSize; y < std::min((by + 1) * blockSize, layout.height); y += stride) {
                    for (int x = bx * blockSize; x < std::min((bx + 1) * blockSize, layout.width); x += stride) {
                        const size_t i = layout.index(y, x);
                        const auto& measurement = plainMerge.sites[i];
                        const bool covered = behindFilter(i);
                        if (covered ? !include[layout.colour(y, x)] : !collectDrift) continue;
                        // Thousands of uncovered samples can reveal large scene changes
                        // even when each comparison is less precise than an ND calibration sample.
                        const double maxRelativeVariance = covered ? 0.0025 : 0.04;
                        if (measurement.flags & (lowerBound | missing | moving | adjacentClipFallback) ||
                            measurement.signal <= 0 ||
                            measurement.variance >= maxRelativeVariance * measurement.signal * measurement.signal)
                            continue;
                        const float black = plain.front().black(y, x);
                        const int colour = layout.colour(y, x);
                        for (const auto& frame : nd) {
                            const auto s = sample(frame, y, x, i, black, o);
                            const double minSignal = covered ? o.ndMinSignal : o.low;
                            if (s.state != SampleState::usable || s.fraction <= 0 || s.fraction < minSignal)
                                continue;
                            if (!covered) {
                                const auto noise = frame.noise[colour];
                                const double variance = noise.shot * s.fraction + noise.read + s.extra;
                                if (variance >= maxRelativeVariance * s.fraction * s.fraction) continue;
                            }
                            const float ev = float(std::log2(frame.effective * measurement.signal / s.fraction));
                            if (covered) {
                                constexpr size_t maxBlockSamples = 2048;
                                auto& values = block.ev[colour];
                                if (values.size() < maxBlockSamples) values.push_back(ev);
                            }
                            else block.drift.push_back(-ev);
                        }
                    }
                }
            }
        });
    };
    gather(stride, {true, true, true}, checkUncovered);

    NDResult result;
    result.measured.fill(std::numeric_limits<double>::quiet_NaN());
    result.mad = result.measured;
    std::vector<double> drift;
    for (const auto& b : blocks) drift.insert(drift.end(), b.drift.begin(), b.drift.end());
    if (checkUncovered && stride > 1 && !coverage.empty() && drift.size() < 100) {
        gather(1, {}, true);
        drift.clear();
        for (const auto& b : blocks) drift.insert(drift.end(), b.drift.begin(), b.drift.end());
    }
    result.sceneDriftSamples = drift.size();
    std::string driftWarning;
    if (drift.size() >= 100) {
        result.sceneDrift = median(drift);
        if (std::abs(*result.sceneDrift) > 0.2) {
            std::ostringstream message;
            message << "ND scene drift of " << std::fixed << std::setprecision(2) << *result.sceneDrift
                    << " EV exceeds 0.20 EV in the uncovered area";
            driftWarning = message.str();
            if (manualEV && !o.allowNDSceneDrift)
                throw std::runtime_error(driftWarning +
                    "; capture matching sweeps or use --allow-nd-scene-drift with manual --nd-ev");
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
            if (values.size() >= minimumSamples) {
                result.measured[c] = median(values);
                deviations.resize(values.size());
                std::transform(values.begin(), values.end(), deviations.begin(),
                               [&](double value) { return std::abs(value - result.measured[c]); });
                result.mad[c] = median(deviations);
            }
        }
    };
    std::array<bool, 3> retry{};
    for (int c = 0; c < 3; ++c) {
        estimate(c);
        retry[c] = result.samples[c] < minimumSamples;
    }
    if (stride > 1 && std::any_of(retry.begin(), retry.end(), [](bool value) { return value; })) {
        gather(1, retry, false);
        for (int c = 0; c < 3; ++c) if (retry[c]) estimate(c);
    }
    const bool enough = std::all_of(result.samples.begin(), result.samples.end(),
        [](size_t count) { return count >= minimumSamples; });

    if (manualEV) {
        result.ev = *manualEV;
        for (int c = 0; c < 3; ++c)
            if (std::isfinite(result.measured[c]) && std::abs(result.measured[c] - (*manualEV)[c]) > 0.2) {
                warnings.push_back("measured ND strength differs from supplied ND EV by more than 0.2 EV");
                break;
            }
    } else {
        if (!enough) {
            const auto error = overlapError(plain, nd, result, o);
            throw std::runtime_error(driftWarning.empty() ? error : error + ". " + driftWarning);
        }
        result.ev = result.measured;
        for (const double ev : result.ev)
            if (ev < -0.2 || ev > 64)
                throw std::runtime_error("measured ND transmission is outside the expected range; check coverage and scene overlap");
    }

    result.rejectedBlocks.resize(blocks.size());
    if (o.aux && o.groups.empty()) result.blockImage.assign(layout.size(), 0);
    for (size_t i = 0; i < blocks.size(); ++i) {
        const auto& b = blocks[i];
        const size_t total = b.ev[0].size() + b.ev[1].size() + b.ev[2].size();
        if (total < 50) continue;
        BlockCheck block{int(i / across) * blockSize, int(i % across) * blockSize, {}, {}, false};
        block.ev.fill(std::numeric_limits<double>::quiet_NaN());
        double deviation = 0;
        for (int c = 0; c < 3; ++c) {
            block.samples[c] = b.ev[c].size();
            if (b.ev[c].empty()) continue;
            block.ev[c] = floatMedian(b.ev[c]);
            if (!std::isfinite(result.measured[c]) || block.samples[c] < 12) continue;
            const double candidate = block.ev[c] - result.measured[c];
            if (std::abs(candidate) > std::abs(deviation)) deviation = candidate;
        }
        block.removed = std::abs(deviation) > 0.25;
        if (block.removed)
            result.rejectedBlocks[i] = 1;
        if (!result.blockImage.empty())
            for (int y = block.row; y < std::min(block.row + blockSize, layout.height); ++y)
                for (int x = block.col; x < std::min(block.col + blockSize, layout.width); ++x)
                    if (behindFilter(layout.index(y, x))) result.blockImage[layout.index(y, x)] = float(deviation);
        result.blocks.push_back(block);
    }
    return result;
}

NDResult mergeND(MergeResult& plainMerge, const std::vector<RawFrame>& plain,
                 const std::vector<RawFrame>& nd, const Mask& coverage, const Options& o,
                 Warnings& warnings) {
    const auto& layout = plain.front().layout;
    const int across = (layout.width + blockSize - 1) / blockSize;
    std::optional<RGB> manual;
    if (o.ndEv) manual = RGB{*o.ndEv, *o.ndEv, *o.ndEv};
    auto result = calibrateND(plainMerge, plain, nd, coverage, manual, !coverage.empty(), o, warnings);
    const RGB transmission{std::exp2(-result.ev[0]), std::exp2(-result.ev[1]), std::exp2(-result.ev[2])};
    const bool fillOnly = o.ndUse == "fill";
    const std::vector<RawFrame> noPlain;
    std::atomic<size_t> blocked{0};
    std::atomic<int64_t> shadowDelta{0};
    parallelRows(layout.height, o.threads, [&](int y) {
        size_t blockedRow = 0;
        int shadowRowDelta = 0;
        for (int x = 0; x < layout.width; ++x) {
            const size_t i = layout.index(y, x);
            auto& site = plainMerge.sites[i];
            if (!coverage.empty() && !coverage[i]) continue;
            if ((fillOnly && !(site.flags & lowerBound)) || (!fillOnly && (site.flags & moving)))
                continue;
            if (result.rejectedBlocks[size_t(y / blockSize) * across + x / blockSize]) {
                blockedRow += (site.flags & lowerBound) != 0;
                continue;
            }
            const auto merged = mergeSite(fillOnly ? noPlain : plain, &nd, transmission, y, x, o);
            auto replacement = merged.site;
            if (!(replacement.flags & fromND)) continue;
            if (site.flags & lowerBound) {
                if (replacement.flags & adjacentClipFallback) continue;
                if (!(replacement.flags & lowerBound) && replacement.signal < 0.9f * site.signal) {
                    site.flags |= ndInconsistent;
                    continue;
                }
                if ((replacement.flags & lowerBound) && replacement.signal <= site.signal) continue;
                replacement.signal = std::max(replacement.signal, site.signal);
            }
            replacement.flags |= site.flags & moving;
            site = replacement;
            if (!plainMerge.shadowSites.empty()) {
                shadowRowDelta += int(merged.shadow) - int(plainMerge.shadowSites[i]);
                plainMerge.shadowSites[i] = merged.shadow;
            }
        }
        blocked.fetch_add(blockedRow, std::memory_order_relaxed);
        shadowDelta.fetch_add(shadowRowDelta, std::memory_order_relaxed);
    });
    result.blockedLowerBounds = blocked;
    plainMerge.shadow = size_t(int64_t(plainMerge.shadow) + shadowDelta.load());
    return result;
}
}

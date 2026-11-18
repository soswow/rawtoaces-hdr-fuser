#include "pipeline.h"

#include <cmath>
#include <limits>

namespace rawfuse {
namespace {

struct Measurement {
    const RawFrame* frame;
    double exposure;
    float fraction, extra, spread = 1, detail = 0;
};

double variance(const Noise& noise, double fraction, double extra, double spread = 1) {
    return std::max(1e-16, spread * (noise.shot * std::max(fraction, 0.0) + noise.read) + extra);
}

Site one(const Measurement& m, int colour) {
    const double rate = m.fraction / m.exposure;
    const double var = variance(m.frame->noise[colour], m.fraction, m.extra, m.spread) / (m.exposure * m.exposure);
    return {float(rate), float(var), 0, 1};
}

double quad(const RawFrame& frame, int row, int col, const Options& o) {
    const auto& layout = frame.layout;
    double sum = 0;
    for (int y = row; y < row + 2; ++y) {
        for (int x = col; x < col + 2; ++x) {
            const auto s = sample(frame, y, x, layout.index(y, x), o);
            if (s.state != SampleState::usable || s.fraction < o.low) return 0;
            sum += s.fraction;
        }
    }
    return sum / (4 * frame.effective);
}

} // namespace

size_t findReference(const std::vector<RawFrame>& frames, const Options& o) {
    if (o.reference) {
        const auto it = std::find_if(frames.begin(), frames.end(), [&](const RawFrame& f) {
            return f.inputIndex + 1 == size_t(o.reference);
        });
        return size_t(it - frames.begin());
    }
    const auto& layout = frames.front().layout;
    size_t selected = (frames.size() - 1) / 2;
    size_t maximum = 0;
    for (size_t i = 0; i < frames.size(); ++i) {
        size_t count = 0;
        for (int y = 0; y + 1 < layout.height; y += 16)
            for (int x = 0; x + 1 < layout.width; x += 16)
                count += quad(frames[i], y, x, o) > 0;
        if (count > maximum || (count == maximum && i == (frames.size() - 1) / 2)) {
            maximum = count;
            selected = i;
        }
    }
    return selected;
}

namespace {

// The merge of one site or pixel colour over a bracket; read(frame) gives each
// frame's sample there. Frames are weighted by the inverse of their noise at a
// common seed rate, so the weights do not follow the noise in each reading.
//
// An interpolated sample (a moved frame, or a colour the pixel does not
// measure) also has an interpolation error, which depends on where its sites
// land and on the scene's detail there. The detail D (rate² per unit of 1 − Σw²)
// comes from the least noisy interpolated frame: the spread of its sites,
// Σw·(f − f̄)², beyond interpolationGate times what their noise alone gives,
// (1 − Σw²)·σ²; noise alone passes that gate about 0.2% of the time for two
// sites. A frame's error is then (1 − Σw²)·D, zero for a frame whose site lies
// on the pixel. Frames whose sites land alike share that error, and weighting
// cannot remove it, so it must not shift weight between them: the noise weight
// is multiplied by 1 / (1 + error / σ²ₘ), with σ²ₘ the noise of the merge. Where
// detail is far above the noise, a frame that measured the pixel dominates;
// in flat or noise-limited areas the weights stay those of the noise. Without
// interpolated frames the merge is the plain inverse-noise one. The variance
// reported is noise only.
constexpr double interpolationGate = 10;

template<class Read>
SiteMerge accumulate(const std::vector<RawFrame>& frames, int colour, Read read, const Options& o) {
    double seed = 0, seedWeight = 0;
    double bound = -1;  // Strongest clipped lower bound, in signal per unit exposure
    std::optional<Measurement> aboveHigh, adjacent;
    size_t longest = 0;
    thread_local std::vector<Measurement> usable;
    usable.clear();
    usable.reserve(frames.size());
    for (const auto& frame : frames) {
        const auto s = read(frame);
        const double exposure = frame.effective;
        if (s.state == SampleState::clipped) bound = std::max(bound, std::max(0.0, double(s.fraction)) / exposure);
        if (s.state == SampleState::aboveHigh && (!aboveHigh || exposure < aboveHigh->exposure))
            aboveHigh = {&frame, exposure, s.fraction, s.extra, s.spread, s.detail};
        if (s.state == SampleState::adjacentClip && (!adjacent || exposure < adjacent->exposure))
            adjacent = {&frame, exposure, s.fraction, s.extra, s.spread, s.detail};
        if (s.state != SampleState::usable) continue;
        usable.push_back({&frame, exposure, s.fraction, s.extra, s.spread, s.detail});
        const double weight = std::max(0.0, std::min(s.fraction - o.low, o.high - s.fraction));
        seed += weight * s.fraction / exposure;
        seedWeight += weight;
        if (usable.size() == 1 || exposure > usable[longest].exposure) longest = usable.size() - 1;
    }
    if (usable.empty()) {
        if (aboveHigh) return {one(*aboveHigh, colour), false};
        if (adjacent) {
            auto fallback = one(*adjacent, colour);
            fallback.flags |= adjacentClipFallback;
            return {fallback, adjacent->fraction < o.low};
        }
        if (bound >= 0) return {{float(bound), 0, lowerBound, 0}, false};
        return {{0, 0, missing, 0}, false};
    }
    if (usable.size() == 1)
        return {one(usable.front(), colour), usable.front().fraction < o.low};
    const auto& reach = usable[longest];
    seed = std::max(0.0, seedWeight > 0 ? seed / seedWeight : reach.fraction / reach.exposure);
    // Scene detail from the least noisy interpolated frame, and the noise of a
    // noise-weighted merge, both in rate² units.
    double detail = 0, best = INFINITY, precision = 0;
    for (const auto& m : usable) {
        const auto noise = m.frame->noise[colour];
        const double var = variance(noise, seed * m.exposure, m.extra, m.spread) / (m.exposure * m.exposure);
        precision += 1 / var;
        const double open = 1 - m.spread;
        if (open <= 1e-6 || var >= best) continue;
        best = var;
        const double expected = interpolationGate * open * variance(noise, seed * m.exposure, 0);
        detail = std::max(0.0, double(m.detail) - expected) / (open * m.exposure * m.exposure);
    }
    double sum = 0, total = 0, noise2 = 0;
    bool shadow = true;
    for (const auto& m : usable) {
        const auto noise = m.frame->noise[colour];
        const double var = variance(noise, seed * m.exposure, m.extra, m.spread);
        const double geometry = detail > 0 ? 1 / (1 + (1 - m.spread) * detail * precision) : 1;
        const double scaledWeight = geometry * m.exposure / var;
        sum += scaledWeight * m.fraction;
        total += scaledWeight * m.exposure;
        noise2 += scaledWeight * scaledWeight * var;
        shadow &= m.fraction < o.low;
    }
    const double merged = detail > 0 ? noise2 / (total * total) : 1 / total;
    return {{float(sum / total), float(merged), 0, uint16_t(usable.size())}, shadow};
}

} // namespace

SiteMerge mergeSite(const std::vector<RawFrame>& frames, int row, int col, const Options& o) {
    const auto& layout = frames.front().layout;
    const size_t index = layout.index(row, col);
    return accumulate(frames, layout.colour(row, col),
                      [&](const RawFrame& frame) { return sample(frame, row, col, index, o); }, o);
}

Mask motionQuads(const std::vector<RawFrame>& frames, size_t reference, const Options& o) {
    const auto& layout = frames.front().layout;
    const auto& ref = frames[reference];
    const int quadWidth = (layout.width + 1) / 2, quadHeight = (layout.height + 1) / 2;
    Mask motion;
    if (!o.deghost || frames.size() < 3) return motion;
    motion.resize(size_t(quadWidth) * quadHeight);
    const double threshold = std::exp2(o.deghostEv);
    parallelRows(layout.height / 2, o.threads, [&](int y) {
        for (int x = 0; x < layout.width / 2; ++x) {
            const double value0 = quad(ref, 2 * y, 2 * x, o);
            if (value0 <= 0) continue;
            for (const auto& frame : frames) {
                if (&frame == &ref) continue;
                const auto value = quad(frame, 2 * y, 2 * x, o);
                if (value > value0 * threshold || (value > 0 && value * threshold < value0)) {
                    motion[size_t(y) * quadWidth + x] = 1;
                    break;
                }
            }
        }
    });
    return motion;
}

bool nearMotion(const Mask& motion, const Layout& layout, int row, int col) {
    if (motion.empty()) return false;
    const int quadWidth = (layout.width + 1) / 2, quadHeight = (layout.height + 1) / 2;
    for (int y = std::max(0, row / 2 - 1); y <= std::min(quadHeight - 1, row / 2 + 1); ++y)
        for (int x = std::max(0, col / 2 - 1); x <= std::min(quadWidth - 1, col / 2 + 1); ++x)
            if (motion[size_t(y) * quadWidth + x]) return true;
    return false;
}

MergeResult mergePlain(const std::vector<RawFrame>& frames, const Options& o, const Mask& selection,
                       bool trackShadow) {
    const auto& layout = frames.front().layout;
    MergeResult result;
    result.reference = findReference(frames, o);
    const auto& ref = frames[result.reference];
    const auto motion = motionQuads(frames, result.reference, o);
    result.sites.resize(layout.size());
    if (trackShadow) result.shadowSites.resize(layout.size());
    std::atomic<size_t> shadow{0};
    parallelRows(layout.height, o.threads, [&](int row) {
        size_t shadows = 0;
        for (int col = 0; col < layout.width; ++col) {
            const auto index = layout.index(row, col);
            if (!selection.empty() && !selection[index]) {
                result.sites[index].flags = missing | uncovered;
                continue;
            }
            const bool isMoving = nearMotion(motion, layout, row, col);
            if (isMoving) {
                const auto s = sample(ref, row, col, index, o);
                if (s.state == SampleState::usable) {
                    result.sites[index] = one({&ref, ref.effective, s.fraction, s.extra, s.spread}, layout.colour(row, col));
                    result.sites[index].flags |= moving;
                    shadows += s.fraction < o.low;
                    if (trackShadow) result.shadowSites[index] = s.fraction < o.low;
                    continue;
                }
            }
            auto value = mergeSite(frames, row, col, o);
            if (isMoving) value.site.flags |= moving;
            result.sites[index] = value.site;
            shadows += value.shadow;
            if (trackShadow) result.shadowSites[index] = value.shadow;
        }
        shadow.fetch_add(shadows, std::memory_order_relaxed);
    });
    result.shadow = shadow;
    return result;
}

MergeResult mergeSplat(const std::vector<RawFrame>& frames, const Options& o, const Mask& selection,
                       bool trackShadow) {
    const auto& layout = frames.front().layout;
    MergeResult result;
    result.reference = findReference(frames, o);
    const auto& ref = frames[result.reference];
    const auto motion = motionQuads(frames, result.reference, o);
    // Offsets, within a pixel's quad, of a site of each colour: it fixes the
    // lattice a colour the pixel does not measure is read on.
    std::array<Point, 3> anchor{};
    for (int k = 3; k >= 0; --k) anchor[size_t(layout.colour(k / 2, k % 2))] = {k / 2, k % 2};
    // Sites are read relative to the gain each full-colour site carries, which
    // the flat field applies after the merge.
    const auto* gains = ref.gains.get();
    result.sites.resize(layout.size() * 3);
    if (trackShadow) result.shadowSites.resize(layout.size() * 3);
    std::atomic<size_t> shadow{0};
    parallelRows(layout.height, o.threads, [&](int row) {
        size_t shadows = 0;
        std::vector<std::array<double, 2>> position(frames.size());
        for (int col = 0; col < layout.width; ++col) {
            const auto index = layout.index(row, col);
            if (!selection.empty() && !selection[index]) {
                for (int c = 0; c < 3; ++c) result.sites[index * 3 + c].flags = missing | uncovered;
                continue;
            }
            for (size_t k = 0; k < frames.size(); ++k) {
                const auto d = frames[k].motion ? frames[k].motion->at(row, col) : std::array<double, 2>{};
                position[k] = {row + d[0], col + d[1]};
            }
            const bool isMoving = nearMotion(motion, layout, row, col);
            const int own = layout.colour(row, col);
            for (int c = 0; c < 3; ++c) {
                const size_t j = index * 3 + size_t(c);
                const Point site = c == own ? Point{row, col}
                                            : Point{(row & ~1) + anchor[size_t(c)].row, (col & ~1) + anchor[size_t(c)].col};
                const double gain = gains ? latticeGain(*gains, layout, row, col, c) : 0;
                auto read = [&](const RawFrame& frame) {
                    if (c == own && !frame.motion) return sample(frame, row, col, index, o);
                    const auto& p = position[size_t(&frame - frames.data())];
                    return sampleAt(frame, p[0], p[1], site.row, site.col, o, gain);
                };
                if (isMoving) {
                    const auto s = read(ref);
                    if (s.state == SampleState::usable) {
                        result.sites[j] = one({&ref, ref.effective, s.fraction, s.extra, s.spread}, c);
                        result.sites[j].flags |= moving;
                        shadows += s.fraction < o.low;
                        if (trackShadow) result.shadowSites[j] = s.fraction < o.low;
                        continue;
                    }
                }
                auto value = accumulate(frames, c, read, o);
                if (isMoving) value.site.flags |= moving;
                result.sites[j] = value.site;
                shadows += value.shadow;
                if (trackShadow) result.shadowSites[j] = value.shadow;
            }
        }
        shadow.fetch_add(shadows, std::memory_order_relaxed);
    });
    result.shadow = shadow;
    return result;
}
}

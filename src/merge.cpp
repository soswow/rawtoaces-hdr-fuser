#include "pipeline.h"

#include <cmath>
#include <limits>

namespace rawfuse {
namespace {

struct Measurement {
    const RawFrame* frame;
    double exposure;
    float fraction, extra, spread = 1;
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

SiteMerge mergeSite(const std::vector<RawFrame>& frames, int row, int col, const Options& o) {
    const auto& first = frames.front();
    const int colour = first.layout.colour(row, col);
    const size_t index = first.layout.index(row, col);
    double seed = 0, seedWeight = 0;
    double bound = -1;  // Strongest clipped lower bound, in signal per unit exposure
    std::optional<Measurement> aboveHigh, adjacent;
    size_t longest = 0;
    thread_local std::vector<Measurement> usable;
    usable.clear();
    usable.reserve(frames.size());
    for (const auto& frame : frames) {
        const auto s = sample(frame, row, col, index, o);
        const double exposure = frame.effective;
        if (s.state == SampleState::clipped) bound = std::max(bound, std::max(0.0, double(s.fraction)) / exposure);
        if (s.state == SampleState::aboveHigh && (!aboveHigh || exposure < aboveHigh->exposure))
            aboveHigh = {&frame, exposure, s.fraction, s.extra, s.spread};
        if (s.state == SampleState::adjacentClip && (!adjacent || exposure < adjacent->exposure))
            adjacent = {&frame, exposure, s.fraction, s.extra, s.spread};
        if (s.state != SampleState::usable) continue;
        usable.push_back({&frame, exposure, s.fraction, s.extra, s.spread});
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
    double sum = 0, total = 0;
    bool shadow = true;
    for (const auto& m : usable) {
        const auto noise = m.frame->noise[colour];
        const double scaledWeight = m.exposure / variance(noise, seed * m.exposure, m.extra, m.spread);
        sum += scaledWeight * m.fraction;
        total += scaledWeight * m.exposure;
        shadow &= m.fraction < o.low;
    }
    return {{float(sum / total), float(1 / total), 0, uint16_t(usable.size())}, shadow};
}

MergeResult mergePlain(const std::vector<RawFrame>& frames, const Options& o, const Mask& selection,
                       bool trackShadow) {
    const auto& layout = frames.front().layout;
    MergeResult result;
    result.reference = findReference(frames, o);
    const auto& ref = frames[result.reference];
    const int quadWidth = (layout.width + 1) / 2, quadHeight = (layout.height + 1) / 2;
    std::vector<uint8_t> motion;
    if (o.deghost && frames.size() >= 3) {
        motion.resize(size_t(quadWidth) * quadHeight);
        const double threshold = std::exp2(o.deghostEv);
        parallelRows(layout.height / 2, o.threads, [&](int y) {
            for (int x = 0; x < layout.width / 2; ++x) {
                const double reference = quad(ref, 2 * y, 2 * x, o);
                if (reference <= 0) continue;
                for (const auto& frame : frames) {
                    if (&frame == &ref) continue;
                    const auto value = quad(frame, 2 * y, 2 * x, o);
                    if (value > reference * threshold || (value > 0 && value * threshold < reference)) {
                        motion[size_t(y) * quadWidth + x] = 1;
                        break;
                    }
                }
            }
        });
    }
    auto nearMotion = [&](int row, int col) {
        for (int y = std::max(0, row / 2 - 1); y <= std::min(quadHeight - 1, row / 2 + 1); ++y)
            for (int x = std::max(0, col / 2 - 1); x <= std::min(quadWidth - 1, col / 2 + 1); ++x)
                if (motion[size_t(y) * quadWidth + x]) return true;
        return false;
    };
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
            const bool isMoving = !motion.empty() && nearMotion(row, col);
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
}

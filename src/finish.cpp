#include "pipeline.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace rawfuse {
namespace {

float interpolateGain(const GainMap& map, int row, int col, const Layout& layout) {
    // DNG GainMap samples the center of each active-area pixel.
    const double y = std::clamp(((row + 0.5) / layout.height - map.originV) / map.spacingV,
                                0.0, double(map.rows - 1));
    const double x = std::clamp(((col + 0.5) / layout.width - map.originH) / map.spacingH,
                                0.0, double(map.cols - 1));
    const int y0 = int(y), y1 = std::min(y0 + 1, map.rows - 1);
    const int x0 = int(x), x1 = std::min(x0 + 1, map.cols - 1);
    const auto* top = map.gains.data() + size_t(y0) * map.cols;
    const auto* bottom = map.gains.data() + size_t(y1) * map.cols;
    const double a = top[x0] + (top[x1] - top[x0]) * (x - x0);
    const double b = bottom[x0] + (bottom[x1] - bottom[x0]) * (x - x0);
    return float(a + (b - a) * (y - y0));
}

constexpr std::array<Point, 4> diagonals{{{-1, -1}, {-1, 1}, {1, -1}, {1, 1}}};
constexpr std::array<Point, 4> axial{{{-1, 0}, {0, -1}, {0, 1}, {1, 0}}};

bool interpolateMissing(const std::vector<Site>& sites, const Layout& layout,
                        int row, int col, const std::array<Point, 4>& offsets,
                        int step, Site& result, const std::vector<uint32_t>& sources,
                        uint32_t& source) {
    double sum = 0, total = 0, boundSum = 0;
    int bounds = 0;
    uint16_t measuredFlags = 0, boundFlags = 0;
    uint32_t measuredSource = 0, boundSource = 0;
    auto provenance = [&](size_t index, uint16_t& flags, uint32_t& combined) {
        flags |= sites[index].flags;
        if (sources.empty() || !sources[index]) return;
        const auto donor = sources[index];
        combined = !combined || combined == donor ? donor : mixedSource;
    };
    for (auto p : offsets) {
        const int y = row + p.row * step, x = col + p.col * step;
        if (y < 0 || x < 0 || y >= layout.height || x >= layout.width) continue;
        const auto index = layout.index(y, x);
        const Site& donor = sites[index];
        if ((donor.flags & uncovered) || ((donor.flags & missing) && !(donor.flags & filled))) continue;
        if (donor.flags & lowerBound) {
            boundSum += donor.signal;
            ++bounds;
            provenance(index, boundFlags, boundSource);
        } else {
            const double weight = 1 / std::max(double(donor.variance), 1e-30);
            sum += donor.signal * weight;
            total += weight;
            provenance(index, measuredFlags, measuredSource);
        }
    }
    if (total) {
        result.signal = float(sum / total);
        result.variance = float(1 / total);
        result.flags |= measuredFlags;
        source = measuredSource;
    } else if (bounds) {
        result.signal = float(boundSum / bounds);
        result.flags |= lowerBound | boundFlags;
        source = boundSource;
    } else return false;
    result.flags |= filled;
    return true;
}

// Mirror about the first and last sites, repeatedly, so any offset lands
// inside the image with its Bayer parity kept (the period 2·(limit − 1) is even).
int mirror(int p, int limit) {
    if (limit < 2) return 0;
    const int period = 2 * (limit - 1);
    p %= period;
    if (p < 0) p += period;
    return p < limit ? p : period - p;
}

struct Tap { size_t index; float weight; };
struct Kernel {
    std::array<Tap, 16> taps{};
    int size = 0;

    void add(size_t index, float weight) {
        for (int i = 0; i < size; ++i) {
            if (taps[i].index == index) {
                taps[i].weight += weight;
                return;
            }
        }
        taps[size++] = {index, weight};
    }

    // Drop taps without a value: uncovered sites and holes repair could not fill.
    void omitUnmeasured(const std::vector<Site>& sites) {
        float total = 0;
        int kept = 0;
        for (int i = 0; i < size; ++i) {
            const auto flags = sites[taps[i].index].flags;
            if ((flags & uncovered) || (flags & (missing | filled)) == missing) continue;
            total += taps[i].weight;
            taps[kept++] = taps[i];
        }
        size = kept;
        for (int i = 0; i < size; ++i) taps[i].weight /= total;
    }

    double evaluate(const std::vector<Site>& sites) const {
        double value = 0;
        for (int i = 0; i < size; ++i) value += sites[taps[i].index].signal * taps[i].weight;
        return value;
    }
};

Kernel bilinear(const Layout& layout, int row, int col, int target) {
    Kernel kernel;
    const int center = layout.colour(row, col);
    auto tap = [&](int y, int x, float weight) {
        kernel.add(layout.index(mirror(row + y, layout.height), mirror(col + x, layout.width)), weight);
    };
    if (center == target) tap(0, 0, 1);
    else if (center == 1) {
        if (layout.colour(row, col ^ 1) == target) {
            tap(0, -1, 0.5f);
            tap(0, 1, 0.5f);
        } else {
            tap(-1, 0, 0.5f);
            tap(1, 0, 0.5f);
        }
    } else {
        for (auto p : target == 1 ? axial : diagonals) tap(p.row, p.col, 0.25f);
    }
    return kernel;
}

// Malvar, He and Cutler, ICASSP 2004, Fig. 2. Preserve half-step coefficients.
constexpr int greenAtRedBlue[5][5] = {
    {0, 0, -1, 0, 0}, {0, 0, 2, 0, 0}, {-1, 2, 4, 2, -1}, {0, 0, 2, 0, 0}, {0, 0, -1, 0, 0}
};
constexpr int opposite[5][5] = {
    {0, 0, -3, 0, 0}, {0, 4, 0, 4, 0}, {-3, 0, 12, 0, -3}, {0, 4, 0, 4, 0}, {0, 0, -3, 0, 0}
};
constexpr int horizontal[5][5] = {
    {0, 0, 1, 0, 0}, {0, -2, 0, -2, 0}, {-2, 8, 10, 8, -2}, {0, -2, 0, -2, 0}, {0, 0, 1, 0, 0}
};

Kernel mhc(const Layout& layout, int row, int col, int target) {
    const int center = layout.colour(row, col);
    const bool transpose = center == 1 && layout.colour(row, col ^ 1) != target;
    const auto& weights = center == 1 ? horizontal : target == 1 ? greenAtRedBlue : opposite;
    const float unit = target == 1 ? 0.125f : 0.0625f;
    Kernel kernel;
    for (int y = -2; y <= 2; ++y) {
        for (int x = -2; x <= 2; ++x) {
            const int weight = transpose ? weights[x + 2][y + 2] : weights[y + 2][x + 2];
            if (weight)
                kernel.add(layout.index(mirror(row + y, layout.height), mirror(col + x, layout.width)),
                           unit * weight);
        }
    }
    return kernel;
}

// Steep-edge interpolation. Beside a source far brighter than its surroundings,
// such as the sun's disc, brightness can change tenfold per site; a straight-line
// average of the neighbours is then dominated by the bright side, so the missing
// colours at a site come out too bright for the one it measured. There the
// estimate works in log space: green along the direction of smaller change
// (Hamilton–Adams style, with the site's own colour correcting the curvature to
// second order, limited to the range of its green neighbours), then red and blue
// as that green times the neighbours' geometric-mean colour ratio (constant hue).
// It blends in between 2× and 4× neighbour contrast and needs every site in the
// 7 × 7 window positive and 8 σ above noise, since logs of noisy values are biased.
constexpr double edgeStart = 0.6931471805599453, edgeFull = 1.3862943611198906;  // log 2, log 4
constexpr double edgeSnr2 = 64;  // 8 σ
constexpr int edgeReach = 3;

// The sites one edge estimate read, for its flags and source IDs.
struct EdgeDependencies {
    std::array<size_t, (2 * edgeReach + 1) * (2 * edgeReach + 1)> index{};
    int size = 0;
    void add(size_t i) {
        for (int k = 0; k < size; ++k) if (index[k] == i) return;
        index[size++] = i;
    }
};

struct EdgeEstimator {
    const std::vector<Site>& sites;
    const Layout& layout;

    size_t indexAt(int y, int x) const { return layout.index(mirror(y, layout.height), mirror(x, layout.width)); }
    const Site& at(int y, int x) const { return sites[indexAt(y, x)]; }
    // Lower bounds (zero variance) pass the noise test and are used at their
    // value, as the linear kernels use them; excluding them falls back to those
    // kernels, which fringe worse beside a clipped core. Their LOWER_BOUND flag
    // reaches the output through the estimate's dependencies.
    bool usable(const Site& s) const {
        return !(s.flags & uncovered) && (s.flags & (missing | filled)) != missing && s.signal > 0 &&
               double(s.signal) * s.signal >= edgeSnr2 * s.variance;
    }
    double log(int y, int x, EdgeDependencies& used) const {
        const size_t i = indexAt(y, x);
        used.add(i);
        return std::log(double(sites[i].signal));
    }

    // Log-range of the taps a straight average would use, or a negative value
    // when any of them is too dark or noisy for log-space estimates.
    double contrast(const std::array<Kernel, 3>& kernels, int own) const {
        double range = 0;
        for (int c = 0; c < 3; ++c) {
            if (c == own) continue;
            double lo = INFINITY, hi = -INFINITY;
            for (int i = 0; i < kernels[c].size; ++i) {
                const auto& s = sites[kernels[c].taps[i].index];
                if (!usable(s)) return -1;
                lo = std::min(lo, std::log(double(s.signal)));
                hi = std::max(hi, std::log(double(s.signal)));
            }
            if (kernels[c].size) range = std::max(range, hi - lo);
        }
        return range;
    }

    bool windowUsable(int y, int x) const {
        for (int dy = -edgeReach; dy <= edgeReach; ++dy)
            for (int dx = -edgeReach; dx <= edgeReach; ++dx)
                if (!usable(at(y + dy, x + dx))) return false;
        return true;
    }

    double logGreen(int y, int x, EdgeDependencies& used) const {
        const double centre = log(y, x, used);
        if (layout.colour(mirror(y, layout.height), mirror(x, layout.width)) == 1) return centre;
        double sum = 0, weights = 0, lo = INFINITY, hi = -INFINITY;
        for (auto [a, b] : {std::pair{0, 1}, std::pair{1, 0}}) {
            const double g0 = log(y - a, x - b, used), g1 = log(y + a, x + b, used);
            const double c0 = log(y - 2 * a, x - 2 * b, used), c1 = log(y + 2 * a, x + 2 * b, used);
            const double estimate = 0.5 * (g0 + g1) + 0.125 * (2 * centre - c0 - c1);
            const double gradient = std::abs(g0 - g1) + std::abs(2 * centre - c0 - c1) + 1e-3;
            sum += estimate / gradient;
            weights += 1 / gradient;
            lo = std::min({lo, g0, g1});
            hi = std::max({hi, g0, g1});
        }
        return std::clamp(sum / weights, lo, hi);
    }

    double estimate(int y, int x, int target, EdgeDependencies& used) const {
        const double green = logGreen(y, x, used);
        if (target == 1) return std::exp(green);
        const int own = layout.colour(y, x);
        double ratio = 0;
        int count = 0;
        auto add = [&](int dy, int dx) {
            ratio += log(y + dy, x + dx, used) - logGreen(y + dy, x + dx, used);
            ++count;
        };
        if (own == 1) {
            if (layout.colour(y, x ^ 1) == target) { add(0, -1); add(0, 1); }
            else { add(-1, 0); add(1, 0); }
        } else {
            for (auto p : diagonals) add(p.row, p.col);
        }
        return std::exp(green + ratio / count);
    }
};

} // namespace

FlatField flatField(std::vector<Site>& sites, const RawFrame& ref, const Options& o) {
    FlatField result;
    result.source = o.vignette;
    const auto& layout = ref.layout;
    int channels = 0;
    std::vector<float> image;
    const bool embedded = o.vignette == "dng";
    if (!embedded && o.vignette != "none") {
        image = readImage(o.vignette, layout, channels);
        if (channels != 1 && channels != 3)
            throw std::runtime_error("flat-field EXR must have one or three channels");
    }
    if (embedded) {
        result.maps = ref.metadata.maps.size();
        if (!result.maps) result.source = "absent";
    }
    if (!result.maps && image.empty()) return result;

    std::vector<RGB> rowMax(layout.height);
    parallelRows(layout.height, o.threads, [&](int y) {
        std::vector<float> gains(layout.width, 1);
        if (embedded) {
            for (const auto& map : ref.metadata.maps) {
                if (y < map.area.top || y >= map.area.bottom || (y - map.area.top) % map.rowPitch) continue;
                for (int x = map.area.left; x < std::min(map.area.right, layout.width);) {
                    gains[x] *= interpolateGain(map, y, x, layout);
                    if (layout.width - x <= map.colPitch) break;
                    x += map.colPitch;
                }
            }
        }
        for (int x = 0; x < layout.width; ++x) {
            const auto i = layout.index(y, x);
            const int colour = layout.colour(y, x);
            const float gain = embedded ? gains[x] : image[i * channels + (channels == 1 ? 0 : colour)];
            if (!(gain > 0) || !std::isfinite(gain)) throw std::runtime_error("flat-field gain must be positive and finite");
            rowMax[y][colour] = std::max(rowMax[y][colour], double(gain));
            sites[i].signal *= gain;
            sites[i].variance *= gain * gain;
        }
    });
    result.maxGain = {};
    for (const auto& row : rowMax)
        for (int c = 0; c < 3; ++c) result.maxGain[c] = std::max(result.maxGain[c], row[c]);
    return result;
}

namespace {

// Opposed-channel highlight reconstruction in white-balanced camera RGB.
// Clipped quads form 8-connected regions. Unclipped quads within two quads of
// a region give its colour as ratios of channel sums, so brighter rim quads
// count more than the steep, aliased edge. A clipped site is estimated from
// the unclipped other-colour sites around it times those ratios, and never
// falls below its own bound or a neighbouring bound carried through the ratios.
// Without a rim the colour is unknown and the estimate is neutral.
void reconstructHighlights(std::vector<Site>& sites, const Layout& layout, const RGB& wb,
                           std::vector<uint32_t>& sources) {
    std::vector<size_t> clipped;
    for (size_t i = 0; i < sites.size(); ++i)
        if (sites[i].flags & lowerBound) clipped.push_back(i);
    if (clipped.empty()) return;

    const int quadWidth = (layout.width + 1) / 2, quadHeight = (layout.height + 1) / 2;
    auto quadOf = [&](size_t i) {
        return size_t(i / layout.width / 2) * quadWidth + (i % layout.width) / 2;
    };
    auto valued = [&](size_t i) {
        const auto flags = sites[i].flags;
        return !(flags & uncovered) && (flags & (missing | filled)) != missing;
    };
    constexpr int32_t unclipped = -1, unlabelled = -2;
    std::vector<int32_t> region(size_t(quadWidth) * quadHeight, unclipped);
    for (size_t i : clipped) region[quadOf(i)] = unlabelled;

    std::vector<std::array<double, 3>> sums;
    std::vector<int> rimQuads;
    // The rim sets the region's colour, so its flags and sources are donors too.
    std::vector<uint16_t> rimFlags;
    std::vector<uint32_t> rimSources;
    auto combine = [](uint32_t current, uint32_t donor) {
        return !donor ? current : !current || current == donor ? donor : mixedSource;
    };
    std::vector<int32_t> seen(region.size(), unclipped);
    std::vector<size_t> queue;
    for (size_t i : clipped) {
        const size_t start = quadOf(i);
        if (region[start] != unlabelled) continue;
        const auto label = int32_t(sums.size());
        sums.push_back({});
        rimQuads.push_back(0);
        rimFlags.push_back(0);
        rimSources.push_back(0);
        region[start] = label;
        queue.assign(1, start);
        for (size_t k = 0; k < queue.size(); ++k) {
            const int qy = int(queue[k] / quadWidth), qx = int(queue[k] % quadWidth);
            for (int dy = -2; dy <= 2; ++dy)
                for (int dx = -2; dx <= 2; ++dx) {
                    const int ny = qy + dy, nx = qx + dx;
                    if (ny < 0 || nx < 0 || ny >= quadHeight || nx >= quadWidth) continue;
                    const size_t n = size_t(ny) * quadWidth + nx;
                    if (region[n] == unlabelled && std::abs(dy) <= 1 && std::abs(dx) <= 1) {
                        region[n] = label;
                        queue.push_back(n);
                    }
                    if (region[n] != unclipped || seen[n] == label) continue;
                    seen[n] = label;
                    if (2 * ny + 1 >= layout.height || 2 * nx + 1 >= layout.width) continue;
                    std::array<double, 3> quad{};
                    bool measured = true;
                    uint16_t quadFlags = 0;
                    uint32_t quadSource = 0;
                    for (int y = 2 * ny; y < 2 * ny + 2 && measured; ++y)
                        for (int x = 2 * nx; x < 2 * nx + 2; ++x) {
                            const size_t j = layout.index(y, x);
                            if (!valued(j) || sites[j].signal <= 0) {
                                measured = false;
                                break;
                            }
                            const int c = layout.colour(y, x);
                            quad[c] += sites[j].signal * wb[c] * (c == 1 ? 0.5 : 1);
                            quadFlags |= sites[j].flags;
                            if (!sources.empty()) quadSource = combine(quadSource, sources[j]);
                        }
                    if (!measured) continue;
                    for (int c = 0; c < 3; ++c) sums[label][c] += quad[c];
                    ++rimQuads[label];
                    rimFlags[label] |= quadFlags;
                    rimSources[label] = combine(rimSources[label], quadSource);
                }
        }
    }

    struct Update { size_t index; float signal, variance; uint16_t flags; uint32_t source; };
    std::vector<Update> updates;
    updates.reserve(clipped.size());
    for (size_t i : clipped) {
        const int y = int(i / layout.width), x = int(i % layout.width);
        const int c = layout.colour(y, x);
        const auto label = region[quadOf(i)];
        const auto& sum = sums[label];
        const bool coloured = rimQuads[label] >= 4 && sum[0] > 0 && sum[1] > 0 && sum[2] > 0;
        double lower = sites[i].signal * wb[c], estimate = 0, estimateVariance = 0;
        int estimates = 0;
        uint16_t donorFlags = 0, lowerFlags = 0;
        uint32_t source = sources.empty() ? 0 : sources[i], lowerSource = 0;
        bool raisedByNeighbour = false;
        for (int k = 0; k < 3; ++k) {
            if (k == c) continue;
            const double ratio = coloured ? sum[c] / sum[k] : 1;
            double total = 0, variance = 0, bound = 0;
            size_t boundSite = 0;
            int count = 0;
            uint16_t flags = 0;
            uint32_t donor = 0;
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    const int ny = y + dy, nx = x + dx;
                    if (ny < 0 || nx < 0 || ny >= layout.height || nx >= layout.width) continue;
                    if (layout.colour(ny, nx) != k) continue;
                    const size_t j = layout.index(ny, nx);
                    if (!valued(j)) continue;
                    const double value = sites[j].signal * wb[k];
                    if (sites[j].flags & lowerBound) {
                        if (value > bound) {
                            bound = value;
                            boundSite = j;
                        }
                        continue;
                    }
                    total += value;
                    variance += sites[j].variance * wb[k] * wb[k];
                    flags |= sites[j].flags;
                    if (!sources.empty() && sources[j])
                        donor = !donor || donor == sources[j] ? sources[j] : mixedSource;
                    ++count;
                }
            if (ratio * bound > lower) {
                lower = ratio * bound;
                lowerFlags = sites[boundSite].flags;
                lowerSource = sources.empty() ? 0 : sources[boundSite];
                raisedByNeighbour = true;
            }
            if (!count) continue;
            const double value = ratio * total / count;
            // Two channels estimate the same quantity when the colour is known.
            // Without it, keep the neutral behaviour: the brightest one wins.
            if (coloured || value > estimate) {
                estimate = coloured ? estimate + value : value;
                estimateVariance = coloured ? estimateVariance + ratio * ratio * variance / (count * count)
                                            : ratio * ratio * variance / (count * count);
                estimates = coloured ? estimates + 1 : 1;
            }
            donorFlags |= flags;
            if (donor) source = !source || source == donor ? donor : mixedSource;
        }
        if (estimates) {
            estimate /= estimates;
            estimateVariance /= double(estimates) * estimates;
        }
        double value = lower, variance = 0;
        uint16_t flags = sites[i].flags | reconstructed;
        if (estimates && estimate > lower) {
            value = estimate;
            variance = estimateVariance;
            flags |= donorFlags;
        } else {
            // The level comes from a bound: this site's own, or a neighbour's
            // carried through the ratio.
            source = sources.empty() ? 0 : sources[i];
            if (raisedByNeighbour) {
                flags |= lowerFlags;
                source = combine(source, lowerSource);
            }
        }
        if (value <= sites[i].signal * wb[c]) continue;
        if (coloured) {
            flags |= rimFlags[label];
            source = combine(source, rimSources[label]);
        }
        updates.push_back({i, float(value / wb[c]), float(variance / (wb[c] * wb[c])), flags, source});
    }
    for (const auto& u : updates) {
        sites[u.index].signal = u.signal;
        sites[u.index].variance = u.variance;
        sites[u.index].flags = u.flags;
        if (!sources.empty()) sources[u.index] = u.source;
    }
}

} // namespace

void repairSites(std::vector<Site>& sites, const Layout& layout, const RGB& wb, const Options& o,
                 std::vector<uint32_t>& sources) {
    std::vector<size_t> holes;
    for (size_t i = 0; i < sites.size(); ++i)
        if ((sites[i].flags & (missing | uncovered)) == missing) holes.push_back(i);
    auto repair = [&](size_t i, Site& value, uint32_t& source) {
        const int y = int(i / layout.width), x = int(i % layout.width);
        const bool green = layout.colour(y, x) == 1;
        return interpolateMissing(sites, layout, y, x, green ? diagonals : axial, green ? 1 : 2, value, sources, source) ||
               interpolateMissing(sites, layout, y, x, green ? axial : diagonals, 2, value, sources, source);
    };

    // Use only original measurements in the first pass, then propagate into wider defects.
    std::vector<size_t> frontier;
    std::vector<Site> repairs(holes.size());
    std::vector<uint32_t> repairSources(sources.empty() ? 0 : holes.size());
    for (size_t k = 0; k < holes.size(); ++k) {
        const auto index = holes[k];
        repairs[k] = sites[index];
        uint32_t source = 0;
        if (repair(index, repairs[k], source)) frontier.push_back(index);
        if (!sources.empty()) repairSources[k] = source;
    }
    for (size_t k = 0; k < holes.size(); ++k) {
        sites[holes[k]] = repairs[k];
        if (!sources.empty()) sources[holes[k]] = repairSources[k];
    }
    for (size_t k = 0; k < frontier.size(); ++k) {
        const size_t index = frontier[k];
        const int y = int(index / layout.width), x = int(index % layout.width);
        const bool green = layout.colour(y, x) == 1;
        auto propagate = [&](const std::array<Point, 4>& offsets, int step) {
            for (auto p : offsets) {
                const int row = y + p.row * step, col = x + p.col * step;
                if (row < 0 || col < 0 || row >= layout.height || col >= layout.width) continue;
                const size_t i = layout.index(row, col);
                auto& site = sites[i];
                uint32_t source = 0;
                if ((site.flags & (missing | filled | uncovered)) == missing && repair(i, site, source)) {
                    if (!sources.empty()) sources[i] = source;
                    frontier.push_back(i);
                }
            }
        };
        propagate(green ? diagonals : axial, green ? 1 : 2);
        propagate(green ? axial : diagonals, 2);
    }
    if (o.highlights == "bound") return;
    if (o.highlights == "reconstruct") {
        reconstructHighlights(sites, layout, wb, sources);
        return;
    }
    parallelRows((layout.height + 1) / 2, o.threads, [&](int row) {
        const int y = row * 2;
        for (int x = 0; x < layout.width; x += 2) {
            double brightest = 0;
            bool clipped = false;
            uint16_t brightestFlags = 0;
            uint32_t brightestSource = 0;
            for (int h = y; h < std::min(y + 2, layout.height); ++h)
                for (int w = x; w < std::min(x + 2, layout.width); ++w) {
                    const auto i = layout.index(h, w);
                    const auto& s = sites[i];
                    clipped |= (s.flags & lowerBound) != 0;
                    const double value = s.signal * wb[layout.colour(h, w)];
                    if (value > brightest) {
                        brightest = value;
                        brightestFlags = s.flags;
                        if (!sources.empty()) brightestSource = sources[i];
                    }
                }
            if (!clipped) continue;
            for (int h = y; h < std::min(y + 2, layout.height); ++h)
                for (int w = x; w < std::min(x + 2, layout.width); ++w) {
                    const auto i = layout.index(h, w);
                    auto& s = sites[i];
                    const double raised = brightest / wb[layout.colour(h, w)];
                    if ((s.flags & lowerBound) && raised > s.signal) {
                        s.signal = float(raised);
                        s.flags |= brightestFlags;
                        if (!sources.empty() && brightestSource)
                            sources[i] = !sources[i] || sources[i] == brightestSource ? brightestSource : mixedSource;
                    }
                }
        }
    });
}

Output render(const std::vector<Site>& sites, const Layout& layout, const Colour& colour,
              double scale, const Options& o, const std::vector<uint32_t>& sources) {
    Output output;
    output.rgb.resize(layout.size() * 3);
    if (o.aux) {
        output.variance.resize(layout.size() * 3);
        output.frames.resize(layout.size());
        output.flags.resize(layout.size());
        if (!sources.empty()) output.sources.resize(layout.size());
    }
    Matrix transform = identity;
    if (o.space == "aces") {
        transform = *colour.matrix;
        for (auto& row : transform)
            for (int c = 0; c < 3; ++c) row[c] *= colour.wb[c];
    }

    const bool useMhc = o.demosaic == "mhc", useEdges = o.edges == "ratio";
    const EdgeEstimator edges{sites, layout};
    std::atomic<size_t> floors{0}, edgePixels{0};
    parallelRows(layout.height, o.threads, [&](int y) {
        size_t flooredRow = 0, edgeRow = 0;
        for (int x = 0; x < layout.width; ++x) {
            const auto i = layout.index(y, x);
            if (sites[i].flags & uncovered) {
                const float floor = float(o.outputFloor);
                std::fill_n(output.rgb.data() + i * 3, 3, floor);
                flooredRow += floor > 0;
                if (o.aux) output.flags[i] = sites[i].flags | (floor > 0 ? floored : 0);
                continue;
            }
            std::array<Kernel, 3> kernels;
            RGB camera{};
            bool sharp = useMhc;
            if (sharp) {
                // MHC's gradient terms need real measurements: no bounds, no
                // uncovered sites and no holes that repair could not fill.
                for (int dy = -2; dy <= 2 && sharp; ++dy)
                    for (int dx = -2; dx <= 2 && sharp; ++dx) {
                        const auto flags = sites[layout.index(mirror(y + dy, layout.height),
                                                              mirror(x + dx, layout.width))].flags;
                        sharp = !(flags & (lowerBound | uncovered)) &&
                                ((flags & (missing | filled)) != missing);
                    }
            }
            for (int c = 0; c < 3; ++c) {
                auto k = bilinear(layout, y, x, c);
                k.omitUnmeasured(sites);
                camera[c] = k.evaluate(sites);
                if (sharp && c != layout.colour(y, x)) {
                    auto corrected = mhc(layout, y, x, c);
                    const double value = corrected.evaluate(sites);
                    if (std::abs(value - camera[c]) <= 0.5 * std::max(camera[c], 0.0)) {
                        k = corrected;
                        camera[c] = value;
                    }
                }
                kernels[c] = k;
            }

            // Steep edges: blend toward the log-space estimate and rescale the
            // kernel to it, so variance keeps its relative size (an approximation:
            // the estimate is nonlinear in sites beyond the kernel). Flags and
            // source IDs also include every site the estimate read.
            bool edge = false;
            std::array<EdgeDependencies, 3> edgeUsed;
            if (useEdges) {
                const int own = layout.colour(y, x);
                const double range = edges.contrast(kernels, own);
                if (range > edgeStart && edges.windowUsable(y, x)) {
                    const double weight = std::min(1.0, (range - edgeStart) / (edgeFull - edgeStart));
                    for (int c = 0; c < 3; ++c) {
                        if (c == own || camera[c] <= 0) continue;
                        const double value = camera[c] + weight * (edges.estimate(y, x, c, edgeUsed[c]) - camera[c]);
                        const float factor = float(value / camera[c]);
                        for (int j = 0; j < kernels[c].size; ++j) kernels[c].taps[j].weight *= factor;
                        camera[c] = value;
                    }
                    edge = true;
                    ++edgeRow;
                }
            }

            bool anyFloor = false;
            uint16_t frames = 0;
            uint32_t flags = 0, source = 0;
            for (int c = 0; c < 3; ++c) {
                double value = 0;
                for (int k = 0; k < 3; ++k) value += transform[c][k] * camera[k];
                value *= scale;
                if (!std::isfinite(float(value))) throw std::runtime_error("output radiance overflows float32");
                if (o.outputFloor > 0 && value < o.outputFloor) {
                    value = o.outputFloor;
                    anyFloor = true;
                }
                output.rgb[i * 3 + c] = float(value);
                if (!o.aux) continue;
                std::array<std::pair<size_t, double>, 48> combined{};
                int count = 0;
                for (int k = 0; k < 3; ++k) {
                    if (transform[c][k] == 0) continue;
                    if (!kernels[k].size) flags |= missing;
                    for (int j = 0; j < kernels[k].size; ++j) {
                        const auto tap = kernels[k].taps[j];
                        if (tap.weight == 0) continue;
                        auto it = std::find_if(combined.begin(), combined.begin() + count,
                                               [&](const auto& entry) { return entry.first == tap.index; });
                        if (it == combined.begin() + count) {
                            combined[count++] = {tap.index, 0};
                        }
                        it->second += transform[c][k] * tap.weight * scale;
                    }
                }
                auto depend = [&](size_t index) {
                    const auto& site = sites[index];
                    frames = std::max(frames, site.frames);
                    flags |= site.flags;
                    if (!sources.empty() && sources[index]) {
                        const uint32_t contributor = sources[index];
                        source = !source || source == contributor ? contributor : mixedSource;
                    }
                };
                double variance = 0;
                for (auto it = combined.begin(); it != combined.begin() + count; ++it) {
                    if (!it->second) continue;
                    variance += it->second * it->second * sites[it->first].variance;
                    depend(it->first);
                }
                for (int k = 0; k < 3; ++k)
                    if (transform[c][k] != 0)
                        for (int j = 0; j < edgeUsed[k].size; ++j) depend(edgeUsed[k].index[j]);
                if (!std::isfinite(float(variance))) throw std::runtime_error("output variance overflows float32");
                output.variance[i * 3 + c] = float(variance);
            }
            flooredRow += anyFloor;
            if (o.aux) {
                output.frames[i] = frames;
                output.flags[i] = flags | (anyFloor ? floored : 0) | (edge ? edgeRatio : 0);
                if (!sources.empty()) output.sources[i] = source;
            }
        }
        floors.fetch_add(flooredRow, std::memory_order_relaxed);
        edgePixels.fetch_add(edgeRow, std::memory_order_relaxed);
    });
    output.flooredPixels = floors;
    output.edgePixels = edgePixels;
    return output;
}
}

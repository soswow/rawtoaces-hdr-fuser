#include "pipeline.h"

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

int mirror(int p, int limit) {
    if (p < 0) p = -p;
    if (p >= limit) p = 2 * (limit - 1) - p;
    return std::abs(p);
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

    const bool useMhc = o.demosaic == "mhc";
    std::atomic<size_t> floors{0};
    parallelRows(layout.height, o.threads, [&](int y) {
        size_t flooredRow = 0;
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
                double variance = 0;
                for (auto it = combined.begin(); it != combined.begin() + count; ++it) {
                    if (!it->second) continue;
                    const auto& site = sites[it->first];
                    variance += it->second * it->second * site.variance;
                    frames = std::max(frames, site.frames);
                    flags |= site.flags;
                    if (!sources.empty() && sources[it->first]) {
                        const uint32_t contributor = sources[it->first];
                        source = !source || source == contributor ? contributor : mixedSource;
                    }
                }
                if (!std::isfinite(float(variance))) throw std::runtime_error("output variance overflows float32");
                output.variance[i * 3 + c] = float(variance);
            }
            flooredRow += anyFloor;
            if (o.aux) {
                output.frames[i] = frames;
                output.flags[i] = flags | (anyFloor ? floored : 0);
                if (!sources.empty()) output.sources[i] = source;
            }
        }
        floors.fetch_add(flooredRow, std::memory_order_relaxed);
    });
    output.flooredPixels = floors;
    return output;
}
}

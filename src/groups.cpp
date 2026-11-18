#include "pipeline.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

namespace rawfuse {
namespace {

using Json = nlohmann::json;

void validateKeys(const Json& object, std::initializer_list<const char*> allowed) {
    if (!object.is_object()) throw std::runtime_error("expected a JSON object");
    for (const auto& entry : object.items())
        if (std::find(allowed.begin(), allowed.end(), entry.key()) == allowed.end())
            throw std::runtime_error("unknown key: " + entry.key());
}

fs::path inputPath(const Json& value, const fs::path& directory) {
    const auto path = fs::path(value.get<std::string>());
    if (path.empty()) throw std::runtime_error("expected a nonempty input path");
    return (directory / path).lexically_normal();
}

Mask intersection(const Mask& a, Mask b) {
    if (a.empty()) return b;
    if (b.empty()) return a;
    for (size_t i = 0; i < a.size(); ++i) b[i] &= a[i];
    return b;
}

// Masks are empty when they select every site.
Mask unite(const Mask& a, const Mask& b) {
    if (a.empty() || b.empty()) return {};
    Mask result = a;
    for (size_t i = 0; i < b.size(); ++i) result[i] |= b[i];
    return result;
}

Mask subtract(const Mask& a, const Mask& b, size_t size) {
    if (b.empty()) return a;
    Mask result = a.empty() ? Mask(size, 1) : a;
    for (size_t i = 0; i < size; ++i) result[i] &= !b[i];
    return result;
}

int quality(const Site& site) {
    if (site.flags & missing) return 0;
    if (site.flags & lowerBound) return 1;
    if (site.flags & adjacentClipFallback) return 2;
    return 3;
}

// Offset matching works on cells of matchCell² sites, per colour.
constexpr int matchCell = 32;
constexpr size_t matchMinimum = 16;
constexpr double matchRelativeVariance = 0.01;  // SNR of at least 10

bool precise(const Site& site) {
    return quality(site) == 3 && site.signal > 0 &&
           site.variance < matchRelativeVariance * double(site.signal) * site.signal;
}

struct OffsetField {
    int rows = 0, cols = 0;
    std::array<std::vector<double>, 3> value;
    std::array<std::vector<uint8_t>, 3> reached;

    // Bilinear between cell centres; false where no cell nearby was reached.
    bool at(int y, int x, int colour, double& offset) const {
        const double fy = std::clamp((y + 0.5) / matchCell - 0.5, 0.0, double(rows - 1));
        const double fx = std::clamp((x + 0.5) / matchCell - 0.5, 0.0, double(cols - 1));
        const int r0 = int(fy), c0 = int(fx);
        const int r1 = std::min(r0 + 1, rows - 1), c1 = std::min(c0 + 1, cols - 1);
        const double wy = fy - r0, wx = fx - c0;
        const int rs[] = {r0, r0, r1, r1}, cs[] = {c0, c1, c0, c1};
        const double ws[] = {(1 - wy) * (1 - wx), (1 - wy) * wx, wy * (1 - wx), wy * wx};
        double sum = 0, weight = 0;
        for (int k = 0; k < 4; ++k) {
            const size_t i = size_t(rs[k]) * cols + cs[k];
            if (!reached[colour][i]) continue;
            sum += ws[k] * value[colour][i];
            weight += ws[k];
        }
        if (weight <= 1e-9) return false;
        offset = sum / weight;
        return true;
    }
};

} // namespace

// Membrane (harmonic) interpolation: unknown cells take the average of their
// neighbours, with known cells fixed. Cells in a component without a known cell
// stay unreached.
void membrane(std::vector<double>& value, const std::vector<uint8_t>& known, int rows, int cols,
              std::vector<uint8_t>& reached) {
    reached = known;
    std::vector<size_t> queue;
    for (size_t i = 0; i < known.size(); ++i) if (known[i]) queue.push_back(i);
    // Nearest-known start, then Gauss-Seidel with over-relaxation.
    for (size_t head = 0; head < queue.size(); ++head) {
        const int r = int(queue[head] / cols), c = int(queue[head] % cols);
        const int dr[] = {-1, 1, 0, 0}, dc[] = {0, 0, -1, 1};
        for (int k = 0; k < 4; ++k) {
            const int rr = r + dr[k], cc = c + dc[k];
            if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) continue;
            const size_t j = size_t(rr) * cols + cc;
            if (reached[j]) continue;
            reached[j] = 1;
            value[j] = value[queue[head]];
            queue.push_back(j);
        }
    }
    double scale = 0;
    for (size_t i = 0; i < known.size(); ++i) if (known[i]) scale = std::max(scale, std::abs(value[i]));
    const double tolerance = 1e-6 * std::max(scale, 1e-12);
    for (int iteration = 0; iteration < 5000; ++iteration) {
        double change = 0;
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
                const size_t i = size_t(r) * cols + c;
                if (known[i] || !reached[i]) continue;
                double sum = 0;
                int n = 0;
                if (r > 0) { sum += value[i - cols]; ++n; }
                if (r + 1 < rows) { sum += value[i + cols]; ++n; }
                if (c > 0) { sum += value[i - 1]; ++n; }
                if (c + 1 < cols) { sum += value[i + 1]; ++n; }
                const double next = value[i] + 1.8 * (sum / n - value[i]);
                change = std::max(change, std::abs(next - value[i]));
                value[i] = next;
            }
        if (change < tolerance) break;
    }
}

std::vector<InputGroup> readGroups(const fs::path& path) {
    try {
        std::ifstream file(path);
        if (!file) throw std::runtime_error("could not open manifest");
        const auto root = Json::parse(file);
        validateKeys(root, {"groups"});
        const auto& entries = root.at("groups");
        if (!entries.is_array() || entries.empty() || entries.size() >= mixedSource)
            throw std::runtime_error("groups must be a nonempty array of up to 65534 input groups");

        const auto directory = fs::absolute(path).parent_path();
        std::vector<InputGroup> groups;
        std::unordered_map<std::string, size_t> indices;
        size_t totalFrames = 0;
        for (const auto& entry : entries) {
            validateKeys(entry, {"name", "frames", "valid_mask", "contribution_mask",
                                 "calibration_mask", "drift_mask", "calibration_reference", "ev",
                                 "match"});
            InputGroup group;
            group.name = entry.at("name").get<std::string>();
            if (group.name.empty() || indices.count(group.name))
                throw std::runtime_error("group names must be nonempty and unique");
            const auto& frames = entry.at("frames");
            if (!frames.is_array() || frames.empty())
                throw std::runtime_error(group.name + ": frames must be a nonempty array of paths");
            totalFrames += frames.size();
            if (totalFrames > 65535) throw std::runtime_error("too many frames");
            for (const auto& frame : frames) group.files.push_back(inputPath(frame, directory));
            if (entry.contains("valid_mask")) group.validMask = inputPath(entry.at("valid_mask"), directory);
            if (entry.contains("contribution_mask"))
                group.contributionMask = inputPath(entry.at("contribution_mask"), directory);
            if (entry.contains("calibration_mask"))
                group.calibrationMask = inputPath(entry.at("calibration_mask"), directory);
            if (entry.contains("drift_mask")) group.driftMask = inputPath(entry.at("drift_mask"), directory);

            if (entry.contains("ev")) {
                const auto& ev = entry.at("ev");
                if (ev == "auto") group.autoEV = true;
                else if (ev.is_array() && ev.size() == 3) group.ev = ev.get<RGB>();
                else if (ev.is_number()) {
                    const double neutral = ev.get<double>();
                    group.ev = RGB{neutral, neutral, neutral};
                } else throw std::runtime_error(group.name + ": ev must be auto, a number, or three RGB values");
            }
            if (entry.contains("match")) {
                const auto& match = entry.at("match");
                if (match != "offset" && match != "none")
                    throw std::runtime_error(group.name + ": match must be offset or none");
                group.matchOffset = match == "offset";
            }
            if (group.ev)
                for (const double ev : *group.ev)
                    if (!std::isfinite(ev) || ev < 0 || ev > 64)
                        throw std::runtime_error(group.name + ": ev must be between 0 and 64");
            const bool hasCalibration = !group.calibrationMask.empty() || !group.driftMask.empty() ||
                                        entry.contains("calibration_reference");
            if (hasCalibration && !group.ev && !group.autoEV)
                throw std::runtime_error(group.name + ": calibration settings require ev");
            if (group.autoEV || hasCalibration) {
                group.calibrationReference = entry.value("calibration_reference", groups.empty() ?
                                                         group.name : groups.front().name);
                if (group.calibrationReference.empty())
                    throw std::runtime_error(group.name + ": calibration_reference must be nonempty");
            }
            indices.emplace(group.name, groups.size());
            groups.push_back(std::move(group));
        }
        for (const auto& group : groups) {
            if (group.calibrationReference.empty()) continue;
            const auto ref = indices.find(group.calibrationReference);
            if (ref == indices.end() || groups[ref->second].autoEV || groups[ref->second].ev)
                throw std::runtime_error(group.name + ": calibration_reference must name an unattenuated group");
        }
        return groups;
    } catch (const std::exception& error) {
        throw std::runtime_error(path.string() + ": " + error.what());
    }
}

GroupMerge mergeGroups(const FrameGroups& frames, const Sidecars& sidecars,
                       const Options& o, Warnings& warnings) {
    const auto& validity = sidecars.valid;
    const auto layout = frames.front()->front().layout;
    // A group's mosaics are needed by its own contribution merge and by any
    // later group that calibrates against it. Release them after that use.
    std::vector<size_t> lastUse(frames.size());
    for (size_t g = 0; g < frames.size(); ++g) {
        lastUse[g] = std::max(lastUse[g], g);
        if (o.groups[g].calibrationReference.empty()) continue;
        for (size_t r = 0; r < frames.size(); ++r)
            if (o.groups[r].name == o.groups[g].calibrationReference) lastUse[r] = std::max(lastUse[r], g);
    }
    // The mosaic has one site per pixel; the full-colour merge (--demosaic
    // splat) three, pixel-major. Masks are per pixel.
    const size_t channels = o.demosaic == "splat" ? 3 : 1;
    auto colourOf = [&](int y, int x, size_t channel) {
        return channels == 1 ? layout.colour(y, x) : int(channel);
    };
    GroupMerge result;
    result.groups.resize(frames.size());
    result.sources.resize(layout.size() * channels);
    result.mosaic.sites.resize(layout.size() * channels, {0, 0, uint16_t(missing | uncovered), 0});
    Mask shadows(layout.size() * channels);
    Options automaticReference = o;
    automaticReference.reference = 0;
    for (size_t g = 0; g < frames.size(); ++g) {
        const auto& group = o.groups[g];
        auto& groupResult = result.groups[g];
        const bool attenuated = group.autoEV || group.ev.has_value();
        if (group.ev) groupResult.ev = *group.ev;
        if (!group.calibrationReference.empty()) {
            const size_t ref = size_t(std::find_if(o.groups.begin(), o.groups.end(), [&](const InputGroup& a) {
                return a.name == group.calibrationReference;
            }) - o.groups.begin());
            const auto& drift = sidecars.drift[g];
            const auto mask = subtract(intersection(validity[g],
                intersection(validity[ref], sidecars.calibration[g])), drift, layout.size());
            const auto driftSites = drift.empty() ? Mask() : intersection(validity[ref], drift);
            const auto selection = driftSites.empty() ? mask : unite(mask, driftSites);
            auto reference = mergePlain(*frames[ref], ref ? automaticReference : o, selection);
            try {
                groupResult.calibration = calibrateND(reference, *frames[ref], *frames[g], mask,
                                                     driftSites, group.ev, o, warnings);
            } catch (const std::exception& error) {
                throw std::runtime_error(group.name + ": " + error.what());
            }
            groupResult.ev = groupResult.calibration->ev;
        }

        // The group is merged wherever it is valid. Inside its contribution mask it
        // replaces equal-quality samples; everywhere else it only fills sites that
        // earlier groups measured worse.
        const auto& contributionMask = sidecars.contribution[g];
        const auto& selection = validity[g];
        const auto& options = g ? automaticReference : o;
        auto contribution = channels == 1 ? mergePlain(*frames[g], options, selection, true)
                                          : mergeSplat(*frames[g], options, selection, true);
        if (!g) result.mosaic.reference = contribution.reference;
        RGB gain = unitWB;
        for (int c = 0; c < 3; ++c) gain[c] = std::exp2(groupResult.ev[c]);
        auto replacesAt = [&](size_t i) { return !g || (!contributionMask.empty() && contributionMask[i]); };

        // An additive offset between this group and the composite so far (the open
        // bracket's glare veil against a blocked bracket, for example), measured
        // outside the replacement area where both are precise.
        std::optional<OffsetField> field;
        if (g && group.matchOffset) {
            OffsetField f;
            f.rows = (layout.height + matchCell - 1) / matchCell;
            f.cols = (layout.width + matchCell - 1) / matchCell;
            const size_t cells = size_t(f.rows) * f.cols;
            std::array<std::vector<uint8_t>, 3> known;
            for (int c = 0; c < 3; ++c) {
                f.value[c].assign(cells, 0);
                known[c].assign(cells, 0);
            }
            std::vector<size_t> rowSamples(f.rows);
            parallelRows(f.rows, o.threads, [&](int cellRow) {
                std::vector<std::array<std::vector<float>, 3>> samples(f.cols);
                const int end = std::min(layout.height, (cellRow + 1) * matchCell);
                for (int y = cellRow * matchCell; y < end; ++y)
                    for (int x = 0; x < layout.width; ++x) {
                        const size_t i = layout.index(y, x);
                        if ((!selection.empty() && !selection[i]) || replacesAt(i)) continue;
                        for (size_t channel = 0; channel < channels; ++channel) {
                            const size_t j = i * channels + channel;
                            const auto& before = result.mosaic.sites[j];
                            const auto& site = contribution.sites[j];
                            if (!precise(before) || !precise(site)) continue;
                            const int colour = colourOf(y, x, channel);
                            const double factor = attenuated ? gain[colour] : 1;
                            samples[x / matchCell][colour].push_back(float(site.signal * factor - before.signal));
                        }
                    }
                for (int c = 0; c < f.cols; ++c)
                    for (int colour = 0; colour < 3; ++colour) {
                        auto& values = samples[c][colour];
                        rowSamples[cellRow] += values.size();
                        // A full-colour merge has a sample of each colour at every pixel,
                        // interpolated from a quarter (red, blue) or half (green) as many
                        // sites; the minimum counts sites, as in the mosaic.
                        const size_t perSite = channels == 1 ? 1 : colour == 1 ? 2 : 4;
                        if (values.size() < matchMinimum * perSite) continue;
                        const size_t cell = size_t(cellRow) * f.cols + c;
                        auto middle = values.begin() + values.size() / 2;
                        std::nth_element(values.begin(), middle, values.end());
                        f.value[colour][cell] = *middle;
                        known[colour][cell] = 1;
                    }
            });
            OffsetMatch match;
            for (const size_t n : rowSamples) match.samples += n;
            for (int c = 0; c < 3; ++c) {
                for (const auto k : known[c]) match.knownCells += k;
                match.anchored[c] = std::find(known[c].begin(), known[c].end(), 1) != known[c].end();
            }
            if (std::all_of(match.anchored.begin(), match.anchored.end(), [](bool a) { return a; })) {
                for (int c = 0; c < 3; ++c) membrane(f.value[c], known[c], f.rows, f.cols, f.reached[c]);
                field = std::move(f);
            }
            groupResult.match = match;
        }

        struct RowMatch {
            size_t corrected = 0, unresolved = 0, capped = 0;
            std::array<std::vector<float>, 3> offset, share;
        };
        // Inserted sites are counted whenever matching is on, also when no field
        // could be built, so the report shows what stayed uncorrected.
        const bool matching = groupResult.match.has_value();
        std::vector<RowMatch> rowMatch(matching ? layout.height : 0);
        parallelRows(layout.height, o.threads, [&](int y) {
            for (int x = 0; x < layout.width; ++x)
            for (size_t channel = 0; channel < channels; ++channel) {
                const size_t i = layout.index(y, x), j = i * channels + channel;
                if (!selection.empty() && !selection[i]) continue;
                auto& current = result.mosaic.sites[j];
                auto site = contribution.sites[j];
                const bool replaces = replacesAt(i);
                const int rank = quality(site), held = quality(current);
                if (rank < held || (rank == held && rank > 1 && !replaces)) continue;
                if (rank == 0) {
                    current = site;
                    continue;
                }
                const int colour = colourOf(y, x, channel);
                if (attenuated) {
                    const double factor = gain[colour];
                    site.signal = float(site.signal * factor);
                    site.variance = float(site.variance * factor * factor);
                    if (!std::isfinite(site.signal) || !std::isfinite(site.variance))
                        throw std::runtime_error(group.name + ": attenuation correction overflows float32");
                    site.flags |= fromND;
                }
                double offset = 0;
                const bool resolved = field && field->at(y, x, colour, offset);
                // The veil is part of what the group measured, so an interpolated offset
                // larger than the measurement is wrong there; remove at most all of it.
                bool capped = false;
                if (offset > site.signal && site.signal > 0) {
                    offset = site.signal;
                    capped = true;
                }
                const float uncorrected = site.signal;
                site.signal = float(site.signal - offset);
                // Between two clipping bounds, keep the more informative (higher) one.
                if (rank == 1 && held == 1 && site.signal <= current.signal) continue;
                // A measurement well below an earlier clipping bound contradicts it: the
                // scene, glare, masks or attenuation differ. A contribution mask asks for
                // this group regardless; a filling group keeps the bound.
                if (rank > 1 && held == 1 && site.signal < 0.9f * current.signal) {
                    if (!replaces) {
                        current.flags |= inconsistent;
                        continue;
                    }
                    site.flags |= inconsistent;
                }
                current = site;
                result.sources[j] = uint32_t(g + 1);
                shadows[j] = contribution.shadowSites[j];
                if (!matching) continue;
                auto& row = rowMatch[y];
                if (!resolved) {
                    ++row.unresolved;
                    continue;
                }
                ++row.corrected;
                row.capped += capped;
                row.offset[colour].push_back(float(offset));
                if (uncorrected > 0) row.share[colour].push_back(float(offset / uncorrected));
            }
        });
        if (matching) {
            auto& match = *groupResult.match;
            std::array<std::vector<double>, 3> offsets, shares;
            for (auto& row : rowMatch) {
                match.corrected += row.corrected;
                match.unresolved += row.unresolved;
                match.capped += row.capped;
                for (int c = 0; c < 3; ++c) {
                    offsets[c].insert(offsets[c].end(), row.offset[c].begin(), row.offset[c].end());
                    shares[c].insert(shares[c].end(), row.share[c].begin(), row.share[c].end());
                }
            }
            for (int c = 0; c < 3; ++c) {
                if (!offsets[c].empty()) match.median[c] = median(offsets[c]);
                if (!shares[c].empty()) match.share[c] = median(shares[c]);
            }
            match.status = !match.corrected && !match.unresolved ? "no_sites" :
                           !field ? "no_overlap" : match.unresolved ? "partial" : "corrected";
        }
        for (size_t r = 0; r <= g; ++r)
            if (lastUse[r] == g) releaseMosaics(*frames[r]);
    }
    for (size_t i = 0; i < result.sources.size(); ++i) {
        result.mosaic.shadow += shadows[i];
        if (result.sources[i]) {
            auto& group = result.groups[result.sources[i] - 1];
            if (result.mosaic.sites[i].flags & lowerBound) ++group.boundSites;
            else ++group.measuredSites;
        }
    }
    return result;
}
}

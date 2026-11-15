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

int quality(const Site& site) {
    if (site.flags & missing) return 0;
    if (site.flags & lowerBound) return 1;
    if (site.flags & adjacentClipFallback) return 2;
    return 3;
}

} // namespace

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
                                 "calibration_mask", "calibration_reference", "ev"});
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

            if (entry.contains("ev")) {
                const auto& ev = entry.at("ev");
                if (ev == "auto") group.autoEV = true;
                else if (ev.is_array() && ev.size() == 3) group.ev = ev.get<RGB>();
                else if (ev.is_number()) {
                    const double neutral = ev.get<double>();
                    group.ev = RGB{neutral, neutral, neutral};
                } else throw std::runtime_error(group.name + ": ev must be auto, a number, or three RGB values");
            }
            if (group.ev)
                for (const double ev : *group.ev)
                    if (!std::isfinite(ev) || ev < 0 || ev > 64)
                        throw std::runtime_error(group.name + ": ev must be between 0 and 64");
            const bool hasCalibration = !group.calibrationMask.empty() || entry.contains("calibration_reference");
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
    GroupMerge result;
    result.groups.resize(frames.size());
    result.sources.resize(layout.size());
    result.mosaic.sites.resize(layout.size(), {0, 0, uint16_t(missing | uncovered), 0});
    Mask shadows(layout.size());
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
            const auto mask = intersection(validity[g],
                intersection(validity[ref], sidecars.calibration[g]));
            auto reference = mergePlain(*frames[ref], ref ? automaticReference : o, mask);
            try {
                groupResult.calibration = calibrateND(reference, *frames[ref], *frames[g], mask,
                                                     group.ev, false, o, warnings);
            } catch (const std::exception& error) {
                throw std::runtime_error(group.name + ": " + error.what());
            }
            groupResult.ev = groupResult.calibration->ev;
        }

        const auto selection = intersection(validity[g], sidecars.contribution[g]);
        auto contribution = mergePlain(*frames[g], g ? automaticReference : o, selection, true);
        if (!g) result.mosaic.reference = contribution.reference;
        RGB gain = unitWB;
        for (int c = 0; c < 3; ++c) gain[c] = std::exp2(groupResult.ev[c]);
        parallelRows(layout.height, o.threads, [&](int y) {
            for (int x = 0; x < layout.width; ++x) {
                const size_t i = layout.index(y, x);
                if (!selection.empty() && !selection[i]) continue;
                auto& current = result.mosaic.sites[i];
                auto site = contribution.sites[i];
                const int rank = quality(site);
                if (rank < quality(current)) continue;
                if (rank == 0) {
                    current = site;
                    continue;
                }
                if (attenuated) {
                    const double factor = gain[layout.colour(y, x)];
                    site.signal = float(site.signal * factor);
                    site.variance = float(site.variance * factor * factor);
                    if (!std::isfinite(site.signal) || !std::isfinite(site.variance))
                        throw std::runtime_error(group.name + ": attenuation correction overflows float32");
                    site.flags |= fromND;
                }
                // Between two clipping bounds, keep the more informative (higher) one.
                if (rank == 1 && quality(current) == 1 && site.signal <= current.signal) continue;
                current = site;
                result.sources[i] = uint32_t(g + 1);
                shadows[i] = contribution.shadowSites[i];
            }
        });
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

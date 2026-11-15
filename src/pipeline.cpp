#include "pipeline.h"

#include <OpenEXR/ImfThreading.h>
#include <OpenImageIO/imageio.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace rawfuse {
namespace {

using Json = nlohmann::ordered_json;
using Clock = std::chrono::steady_clock;

struct Timings {
    Json values = Json::object();
    Clock::time_point started = Clock::now(), last = started;

    void mark(const char* stage) {
        const auto now = Clock::now();
        values[stage] = std::chrono::duration<double>(now - last).count();
        values["total"] = std::chrono::duration<double>(now - started).count();
        last = now;
    }
};

bool sameFile(const fs::path& a, const fs::path& b) {
    if (fs::weakly_canonical(a) == fs::weakly_canonical(b)) return true;
    std::error_code ec;
    return fs::equivalent(a, b, ec);
}

void validatePaths(const Options& o) {
    std::vector<fs::path> inputs, outputs{o.output, o.report};
    if (!o.writeCfa.empty()) outputs.push_back(o.writeCfa);
    auto addFrames = [&](const std::vector<fs::path>& frames) {
        for (const auto& path : frames) {
            if (!fs::is_regular_file(path)) throw std::runtime_error("input is not a regular file: " + path.string());
            if (std::any_of(inputs.begin(), inputs.end(), [&](const auto& input) { return sameFile(input, path); }))
                throw std::runtime_error("duplicate input frame: " + path.string());
            inputs.push_back(path);
        }
    };
    addFrames(o.plain);
    addFrames(o.nd);
    for (const auto& group : o.groups) addFrames(group.files);
    for (const auto& group : o.groups) {
        for (const auto& mask : {group.validMask, group.contributionMask, group.calibrationMask})
            if (!mask.empty()) inputs.push_back(mask);
    }
    for (const auto& sidecar : {o.darkModel, o.ndMask, o.manifest})
        if (!sidecar.empty()) inputs.push_back(sidecar);
    if (o.vignette != "none" && o.vignette != "dng") inputs.emplace_back(o.vignette);
    if (o.badPixels != "none" && o.badPixels != "dng") inputs.emplace_back(o.badPixels);
    if (o.matrix.rfind("custom:", 0) == 0) inputs.emplace_back(o.matrix.substr(7));

    for (size_t i = 0; i < outputs.size(); ++i) {
        const auto& output = outputs[i];
        for (const auto& input : inputs)
            if (sameFile(output, input)) throw std::runtime_error("output aliases input: " + output.string());
        for (size_t j = 0; j < i; ++j)
            if (sameFile(output, outputs[j])) throw std::runtime_error("output paths must be distinct");
        const auto status = fs::symlink_status(output);
        if (fs::exists(status) && !(o.force && fs::is_regular_file(status))) {
            throw std::runtime_error("output already exists: " + output.string() +
                                     (fs::is_regular_file(status) ? " (use --force)" : ""));
        }
        if (!fs::is_directory(fs::absolute(output).parent_path()))
            throw std::runtime_error("output directory does not exist: " + output.parent_path().string());
    }
}

class StagedFile {
public:
    explicit StagedFile(fs::path target) : target(fs::absolute(std::move(target))) {
        std::random_device random;
        for (int i = 0; i < 32; ++i) {
            directory = this->target.parent_path() / (".rawfuse-" + std::to_string(random()));
            if (fs::create_directory(directory)) return;
        }
        throw std::runtime_error("could not reserve a temporary output directory");
    }

    StagedFile(const StagedFile&) = delete;
    StagedFile& operator=(const StagedFile&) = delete;

    ~StagedFile() {
        std::error_code ignored;
        fs::remove_all(directory, ignored);
    }

    fs::path path() const { return directory / "output"; }

    void publish(bool force) {
        const auto status = fs::symlink_status(target);
        if (force && fs::is_regular_file(status)) {
            fs::create_hard_link(target, directory / "previous");
            fs::rename(path(), target);
            replaced = true;
        } else {
            fs::create_hard_link(path(), target); // Atomically refuses a concurrent writer's file.
        }
        published = true;
    }

    void rollback() {
        if (!published) return;
        if (replaced) fs::rename(directory / "previous", target);
        else if (sameFile(path(), target)) fs::remove(target);
    }

private:
    fs::path target, directory;
    bool published = false, replaced = false;
};

// run() configures library-wide thread pools; restore the caller's settings.
class ThreadSettings {
public:
    explicit ThreadSettings(unsigned threads) : exr(OPENEXR_IMF_NAMESPACE::globalThreadCount()) {
        hasOIIO = OIIO::getattribute("threads", oiio);
        OPENEXR_IMF_NAMESPACE::setGlobalThreadCount(int(threads));
        OIIO::attribute("threads", int(threads));
    }
    ThreadSettings(const ThreadSettings&) = delete;
    ThreadSettings& operator=(const ThreadSettings&) = delete;
    ~ThreadSettings() {
        OPENEXR_IMF_NAMESPACE::setGlobalThreadCount(exr);
        if (hasOIIO) OIIO::attribute("threads", oiio);
    }

private:
    int exr, oiio = 0;
    bool hasOIIO = false;
};

Json rgbReport(const RGB& rgb) { return {{"R", rgb[0]}, {"G", rgb[1]}, {"B", rgb[2]}}; }

std::string noiseSource(const RawFrame& frame, const Options& o) {
    if (o.noise != "dng") return "custom";
    return frame.metadata.noise.empty() ? "fallback" : "dng";
}

Json noiseReport(const RawFrame& f, const Options& o) {
    if (f.noise[0].shot == f.noise[1].shot && f.noise[0].shot == f.noise[2].shot &&
        f.noise[0].read == f.noise[1].read && f.noise[0].read == f.noise[2].read)
        return {{"source", noiseSource(f, o)}, {"S", f.noise[0].shot}, {"O", f.noise[0].read}};
    return {{"source", noiseSource(f, o)}, {"planes", {
        {"R", {{"S", f.noise[0].shot}, {"O", f.noise[0].read}}},
        {"G", {{"S", f.noise[1].shot}, {"O", f.noise[1].read}}},
        {"B", {{"S", f.noise[2].shot}, {"O", f.noise[2].read}}}}}};
}

Json inputReport(const RawFrame& f, const char* role, const Options& o) {
    const auto& m = f.metadata;
    return {{"path", f.path.string()}, {"role", role}, {"index", f.inputIndex + 1},
            {"exposure_time", m.exposure}, {"seconds", m.seconds}, {"iso", m.iso},
            {"f_number", m.aperture}, {"effective_exposure", f.effective},
            {"as_shot_neutral", m.neutral ? Json(*m.neutral) : Json(nullptr)},
            {"image_description", m.description}, {"xmp", m.xmp}, {"capture_date", m.captureDate},
            {"baseline_exposure", m.baseline}, {"noise_profile", noiseReport(f, o)},
            {"hot_sites", {{"kept", f.hotKept}, {"missing", f.hotMissing}, {"model_rejected", f.hotRejected}}}};
}

Json whiteBalanceSpread(const std::vector<RawFrame>& frames, Warnings& warnings) {
    RGB low{}, high{}, spread{};
    low.fill(std::numeric_limits<double>::infinity());
    for (const auto& frame : frames) {
        if (!frame.hasWB) continue;
        for (int i = 0; i < 3; ++i) {
            low[i] = std::min(low[i], frame.wb[i]);
            high[i] = std::max(high[i], frame.wb[i]);
        }
    }
    if (!std::isfinite(low[0])) return nullptr;
    for (int i = 0; i < 3; ++i) spread[i] = std::log2(high[i] / low[i]);
    if (*std::max_element(spread.begin(), spread.end()) > 0.05)
        warnings.push_back("AsShotNeutral varies across plain frames; applying one white balance after fusion");
    return {{"min", rgbReport(low)}, {"max", rgbReport(high)}, {"ev", rgbReport(spread)}};
}

Json darkReport(const DarkFit& fit) {
    auto optional = [](const std::optional<double>& value) { return value ? Json(*value) : Json(nullptr); };
    return {{"method", fit.method}, {"offset_scale", optional(fit.offsetScale)},
            {"rate_scale", optional(fit.rateScale)}, {"relative_scatter", optional(fit.relativeScatter)},
            {"relative_uncertainty", fit.uncertainty}, {"sites", fit.sites}, {"observations", fit.observations}};
}

Json ndReport(const NDResult& nd, size_t count, bool manual, const fs::path& mask, const std::string& use) {
    Json blocks = Json::array(), rejected = Json::array();
    for (const auto& block : nd.blocks) {
        blocks.push_back({{"row", block.row}, {"col", block.col}, {"ev", rgbReport(block.ev)},
                          {"samples", block.samples}});
        if (block.removed) rejected.push_back({{"row", block.row}, {"col", block.col}});
    }
    return {{"frames", count}, {"method", manual ? "manual" : "auto"}, {"ev", rgbReport(nd.ev)},
            {"measured_ev", rgbReport(nd.measured)}, {"mad_ev", rgbReport(nd.mad)},
            {"samples", {{"R", nd.samples[0]}, {"G", nd.samples[1]}, {"B", nd.samples[2]}}},
            {use == "calibration" ? "blocks_flagged" : "blocks_removed", std::move(rejected)},
            {"blocks", std::move(blocks)},
            {"scene_drift_ev", nd.sceneDrift ? Json(*nd.sceneDrift) : Json(nullptr)},
            {"scene_drift_samples", nd.sceneDriftSamples},
            {"blocked_lower_bound_sites", nd.blockedLowerBounds},
            {"mask", mask.string()}, {"use", use}};
}

Json mergeReport(const MergeResult& merge, const Output& output) {
    std::array<size_t, 7> counts{};
    size_t unsupported = 0, adjacent = 0, outside = 0;
    for (const auto& site : merge.sites) {
        for (int bit = 0; bit < 6; ++bit) counts[bit] += (site.flags >> bit) & 1;
        unsupported += !site.frames;
        adjacent += (site.flags & adjacentClipFallback) != 0;
        outside += (site.flags & uncovered) != 0;
    }
    const double size = double(merge.sites.size());
    return {{"lower_bound_fraction", counts[0] / size}, {"from_nd_fraction", counts[1] / size},
            {"moving_fraction", counts[2] / size}, {"filled_fraction", counts[3] / size},
            {"missing_fraction", counts[4] / size}, {"nd_inconsistent_fraction", counts[5] / size},
            {"shadow_fraction", merge.shadow / size}, {"unsupported_fraction", unsupported / size},
            {"adjacent_clip_fallback_fraction", adjacent / size},
            {"uncovered_fraction", outside / size},
            {"floored_fraction", output.flooredPixels / size}};
}

void fuse(const Options& o, Warnings& warnings) {
    Timings time;
    validatePaths(o);
    ThreadSettings threads(o.threads);
    const bool grouped = !o.groups.empty();
    std::vector<std::vector<RawFrame>> storage(grouped ? o.groups.size() : 2);
    FrameGroups frames;
    for (auto& group : storage) frames.push_back(&group);
    const auto& plain = storage.front();
    std::optional<Sidecars> sidecars;
    auto read = [&](const std::vector<fs::path>& files, std::vector<RawFrame>& frames) {
        for (const auto& file : files) {
            if (o.verbose) std::cerr << "Reading " << file.string() << '\n';
            auto frame = readRaw(file, warnings);
            frame.inputIndex = frames.size();
            frames.push_back(std::move(frame));
            // Fail on a bad sidecar before decoding the remaining frames.
            if (!sidecars) sidecars = loadSidecars(o, frames.back().layout);
        }
    };
    if (grouped) {
        for (size_t g = 0; g < o.groups.size(); ++g) read(o.groups[g].files, storage[g]);
    } else {
        read(o.plain, storage.front());
        read(o.nd, storage[1]);
    }
    const double referenceAperture = validateFrames(frames, warnings);
    auto byExposure = [](const auto& a, const auto& b) { return a.effective < b.effective; };
    for (auto* group : frames) std::stable_sort(group->begin(), group->end(), byExposure);
    time.mark("decode");

    const auto darkFits = prepareFrames(frames, *sidecars, o, warnings);
    time.mark("prepare");
    if (o.verbose) std::cerr << "Merging Bayer samples\n";
    MergeResult merge;
    std::optional<NDResult> nd;
    std::vector<uint32_t> sources;
    std::vector<GroupResult> groupResults;
    if (grouped) {
        auto composition = mergeGroups(frames, *sidecars, o, warnings);
        merge = std::move(composition.mosaic);
        sources = std::move(composition.sources);
        groupResults = std::move(composition.groups);
        time.mark("group_merge");
    } else {
        merge = mergePlain(plain, o, {}, !storage[1].empty());
        time.mark("plain_merge");
        if (!storage[1].empty()) {
            if (o.verbose) std::cerr << "Calibrating and merging ND exposures\n";
            nd = mergeND(merge, plain, storage[1], sidecars->ndCoverage, o, warnings);
            time.mark("nd_merge");
        }
    }
    sidecars.reset();
    const auto& ref = plain[merge.reference];
    for (auto* group : frames) releaseMosaics(*group);

    const auto colour = configureColour(plain, merge.reference, o, warnings);
    const auto spread = whiteBalanceSpread(plain, warnings);
    time.mark("colour_setup");
    if (o.vignette == "dng" && std::any_of(frames.begin(), frames.end(), [&](const auto* group) {
        return std::any_of(group->begin(), group->end(), [&](const auto& frame) {
            return frame.metadata.maps != ref.metadata.maps;
        });
    })) warnings.push_back("input frames carry differing GainMaps; applying the reference frame maps");
    const auto field = flatField(merge.sites, ref, o);
    for (const auto& site : merge.sites)
        if (!std::isfinite(site.signal) || !std::isfinite(site.variance))
            throw std::runtime_error("merged CFA signal or variance overflows float32");
    time.mark("flat_field");

    std::unique_ptr<StagedFile> cfaFile;
    if (!o.writeCfa.empty()) {
        cfaFile = std::make_unique<StagedFile>(o.writeCfa);
        writeCfa(cfaFile->path(), ref.layout, merge.sites, o, sources);
        time.mark("write_cfa");
    }
    repairSites(merge.sites, ref.layout, colour.wb, o, sources);
    double sum = 0;
    for (const auto& frame : plain) sum += std::log(frame.effective);
    const double scale = o.normalize == "none" ? 1 : std::exp(sum / plain.size());
    if (o.verbose) std::cerr << "Demosaicing and converting colour\n";
    const auto output = render(merge.sites, ref.layout, colour, scale, o, sources);
    time.mark("render");

    const auto stats = mergeReport(merge, output);
    std::vector<Site>().swap(merge.sites);
    StagedFile imageFile(o.output), reportFile(o.report);
    writeImage(imageFile.path(), ref.layout, output, ref, colour, nd, scale, referenceAperture, o);
    time.mark("write_exr");

    Json inputs = Json::array();
    Json groupsReport = Json::array();
    auto append = [&](const std::vector<RawFrame>& group, const char* role, const InputGroup* spec) {
        std::vector<Json> byInput(group.size());
        for (const auto& f : group) {
            auto input = inputReport(f, role, o);
            if (spec) input["group"] = spec->name;
            byInput[f.inputIndex] = std::move(input);
        }
        for (auto& f : byInput) inputs.push_back(std::move(f));
    };
    if (grouped) {
        for (size_t g = 0; g < o.groups.size(); ++g) {
            const auto& spec = o.groups[g];
            const auto& result = groupResults[g];
            append(storage[g], spec.ev || spec.autoEV ? "nd" : "plain", &spec);
            groupsReport.push_back({{"name", spec.name}, {"source_id", g + 1},
                {"valid_mask", spec.validMask.string()}, {"contribution_mask", spec.contributionMask.string()},
                {"ev", rgbReport(result.ev)}, {"calibration_reference", spec.calibrationReference},
                {"calibration", result.calibration ? ndReport(*result.calibration, storage[g].size(),
                    bool(spec.ev), spec.calibrationMask, "calibration") : Json(nullptr)},
                {"measured_sites", result.measuredSites}, {"lower_bound_sites", result.boundSites},
                {"dark_model", o.darkModel.empty() ? Json(nullptr) : darkReport(darkFits[g])}});
        }
    } else {
        append(plain, "plain", nullptr);
        append(storage[1], "nd", nullptr);
    }
    const Json report{
        {"rawfuse_version", RAWFUSE_VERSION}, {"inputs", std::move(inputs)},
        {"reference_frame", ref.inputIndex + 1}, {"width", ref.layout.width}, {"height", ref.layout.height},
        {"cfa", ref.layout.pattern()}, {"make", ref.make}, {"model", ref.model},
        {"space", o.space == "aces" ? "aces-ap0" : "camera"},
        {"white_balance", {{"method", o.wb}, {"multipliers", colour.wb}, {"spread_across_frames", spread}}},
        {"matrix", {{"method", o.matrix}, {"values", colour.matrix ? Json(*colour.matrix) : Json(nullptr)}}},
        {"exposure_reference", scale}, {"reference_f_number", referenceAperture},
        {"noise_profile", noiseReport(ref, o)},
        {"output_floor", o.outputFloor > 0 ? Json(o.outputFloor) : Json(nullptr)},
        {"dark_model", o.darkModel.empty() ? Json(nullptr) : Json{
            {"path", o.darkModel.string()}, {"plain", darkReport(darkFits[0])},
            {"nd", nd ? darkReport(darkFits[1]) : Json(nullptr)}}},
        {"vignette", {{"source", field.source}, {"gain_maps", field.maps}, {"max_gain", rgbReport(field.maxGain)}}},
        {"nd", nd ? ndReport(*nd, storage[1].size(), bool(o.ndEv), o.ndMask, o.ndUse) : Json(nullptr)},
        {"groups", grouped ? std::move(groupsReport) : Json(nullptr)}, {"merge", stats},
        {"demosaic", o.demosaic}, {"warnings", warnings}, {"timing_seconds", time.values}};
    {
        std::ofstream stream(reportFile.path(), std::ios::binary);
        stream.exceptions(std::ios::failbit | std::ios::badbit);
        stream << report.dump(2, ' ', false, Json::error_handler_t::replace) << '\n';
        stream.close();
    }

    try {
        if (cfaFile) cfaFile->publish(o.force);
        reportFile.publish(o.force);
        imageFile.publish(o.force);
    } catch (...) {
        reportFile.rollback();
        if (cfaFile) cfaFile->rollback();
        throw;
    }
}

void printWarnings(const Warnings& warnings) {
    for (const auto& warning : warnings) std::cerr << "rawfuse: warning: " << warning << '\n';
}

} // namespace

void run(const Options& o) {
    Warnings warnings;
    try {
        fuse(o, warnings);
    } catch (...) {
        // Context gathered before a failure often explains it.
        printWarnings(warnings);
        throw;
    }
    printWarnings(warnings);
}
}

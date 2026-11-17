#include "pipeline.h"

#include <OpenEXR/ImfThreading.h>
#include <OpenImageIO/imageio.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <random>
#include <sstream>
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
    for (const auto& group : o.groups) addFrames(group.files);
    for (const auto& group : o.groups) {
        for (const auto& mask : {group.validMask, group.contributionMask, group.calibrationMask, group.driftMask})
            if (!mask.empty()) inputs.push_back(mask);
    }
    for (const auto& sidecar : {o.darkModel, o.manifest})
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

Json noiseReport(const RawFrame& f) {
    if (f.noise[0].shot == f.noise[1].shot && f.noise[0].shot == f.noise[2].shot &&
        f.noise[0].read == f.noise[1].read && f.noise[0].read == f.noise[2].read)
        return {{"source", f.noiseSource}, {"S", f.noise[0].shot}, {"O", f.noise[0].read}};
    return {{"source", f.noiseSource}, {"planes", {
        {"R", {{"S", f.noise[0].shot}, {"O", f.noise[0].read}}},
        {"G", {{"S", f.noise[1].shot}, {"O", f.noise[1].read}}},
        {"B", {{"S", f.noise[2].shot}, {"O", f.noise[2].read}}}}}};
}

Json optionalNumber(double value) { return std::isfinite(value) && value > 0 ? Json(value) : Json(nullptr); }

Json inputReport(const RawFrame& f, const char* role) {
    const auto& m = f.metadata;
    return {{"path", f.path.string()}, {"role", role}, {"index", f.inputIndex + 1},
            {"exposure_time", m.exposure}, {"seconds", m.seconds}, {"iso", m.iso},
            {"f_number", m.aperture}, {"exif_effective_exposure", f.exifEffective},
            {"effective_exposure", f.effective},
            {"exposure_correction_ev", std::log2(f.effective / f.exifEffective)},
            {"as_shot_neutral", m.neutral ? Json(*m.neutral) : Json(nullptr)},
            {"image_description", m.description}, {"xmp", m.xmp}, {"capture_date", m.captureDate},
            {"baseline_exposure", m.baseline}, {"noise_profile", noiseReport(f)},
            {"optical_black_sigma_dn", optionalNumber(f.opticalBlackSigma)},
            {"black_level", f.blackPattern}, {"white_level", f.white},
            {"hot_sites", {{"kept", f.hotKept}, {"missing", f.hotMissing}, {"model_rejected", f.hotRejected}}}};
}

Json stepsReport(const std::vector<ExposureStep>& steps) {
    Json result = Json::array();
    for (const auto& step : steps)
        result.push_back({{"from", step.from}, {"to", step.to}, {"exif_ratio", step.exifRatio},
                          {"measured_ratio", optionalNumber(step.measuredRatio)}, {"samples", step.samples},
                          {"ratio", step.ratio}, {"applied", step.applied},
                          {"borrowed_from", step.borrowedFrom.empty() ? Json(nullptr) : Json(step.borrowedFrom)}});
    return result;
}

Json noiseFitsReport(const std::vector<NoiseFit>& fits) {
    Json result = Json::array();
    for (const auto& fit : fits)
        result.push_back({{"iso", fit.iso}, {"S", fit.shot}, {"O", fit.read}, {"fit_points", fit.bins},
                          {"read_noise", fit.opticalBlack ? "optical_black" : "bracket"}});
    return result;
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
        warnings.push_back("AsShotNeutral varies across the first group's frames; applying one white balance after fusion");
    return {{"min", rgbReport(low)}, {"max", rgbReport(high)}, {"ev", rgbReport(spread)}};
}

Json darkReport(const DarkFit& fit) {
    auto optional = [](const std::optional<double>& value) { return value ? Json(*value) : Json(nullptr); };
    return {{"method", fit.method}, {"offset_scale", optional(fit.offsetScale)},
            {"rate_scale", optional(fit.rateScale)}, {"relative_scatter", optional(fit.relativeScatter)},
            {"relative_uncertainty", fit.uncertainty}, {"sites", fit.sites}, {"observations", fit.observations}};
}

Json calibrationReport(const NDResult& nd, bool manual, const InputGroup& group) {
    Json blocks = Json::array(), flagged = Json::array();
    for (const auto& block : nd.blocks) {
        blocks.push_back({{"row", block.row}, {"col", block.col}, {"ev", rgbReport(block.ev)},
                          {"samples", block.samples}});
        if (block.flagged) flagged.push_back({{"row", block.row}, {"col", block.col}});
    }
    return {{"method", manual ? "manual" : "auto"}, {"ev", rgbReport(nd.ev)},
            {"measured_ev", rgbReport(nd.measured)}, {"mad_ev", rgbReport(nd.mad)},
            {"samples", {{"R", nd.samples[0]}, {"G", nd.samples[1]}, {"B", nd.samples[2]}}},
            {"blocks_flagged", std::move(flagged)}, {"blocks", std::move(blocks)},
            {"scene_drift_ev", nd.sceneDrift ? Json(*nd.sceneDrift) : Json(nullptr)},
            {"scene_drift_samples", nd.sceneDriftSamples},
            {"calibration_mask", group.calibrationMask.string()}, {"drift_mask", group.driftMask.string()}};
}

// Offsets are reported in camera RGB at the output exposure scale.
Json matchReport(const std::optional<OffsetMatch>& match, double scale) {
    if (!match) return nullptr;
    const RGB median{match->median[0] * scale, match->median[1] * scale, match->median[2] * scale};
    return {{"overlap_samples", match->samples}, {"known_cells", match->knownCells},
            {"corrected_sites", match->corrected}, {"unresolved_sites", match->unresolved},
            {"capped_sites", match->capped},
            {"median_offset", rgbReport(median)}, {"median_share", rgbReport(match->share)}};
}

Json mergeReport(const MergeResult& merge, const Output& output) {
    std::array<size_t, 7> counts{};
    size_t unsupported = 0, adjacent = 0, outside = 0, rebuilt = 0;
    for (const auto& site : merge.sites) {
        for (int bit = 0; bit < 6; ++bit) counts[bit] += (site.flags >> bit) & 1;
        unsupported += !site.frames;
        adjacent += (site.flags & adjacentClipFallback) != 0;
        outside += (site.flags & uncovered) != 0;
        rebuilt += (site.flags & reconstructed) != 0;
    }
    const double size = double(merge.sites.size());
    return {{"lower_bound_fraction", counts[0] / size}, {"from_nd_fraction", counts[1] / size},
            {"moving_fraction", counts[2] / size}, {"filled_fraction", counts[3] / size},
            {"missing_fraction", counts[4] / size}, {"inconsistent_fraction", counts[5] / size},
            {"shadow_fraction", merge.shadow / size}, {"unsupported_fraction", unsupported / size},
            {"adjacent_clip_fallback_fraction", adjacent / size},
            {"uncovered_fraction", outside / size}, {"reconstructed_fraction", rebuilt / size},
            {"floored_fraction", output.flooredPixels / size}};
}

std::string grouped(size_t n) {
    auto text = std::to_string(n);
    for (int i = int(text.size()) - 3; i > 0; i -= 3) text.insert(size_t(i), ",");
    return text;
}

std::string fixed(double value, int digits) {
    std::ostringstream text;
    text << std::fixed << std::setprecision(digits) << value;
    return text.str();
}

std::string seconds(double value) {
    if (value >= 0.3) return fixed(value, value >= 10 ? 0 : 1) + " s";
    std::ostringstream text;
    text << "1/" << std::setprecision(4) << 1 / value << " s";
    return text.str();
}

std::string evText(const RGB& ev) {
    if (ev[0] == ev[1] && ev[1] == ev[2]) return fixed(ev[1], 2) + " EV";
    return "R " + fixed(ev[0], 2) + " G " + fixed(ev[1], 2) + " B " + fixed(ev[2], 2) + " EV";
}

std::string percent(double fraction) {
    if (fraction == 0) return "0%";
    return fixed(100 * fraction, fraction < 0.001 ? 4 : fraction < 0.01 ? 3 : 1) + "%";
}

// Progress and outcome lines on stderr; --quiet keeps only warnings and errors.
struct Status {
    const Options& o;
    void operator()(const std::string& line) const {
        if (!o.quiet) std::cerr << "rawfuse: " << line << '\n';
    }
};

void reportBrackets(const Status& status, const FrameGroups& frames, const Preparation& preparation,
                    const Options& o) {
    for (size_t g = 0; g < frames.size(); ++g) {
        const auto& group = *frames[g];
        std::ostringstream line;
        line << o.groups[g].name << ": " << group.size() << (group.size() == 1 ? " frame" : " frames");
        if (group.size() > 1)
            line << ", " << seconds(group.front().metadata.seconds) << " to " << seconds(group.back().metadata.seconds);
        else line << ", " << seconds(group.front().metadata.seconds);
        const auto& steps = preparation.exposureSteps[g];
        if (o.exposures == "exif") line << "; EXIF exposures";
        else if (!steps.empty()) {
            size_t measured = 0, borrowed = 0;
            double largest = 0;
            for (const auto& step : steps) {
                measured += step.applied;
                borrowed += !step.applied && !step.borrowedFrom.empty();
                largest = std::max(largest, std::abs(std::log2(step.ratio / step.exifRatio)));
            }
            line << "; exposure steps measured " << measured << " of " << steps.size();
            if (borrowed) line << ", " << borrowed << " borrowed";
            if (steps.size() > measured + borrowed) line << ", " << steps.size() - measured - borrowed << " EXIF";
            line << " (largest correction " << fixed(largest, 2) << " EV)";
        }
        status(line.str());
    }
    const auto& noise = frames.front()->front();
    std::ostringstream line;
    line << "noise: ";
    if (noise.noiseSource == "dng") line << "DNG NoiseProfile";
    else if (noise.noiseSource == "custom") line << "--noise";
    else if (noise.noiseSource == "generic") line << "generic estimate (nothing to measure)";
    else {
        line << "measured from frame pairs";
        for (const auto& fit : preparation.noiseFits)
            line << "; ISO " << fixed(fit.iso, 0) << " shot " << std::setprecision(2) << fit.shot
                 << ", read " << fit.read << (fit.opticalBlack ? " (optical black)" : "");
    }
    status(line.str());
}

void reportComposition(const Status& status, const std::vector<GroupResult>& results, double scale,
                       const Options& o) {
    for (size_t g = 0; g < results.size(); ++g) {
        const auto& spec = o.groups[g];
        const auto& result = results[g];
        if (result.calibration) {
            const auto& nd = *result.calibration;
            std::ostringstream line;
            line << spec.name << ": " << (spec.ev ? "ev " + evText(nd.ev) + " (manual), checked against "
                                                 : "ev " + evText(nd.ev) + " measured against ")
                 << spec.calibrationReference << " on " << grouped(nd.samples[0]) << "/" << grouped(nd.samples[1])
                 << "/" << grouped(nd.samples[2]) << " R/G/B sites";
            if (spec.ev) line << " (measured " << evText(nd.measured) << ")";
            if (nd.sceneDrift) line << "; scene drift " << fixed(*nd.sceneDrift, 2) << " EV";
            status(line.str());
        } else if (spec.ev) status(spec.name + ": ev " + evText(*spec.ev) + " (manual, unchecked)");
        std::ostringstream line;
        line << spec.name << ": supplies " << grouped(result.measuredSites) << " sites";
        if (result.boundSites) line << " and " << grouped(result.boundSites) << " clipped lower bounds";
        if (!g) line << " (reference)";
        else line << (spec.contributionMask.empty() ? " (fills gaps)" : " (replaces inside its contribution mask, fills elsewhere)");
        status(line.str());
        if (result.match) {
            const auto& m = *result.match;
            std::ostringstream match;
            match << spec.name << ": ";
            if (!m.corrected && !m.unresolved) match << "offset match: nothing to correct";
            else if (!m.knownCells) match << "offset match: no precise overlap with earlier groups; not corrected";
            else {
                match << "offset matched to earlier groups from " << grouped(m.samples)
                      << " overlap sites; subtracted a median " << std::setprecision(3) << m.median[1] * scale
                      << " (G, camera RGB at output scale), " << fixed(100 * m.share[1], 0) << "% of its signal, at "
                      << grouped(m.corrected) << " sites";
                if (m.capped) match << " (" << grouped(m.capped) << " capped at their full value)";
                if (m.unresolved) match << ", " << grouped(m.unresolved) << " beyond reach";
            }
            status(match.str());
        }
    }
}

void fuse(const Options& o, Warnings& warnings) {
    Timings time;
    const Status status{o};
    validatePaths(o);
    ThreadSettings threads(o.threads);
    std::vector<std::vector<RawFrame>> storage(o.groups.size());
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
    for (size_t g = 0; g < o.groups.size(); ++g) read(o.groups[g].files, storage[g]);
    const double referenceAperture = validateFrames(frames, warnings);
    auto byExposure = [](const auto& a, const auto& b) { return a.effective < b.effective; };
    for (auto* group : frames) std::stable_sort(group->begin(), group->end(), byExposure);
    time.mark("decode");
    {
        size_t count = 0;
        for (const auto* group : frames) count += group->size();
        const auto& first = plain.front();
        status("read " + std::to_string(count) + " frames in " + std::to_string(frames.size()) +
               (frames.size() == 1 ? " group" : " groups") + " from " + first.make + " " + first.model + " (" +
               std::to_string(first.layout.width) + "x" + std::to_string(first.layout.height) + " " +
               first.layout.pattern() + ")");
    }

    const auto preparation = prepareFrames(frames, *sidecars, o, warnings);
    time.mark("prepare");
    reportBrackets(status, frames, preparation, o);
    if (o.verbose) std::cerr << "Merging Bayer samples\n";
    auto composition = mergeGroups(frames, *sidecars, o, warnings);
    auto merge = std::move(composition.mosaic);
    auto sources = std::move(composition.sources);
    const auto groupResults = std::move(composition.groups);
    time.mark("merge");
    double sum = 0;
    for (const auto& frame : plain) sum += std::log(frame.effective);
    const double scale = o.normalize == "none" ? 1 : std::exp(sum / plain.size());
    reportComposition(status, groupResults, scale, o);
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
    if (o.verbose) std::cerr << "Demosaicing and converting colour\n";
    const auto output = render(merge.sites, ref.layout, colour, scale, o, sources);
    time.mark("render");

    const auto stats = mergeReport(merge, output);
    std::vector<Site>().swap(merge.sites);
    StagedFile imageFile(o.output), reportFile(o.report);
    writeImage(imageFile.path(), ref.layout, output, ref, colour, scale, referenceAperture, o);
    time.mark("write_exr");

    Json inputs = Json::array();
    Json groupsReport = Json::array();
    for (size_t g = 0; g < o.groups.size(); ++g) {
        const auto& spec = o.groups[g];
        const auto& result = groupResults[g];
        std::vector<Json> byInput(storage[g].size());
        for (const auto& f : storage[g]) {
            auto input = inputReport(f, spec.ev || spec.autoEV ? "attenuated" : "plain");
            input["group"] = spec.name;
            byInput[f.inputIndex] = std::move(input);
        }
        for (auto& f : byInput) inputs.push_back(std::move(f));
        groupsReport.push_back({{"name", spec.name}, {"source_id", g + 1}, {"frames", storage[g].size()},
            {"valid_mask", spec.validMask.string()}, {"contribution_mask", spec.contributionMask.string()},
            {"ev", rgbReport(result.ev)}, {"calibration_reference", spec.calibrationReference},
            {"calibration", result.calibration ?
                calibrationReport(*result.calibration, bool(spec.ev), spec) : Json(nullptr)},
            {"exposure_steps", stepsReport(preparation.exposureSteps[g])},
            {"measured_sites", result.measuredSites}, {"lower_bound_sites", result.boundSites},
            {"offset_match", matchReport(result.match, scale)},
            {"dark_model", o.darkModel.empty() ? Json(nullptr) : darkReport(preparation.darkFits[g])}});
    }
    const Json report{
        {"rawfuse_version", RAWFUSE_VERSION}, {"inputs", std::move(inputs)},
        {"reference_frame", ref.inputIndex + 1}, {"width", ref.layout.width}, {"height", ref.layout.height},
        {"cfa", ref.layout.pattern()}, {"make", ref.make}, {"model", ref.model},
        {"space", o.space == "aces" ? "aces-ap0" : "camera"},
        {"white_balance", {{"method", o.wb}, {"multipliers", colour.wb}, {"spread_across_frames", spread}}},
        {"matrix", {{"method", o.matrix}, {"values", colour.matrix ? Json(*colour.matrix) : Json(nullptr)}}},
        {"exposure_reference", scale}, {"reference_f_number", referenceAperture},
        {"exposures", o.exposures},
        {"noise_profile", noiseReport(ref)}, {"noise_fits", noiseFitsReport(preparation.noiseFits)},
        {"output_floor", o.outputFloor > 0 ? Json(o.outputFloor) : Json(nullptr)},
        {"dark_model", o.darkModel.empty() ? Json(nullptr) : Json(o.darkModel.string())},
        {"vignette", {{"source", field.source}, {"gain_maps", field.maps}, {"max_gain", rgbReport(field.maxGain)}}},
        {"groups", std::move(groupsReport)}, {"merge", stats},
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
    std::ostringstream outcome;
    outcome << "result: " << percent(stats["missing_fraction"].get<double>()) << " of sites missing, "
            << percent(stats["lower_bound_fraction"].get<double>()) << " clipped lower bounds, "
            << percent(stats["moving_fraction"].get<double>()) << " motion-rejected samples, "
            << percent(stats["inconsistent_fraction"].get<double>()) << " inconsistent";
    status(outcome.str());
    status("wrote " + o.output.string() + " (" + (o.space == "aces" ? "ACES AP0" : "camera RGB") + ") and " +
           o.report.string() + (o.writeCfa.empty() ? "" : " and " + o.writeCfa.string()) + " in " +
           fixed(time.values["total"].get<double>(), 1) + " s" +
           (warnings.empty() ? "" : "; " + std::to_string(warnings.size()) +
                                (warnings.size() == 1 ? " warning follows" : " warnings follow")));
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

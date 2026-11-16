#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <memory>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace rawfuse {

namespace fs = std::filesystem;
using RGB = std::array<double, 3>;
using Matrix = std::array<RGB, 3>;
using Warnings = std::vector<std::string>;

constexpr RGB unitWB{1, 1, 1};
constexpr Matrix identity{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

struct InputGroup {
    std::string name, calibrationReference;
    std::vector<fs::path> files;
    fs::path validMask, contributionMask, calibrationMask, driftMask;
    std::optional<RGB> ev;
    bool autoEV = false;
    bool matchOffset = true;  // "match": "offset" (default) or "none"
};

struct Options {
    std::vector<fs::path> plain;  // Positional frames; parsed into a single group
    fs::path manifest;
    std::vector<InputGroup> groups;
    fs::path output, report, writeCfa, darkModel;
    std::string space = "aces", wb = "median", matrix = "metadata";
    std::string normalize = "geomean", badPixels = "dng", noise = "auto";
    std::string exposures = "measured";
    std::optional<Matrix> customMatrix;  // Parsed --matrix custom:FILE
    std::string vignette = "dng", highlights = "neutral", demosaic = "bilinear";
    double outputFloor = 0;
    int reference = 0, clipDilate = 1;
    unsigned threads = 1;
    double deghostEv = 1, darkModelError = 0.10, low = 0.0083, high = 0.88;
    double clipMargin = 0.001;
    std::vector<double> noiseOverride;  // Parsed --noise S,O
    bool deghost = true, aux = false, force = false, verbose = false, quiet = false, allowSceneDrift = false;
};

struct Point { int row, col; };
struct Rect { int top, left, bottom, right; };
struct GainMap {
    Rect area;
    int rowPitch, colPitch, rows, cols;
    double spacingV, spacingH, originV, originH;
    std::vector<float> gains;
    bool operator==(const GainMap& other) const;
};

struct Metadata {
    std::string exposure, description, xmp, captureDate;
    double seconds = 0, iso = 0, aperture = 0, baseline = 0;
    std::optional<RGB> neutral;
    std::optional<Rect> activeArea;
    std::optional<std::array<int, 4>> cfa;
    std::array<int, 3> planeColours{0, 1, 2};
    std::vector<double> noise, black, white, blackRows, blackCols;
    int blackRepeatRows = 1, blackRepeatCols = 1;
    std::vector<Point> badPoints;
    std::vector<Rect> badRects;
    std::vector<GainMap> maps;
};

Metadata readDngMetadata(const fs::path&, Warnings&);

struct Layout {
    int width = 0, height = 0;
    std::array<int, 4> cfa{};
    size_t size() const { return size_t(width) * height; }
    size_t index(int row, int col) const { return size_t(row) * width + col; }
    int colour(int row, int col) const { return cfa[(row & 1) * 2 + (col & 1)]; }
    std::string pattern() const;
};

struct Noise { double shot, read; };
enum FrameMask : uint8_t { saturated = 1, bloom = 2, excluded = 4, modeled = 8, badRect = 16, outsideValidity = 32 };
struct RawFrame {
    fs::path path;
    Layout layout;
    Metadata metadata;
    std::string make, model;
    int top = 0, left = 0, rawWidth = 0, rawHeight = 0, stride = 0;
    int blackRows = 1, blackCols = 1;
    std::vector<float> blackPattern, blackDeltaRows, blackDeltaCols;
    float white = 0;
    double effective = 0, exifEffective = 0;  // Measured (or EXIF) and EXIF-derived exposure
    double opticalBlackSigma = 0;  // Robust DN spread of masked pixels; 0 when unavailable
    std::array<Noise, 3> noise{};
    std::string noiseSource;
    RGB wb = unitWB;
    bool dng = false, hasWB = false;
    Matrix cameraToSRGB{};
    size_t inputIndex = 0;
    long long captureTimestamp = 0;
    std::shared_ptr<const void> mosaicOwner;
    const uint16_t* mosaic = nullptr;
    std::vector<uint8_t> rejected;
    std::vector<std::pair<size_t, float>> dark;
    size_t hotKept = 0, hotMissing = 0, hotRejected = 0, undatedDefects = 0;
    double darkUncertainty = 0;

    uint16_t raw(int row, int col) const { return mosaic[size_t(row) * stride + col]; }
    float black(int row, int col) const {
        if (blackPattern.size() == 1 && blackDeltaRows.empty() && blackDeltaCols.empty())
            return blackPattern.front();
        float value = blackPattern[size_t(row % blackRows) * blackCols + col % blackCols];
        if (!blackDeltaRows.empty()) value += blackDeltaRows[row];
        if (!blackDeltaCols.empty()) value += blackDeltaCols[col];
        return value;
    }
};

RawFrame readRaw(const fs::path&, Warnings&);
using FrameGroups = std::vector<std::vector<RawFrame>*>;
using Mask = std::vector<uint8_t>;
double validateFrames(const FrameGroups&, Warnings&);
using DarkModel = std::map<size_t, std::array<float, 2>>;
struct DarkFit {
    std::string method = "unavailable";
    std::optional<double> offsetScale, rateScale, relativeScatter;
    double uncertainty = 0;
    size_t sites = 0, observations = 0;
};
DarkFit fitDarkModel(std::vector<RawFrame>&, const DarkModel&, const Options&);
struct BadPixel { Point site; long long since; };

// Caller-supplied files that depend on the sensor grid. They are loaded and
// validated after the first frame is decoded, before the remaining frames.
struct Sidecars {
    DarkModel darkModel;
    std::vector<BadPixel> badPixels;
    std::vector<Mask> valid, contribution, calibration, drift;  // Per input group
};
Sidecars loadSidecars(const Options&, const Layout&);

// Relative exposures measured between neighbouring frames of one group.
struct ExposureStep {
    size_t from = 0, to = 0;  // One-based input indices within the group
    double exifRatio = 0, measuredRatio = 0, ratio = 0;  // ratio: the one used
    double fromSeconds = 0, toSeconds = 0;
    size_t samples = 0;
    bool applied = false;      // measuredRatio is used
    std::string borrowedFrom;  // Group whose measurement for the same settings is used
};
// Noise fitted from temporal differences between neighbouring frames, per ISO.
struct NoiseFit {
    double iso = 0, shot = 0, read = 0;
    size_t bins = 0;
    bool opticalBlack = false;
};
struct Preparation {
    std::vector<DarkFit> darkFits;
    std::vector<std::vector<ExposureStep>> exposureSteps;  // Per group
    std::vector<NoiseFit> noiseFits;
};
Preparation prepareFrames(const FrameGroups&, const Sidecars&, const Options&, Warnings&);
// bracket.cpp. Both need clipping and defect masks. measureExposures rescales
// `effective` of exposure-sorted groups; assignNoise uses those exposures for
// its temporal differences.
std::vector<std::vector<ExposureStep>> measureExposures(const FrameGroups&, const Options&, Warnings&);
std::vector<NoiseFit> assignNoise(const FrameGroups&, const Options&, Warnings&);

constexpr uint32_t lowerBound = 1u << 0, fromND = 1u << 1, moving = 1u << 2,
                   filled = 1u << 3, missing = 1u << 4, inconsistent = 1u << 5,
                   floored = 1u << 6, adjacentClipFallback = 1u << 7, uncovered = 1u << 8;
constexpr uint32_t mixedSource = 65535;
struct Site {
    float signal = 0, variance = 0;
    uint16_t flags = 0, frames = 0;
};
enum class SampleState { usable, aboveHigh, clipped, adjacentClip, missing };
struct Sample {
    SampleState state;
    float fraction = 0, extra = 0;
};

// Each frame uses its own black and white levels; they may differ within a bracket.
Sample sample(const RawFrame&, int row, int col, size_t index, const Options&);
struct MergeResult {
    std::vector<Site> sites;
    size_t reference = 0;
    size_t shadow = 0;
    Mask shadowSites;
};
MergeResult mergePlain(const std::vector<RawFrame>&, const Options&, const Mask& selection = {},
                       bool trackShadow = false);
struct SiteMerge {
    Site site;
    bool shadow = false;
};
SiteMerge mergeSite(const std::vector<RawFrame>&, int row, int col, const Options&);

struct BlockCheck { int row, col; RGB ev; std::array<size_t, 3> samples; bool flagged; };
struct NDResult {
    RGB ev{}, measured{}, mad{};
    std::array<size_t, 3> samples{};
    std::optional<double> sceneDrift;
    size_t sceneDriftSamples = 0;
    std::vector<BlockCheck> blocks;
};
// Measures the attenuation of `filtered` against the merged `reference` group on
// `calibration` sites (empty: all), and scene drift on `drift` sites, which the
// filtered frames record without attenuation.
NDResult calibrateND(const MergeResult& merged, const std::vector<RawFrame>& reference,
                     const std::vector<RawFrame>& filtered, const Mask& calibration, const Mask& drift,
                     const std::optional<RGB>& manualEV, const Options&, Warnings&);

// Additive offset between a group and the composite of earlier groups (glare
// veil, flare or a small level mismatch), measured where both are precise and
// interpolated into the sites the group supplies.
struct OffsetMatch {
    size_t samples = 0, knownCells = 0, corrected = 0, unresolved = 0, capped = 0;
    RGB median{}, share{};  // Applied offset and its share of the uncorrected signal (medians)
};

struct GroupResult {
    RGB ev{};
    std::optional<NDResult> calibration;
    std::optional<OffsetMatch> match;
    size_t measuredSites = 0, boundSites = 0;
};

struct GroupMerge {
    MergeResult mosaic;
    std::vector<uint32_t> sources;
    std::vector<GroupResult> groups;
};
std::vector<InputGroup> readGroups(const fs::path&);
GroupMerge mergeGroups(const FrameGroups&, const Sidecars&, const Options&, Warnings&);
Mask readMask(const fs::path&, const Layout&);
void releaseMosaics(std::vector<RawFrame>&);

struct FlatField {
    std::string source;
    size_t maps = 0;
    RGB maxGain = unitWB;
};
FlatField flatField(std::vector<Site>&, const RawFrame&, const Options&);
void repairSites(std::vector<Site>&, const Layout&, const RGB&, const Options&,
                 std::vector<uint32_t>& sources);

struct Colour {
    RGB wb = unitWB;
    std::optional<Matrix> matrix;
};
Colour configureColour(const std::vector<RawFrame>&, size_t, const Options&, Warnings&);
struct Output {
    std::vector<float> rgb, variance;
    std::vector<uint16_t> frames;
    std::vector<uint32_t> flags, sources;
    size_t flooredPixels = 0;
};
Output render(const std::vector<Site>&, const Layout&, const Colour&, double scale,
              const Options&, const std::vector<uint32_t>& sources);
std::vector<float> readImage(const fs::path&, const Layout&, int& channels);
int imageChannels(const fs::path&, const Layout&);
void writeImage(const fs::path&, const Layout&, const Output&, const RawFrame&, const Colour&,
                double scale, double referenceAperture, const Options&);
void writeCfa(const fs::path&, const Layout&, const std::vector<Site>&, const Options&,
              const std::vector<uint32_t>& sources);

void printHelp();
Options parseOptions(int, char**);
void run(const Options&);

double number(const std::string&);
std::vector<double> numberList(const std::string&, size_t expected);
double median(std::vector<double>&);  // Requires a nonempty input
// Variance of the central distribution: a MAD start, then a 4σ-clipped second
// moment, which stays accurate for integer-quantized data where the MAD does not.
double robustVariance(std::vector<double>);  // Requires a nonempty input

template<class Function>
void parallelRows(int rows, unsigned threads, Function function) {
    const unsigned workers = std::min(threads, unsigned((rows + 7) / 8));
    if (workers < 2) {
        for (int row = 0; row < rows; ++row) function(row);
        return;
    }
    std::atomic<int> next{0};
    std::mutex errorMutex;
    std::exception_ptr failure;
    auto work = [&] {
        try {
            for (int start = next.fetch_add(8); start < rows; start = next.fetch_add(8)) {
                for (int row = start; row < std::min(start + 8, rows); ++row) function(row);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(errorMutex);
            if (!failure) failure = std::current_exception();
            next.store(rows);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    try {
        for (unsigned i = 1; i < workers; ++i) pool.emplace_back(work);
    } catch (...) {
        next.store(rows);
        for (auto& thread : pool) thread.join();
        throw;
    }
    work();
    for (auto& thread : pool) thread.join();
    if (failure) std::rethrow_exception(failure);
}
}

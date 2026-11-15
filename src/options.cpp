#include "pipeline.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace rawfuse {

double number(const std::string& value) {
    size_t used = 0;
    double result;
    try { result = std::stod(value, &used); }
    catch (const std::exception&) { throw std::runtime_error("expected a number, got " + value); }
    if (used != value.size() || !std::isfinite(result))
        throw std::runtime_error("expected a finite number, got " + value);
    return result;
}

std::vector<double> numberList(const std::string& text, size_t expected) {
    std::vector<double> result;
    size_t start = 0;
    for (;;) {
        const size_t end = text.find(',', start);
        result.push_back(number(text.substr(start, end - start)));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    if (result.size() != expected)
        throw std::runtime_error("expected " + std::to_string(expected) + " comma-separated numbers");
    return result;
}

double median(std::vector<double>& values) {
    if (values.empty()) throw std::logic_error("median of an empty set");
    auto middle = values.begin() + values.size() / 2;
    std::nth_element(values.begin(), middle, values.end());
    if (values.size() % 2) return *middle;
    return (*std::max_element(values.begin(), middle) + *middle) / 2;
}

void printHelp() {
    std::cout << R"(Usage: rawfuse [options] FRAME... [--nd FRAME...] -o OUT.exr

       rawfuse [options] --groups FILE -o OUT.exr

Merge fixed-camera Bayer exposures into a linear float OpenEXR.
  --groups FILE                       JSON manifest of brackets and masks
  -o, --output PATH                   Output EXR (required)
  --report PATH                       JSON report (default: OUT.json)
  --force                             Replace existing output files
  --space aces|camera                 Output colour space (aces)
  --wb reference|median|frame:N|custom:R,G,B  (median)
  --matrix metadata|auto|spectral|adobe|custom:FILE
  --normalize geomean|none
  --floor F                           Positive minimum in output units (off by default)
  --reference auto|N                  Frame number in the plain bracket or first group
  --no-deghost                        Disable motion rejection
  --deghost-ev EV                     Motion threshold (1)
  --bad-pixels dng|none|FILE           Defect list or DNG opcodes
  --dark-model CSV                    Per-site dark correction model
  --dark-model-error F                Minimum relative model uncertainty (0.10)
  --noise dng|S,O                     Normalized shot/read noise coefficients
  --low F --high F                    Signal range (0.0083, 0.88)
  --clip-margin F --clip-dilate N      Saturation margin and radius (0.001, 1)
  --vignette dng|none|FILE             DNG gain maps or float flat-field EXR
  --highlights neutral|bound
  --demosaic bilinear|mhc
  --nd FRAME...                       ND exposures up to the next option
  --nd-ev auto|EV                     Estimate RGB transmission or set EV
  --nd-mask FILE                      Full-frame coverage image
  --nd-use fill|all                   ND contribution (fill)
  --nd-min-signal F                   Minimum ND signal for estimation (0.05)
  --allow-nd-scene-drift              Permit measured drift with manual ND EV
  --write-cfa PATH                    Diagnostic mosaic EXR
  --aux                               Include variance, frame count and flags
  --threads N                         Worker limit (all logical cores)
  -v                                  Show progress
)";
}

static int integer(const std::string& value, int min, int max) {
    const double n = number(value);
    if (n < min || n > max || n != std::floor(n))
        throw std::runtime_error("invalid integer: " + value);
    return int(n);
}

static void oneOf(const std::string& option, const std::string& value,
                  std::initializer_list<const char*> choices) {
    if (std::none_of(choices.begin(), choices.end(), [&](const char* x) { return value == x; }))
        throw std::runtime_error("invalid " + option + ": " + value);
}

Options parseOptions(int argc, char** argv) {
    Options o;
    o.threads = std::max(1u, std::thread::hardware_concurrency());
    const std::unordered_map<std::string, std::string Options::*> strings{
        {"--space", &Options::space}, {"--wb", &Options::wb}, {"--matrix", &Options::matrix},
        {"--normalize", &Options::normalize}, {"--bad-pixels", &Options::badPixels},
        {"--noise", &Options::noise}, {"--vignette", &Options::vignette},
        {"--highlights", &Options::highlights}, {"--demosaic", &Options::demosaic},
        {"--nd-use", &Options::ndUse}};
    const std::unordered_map<std::string, fs::path Options::*> paths{
        {"-o", &Options::output}, {"--output", &Options::output}, {"--report", &Options::report},
        {"--write-cfa", &Options::writeCfa}, {"--dark-model", &Options::darkModel},
        {"--nd-mask", &Options::ndMask}, {"--groups", &Options::manifest}};
    const std::unordered_map<std::string, double Options::*> numbers{
        {"--deghost-ev", &Options::deghostEv}, {"--dark-model-error", &Options::darkModelError},
        {"--low", &Options::low}, {"--high", &Options::high}, {"--clip-margin", &Options::clipMargin},
        {"--nd-min-signal", &Options::ndMinSignal}, {"--floor", &Options::outputFloor}};
    bool positional = false, ndOption = false, ndUseOption = false;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto argument = [&] {
            if (++i >= argc || !*argv[i]) throw std::runtime_error("missing value for " + key);
            return std::string(argv[i]);
        };
        if (positional || key.empty() || key.front() != '-') o.plain.emplace_back(key);
        else if (key == "--") positional = true;
        else if (const auto it = strings.find(key); it != strings.end()) {
            o.*(it->second) = argument();
            ndUseOption |= key == "--nd-use";
        }
        else if (const auto it = paths.find(key); it != paths.end()) o.*(it->second) = argument();
        else if (const auto it = numbers.find(key); it != numbers.end()) {
            o.*(it->second) = number(argument());
            if (key == "--floor" && (o.outputFloor > std::numeric_limits<float>::max() ||
                                     !(float(o.outputFloor) > 0)))
                throw std::runtime_error("--floor requires a positive float32 value");
        }
        else if (key == "--nd") {
            const auto start = o.nd.size();
            while (i + 1 < argc && argv[i + 1][0] != '-') o.nd.emplace_back(argv[++i]);
            if (o.nd.size() == start) throw std::runtime_error("--nd requires at least one frame");
        } else if (key == "--nd-ev") {
            const auto value = argument();
            o.ndEv = value == "auto" ? std::nullopt : std::optional<double>(number(value));
            ndOption = true;
        } else if (key == "--reference") {
            const auto value = argument();
            o.reference = value == "auto" ? 0 : integer(value, 1, 65535);
        } else if (key == "--threads") o.threads = unsigned(integer(argument(), 1, 1024));
        else if (key == "--clip-dilate") o.clipDilate = integer(argument(), 0, 65535);
        else if (key == "--allow-nd-scene-drift") o.allowNDSceneDrift = true;
        else if (key == "--no-deghost") o.deghost = false;
        else if (key == "--aux") o.aux = true;
        else if (key == "--force") o.force = true;
        else if (key == "-v") o.verbose = true;
        else throw std::runtime_error("unknown option: " + key);
    }

    if (!o.manifest.empty()) {
        if (!o.plain.empty() || !o.nd.empty() || ndOption || !o.ndMask.empty() ||
            ndUseOption || o.allowNDSceneDrift)
            throw std::runtime_error("--groups takes input frames and ND settings from the manifest");
        o.groups = readGroups(o.manifest);
    }
    if ((o.plain.empty() && o.groups.empty()) || o.output.empty())
        throw std::runtime_error("provide RAW frames or --groups FILE, and -o OUT.exr (see --help)");
    const size_t baseFrames = o.groups.empty() ? o.plain.size() : o.groups.front().files.size();
    if (o.plain.size() + o.nd.size() > 65535 || size_t(o.reference) > baseFrames)
        throw std::runtime_error("frame count or reference is out of range");
    if (o.report.empty()) o.report = fs::path(o.output).replace_extension(".json");
    oneOf("--space", o.space, {"aces", "camera"});
    oneOf("--normalize", o.normalize, {"geomean", "none"});
    oneOf("--highlights", o.highlights, {"neutral", "bound"});
    oneOf("--demosaic", o.demosaic, {"bilinear", "mhc"});
    oneOf("--nd-use", o.ndUse, {"fill", "all"});
    if (o.matrix.rfind("custom:", 0) == 0) {
        const auto path = o.matrix.substr(7);
        std::ifstream file(path);
        if (!file) throw std::runtime_error("could not open custom matrix: " + path);
        Matrix matrix{};
        for (auto& row : matrix)
            for (auto& value : row)
                if (!(file >> value) || !std::isfinite(float(value)))
                    throw std::runtime_error("custom matrix requires nine finite float32 numbers in row order");
        std::string extra;
        if (file >> extra) throw std::runtime_error("custom matrix requires exactly nine numbers");
        o.customMatrix = matrix;
    } else oneOf("--matrix", o.matrix, {"metadata", "auto", "spectral", "adobe"});
    if (o.noise != "dng") {
        o.noiseOverride = numberList(o.noise, 2);
        if (o.noiseOverride[0] < 0 || o.noiseOverride[1] < 0 ||
            (o.noiseOverride[0] == 0 && o.noiseOverride[1] == 0))
            throw std::runtime_error("--noise coefficients must be nonnegative and include a positive term");
    }
    if (o.wb.rfind("frame:", 0) == 0) integer(o.wb.substr(6), 1, int(baseFrames));
    else if (o.wb.rfind("custom:", 0) == 0) {
        for (double w : numberList(o.wb.substr(7), 3))
            if (w <= 0 || w > 1e6 || w < 1e-6) throw std::runtime_error("invalid white balance multiplier");
    } else oneOf("--wb", o.wb, {"reference", "median"});
    if (o.low < 0 || o.low >= o.high || o.high > 1 ||
        o.clipMargin < 0 || o.clipMargin >= 1 || o.darkModelError < 0 || o.deghostEv < 0 ||
        o.ndMinSignal <= 0 || o.ndMinSignal > 1 || ((!o.nd.empty() || std::any_of(o.groups.begin(), o.groups.end(), [](const InputGroup& g) {
            return !g.calibrationReference.empty();
        })) && o.ndMinSignal > o.high) ||
        (o.ndEv && (*o.ndEv < 0 || *o.ndEv > 64)))
        throw std::runtime_error("signal thresholds, ND EV or uncertainty are out of range");
    if (o.nd.empty() && (ndOption || !o.ndMask.empty() || o.ndUse != "fill"))
        throw std::runtime_error("ND options require --nd frames");
    if (o.allowNDSceneDrift && (o.nd.empty() || !o.ndEv))
        throw std::runtime_error("--allow-nd-scene-drift requires a manual --nd-ev EV");
    return o;
}
}

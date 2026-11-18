#include "pipeline.h"

#include <OpenImageIO/paramlist.h>
#include <rawtoaces/image_converter.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace rawfuse {
namespace {

Matrix product(const Matrix& left, const Matrix& right) {
    Matrix result{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            for (int k = 0; k < 3; ++k) result[r][c] += left[r][k] * right[k][c];
    return result;
}

// rawtoaces writes informational lines and warnings straight to std::cerr.
// Collect them while it is configured so they follow --quiet and the
// fuser's warning list. configureColour runs on one thread.
class CaptureStderr {
public:
    CaptureStderr() : previous(std::cerr.rdbuf(text.rdbuf())) {}
    CaptureStderr(const CaptureStderr&) = delete;
    CaptureStderr& operator=(const CaptureStderr&) = delete;
    ~CaptureStderr() { std::cerr.rdbuf(previous); }
    std::vector<std::string> lines() const {
        std::vector<std::string> result;
        std::istringstream stream(text.str());
        for (std::string line; std::getline(stream, line);) {
            while (!line.empty() && (line.back() == ' ' || line.back() == '\r')) line.pop_back();
            if (!line.empty()) result.push_back(line);
        }
        return result;
    }

private:
    std::ostringstream text;
    std::streambuf* previous;
};

std::string afterPrefix(const std::string& line, const std::string& prefix) {
    return line.rfind(prefix, 0) == 0 ? line.substr(prefix.size()) : std::string();
}

bool sensible(const RGB& wb) {
    return std::all_of(wb.begin(), wb.end(), [](double w) { return w > 0 && std::isfinite(w); });
}

} // namespace

Colour configureColour(const std::vector<RawFrame>& plain, size_t reference,
                       const Options& o, Warnings& warnings) {
    using Settings = rta::util::ImageConverter::Settings;
    using Method = Settings::MatrixMethod;
    rta::util::ImageConverter converter;
    const auto& ref = plain[reference];
    const auto& planes = ref.metadata.planeColours;
    const bool reordered = planes != std::array<int, 3>{0, 1, 2};
    // LibRaw assumes green on plane 1 and can return an identity rgb_cam for
    // other DNG plane orders, so its Adobe matrix is not trusted for them.
    if (reordered && o.matrix == "adobe")
        throw std::runtime_error("Adobe matrix mode requires RGB DNG plane order; use --matrix metadata or custom:FILE");
    // The DNG solver consumes native camera-plane matrix columns. The fuser's
    // pixels and all caller-supplied multipliers/matrices use canonical RGB.
    bool nativePlanes = reordered && (o.matrix == "metadata" || o.matrix == "auto");
    Colour result;
    bool customWB = reordered;
    if (o.wb.rfind("custom:", 0) == 0) {
        const auto values = numberList(o.wb.substr(7), 3);
        std::copy(values.begin(), values.end(), result.wb.begin());
        customWB = true;
    } else if (o.wb.rfind("frame:", 0) == 0) {
        const size_t input = size_t(number(o.wb.substr(6))) - 1;
        const auto it = std::find_if(plain.begin(), plain.end(), [&](const auto& frame) {
            return frame.inputIndex == input;
        });
        if (!it->hasWB) throw std::runtime_error("selected white balance frame has no camera multipliers");
        result.wb = it->wb;
        customWB = true;
    } else if (o.wb == "median" && plain.size() > 1) {
        std::array<std::vector<double>, 3> values;
        for (const auto& frame : plain)
            if (frame.hasWB)
                for (int c = 0; c < 3; ++c) values[c].push_back(frame.wb[c]);
        if (values[0].empty()) warnings.push_back("first group has no camera white balance; using unit multipliers");
        else for (int c = 0; c < 3; ++c) result.wb[c] = median(values[c]);
        customWB = true;
    } else {
        result.wb = ref.wb;
        if (!ref.hasWB) {
            warnings.push_back("reference frame has no camera white balance; using unit multipliers");
            customWB = true;
        }
    }
    auto setWhiteBalance = [&] {
        converter.settings.WB_method = Settings::WBMethod::Custom;
        for (int i = 0; i < 4; ++i) {
            const int channel = i == 3 ? 1 : i;
            converter.settings.custom_WB[i] = float(result.wb[nativePlanes ? planes[channel] : channel]);
        }
    };
    if (customWB) setWhiteBalance();

    if (o.matrix.rfind("custom:", 0) == 0) {
        converter.settings.matrix_method = Method::Custom;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) converter.settings.custom_matrix[r][c] = float((*o.customMatrix)[r][c]);
    } else if (o.matrix == "adobe") converter.settings.matrix_method = Method::Adobe;
    else if (o.matrix == "spectral") converter.settings.matrix_method = Method::Spectral;
    else if (o.matrix == "auto") converter.settings.matrix_method = Method::Auto;
    else converter.settings.matrix_method = Method::Metadata;

    OIIO::ParamValueList hints;
    auto configure = [&] {
        bool configured;
        std::vector<std::string> lines;
        {
            CaptureStderr capture;
            configured = converter.configure(ref.path.string(), hints);
            lines = capture.lines();
        }
        for (const auto& line : lines) {
            const auto fallback = afterPrefix(line, "Info: Falling back to metadata matrix method");
            if (!fallback.empty() || line == "Info: Falling back to metadata matrix method") {
                const auto start = fallback.find_first_not_of(". ");
                auto reason = start == std::string::npos ? std::string() : fallback.substr(start);
                if (reason.rfind("because ", 0) == 0) reason.erase(0, 8);
                result.fallback = reason.empty() ? "no spectral data" : reason;
            } else if (afterPrefix(line, "Info: ").empty()) {  // Informational lines are dropped
                const auto warning = afterPrefix(line, "Warning: ");
                const auto text = "rawtoaces: " + (warning.empty() ? line : warning);
                // Each configure() repeats its warnings.
                if (std::find(warnings.begin(), warnings.end(), text) == warnings.end()) warnings.push_back(text);
            }
        }
        if (!configured)
            throw std::runtime_error(converter.last_error_message.empty() ?
                                     "could not configure rawtoaces" : converter.last_error_message);
    };
    try {
        configure();
        if (nativePlanes && o.matrix == "auto" && hints.contains("rta:FoundIlluminant")) {
            // Auto selected the RGB spectral solver rather than the native-plane
            // DNG solver. Configure it again with the caller's RGB white balance.
            nativePlanes = false;
            setWhiteBalance();
            hints.clear();
            configure();
        }
        if (o.matrix == "auto") result.method = hints.contains("rta:FoundIlluminant") ? "spectral" : "metadata";
        else result.method = o.matrix.rfind("custom:", 0) == 0 ? "custom" : o.matrix;
        if (result.method == "spectral" || o.matrix == "auto") result.databases = converter.settings.database_directories;
        if (result.method == "spectral") result.fallback.clear();
        const auto& wb = converter.get_WB_multipliers();
        if (wb.size() >= 3 && sensible({wb[0], wb[1], wb[2]}))
            for (int c = 0; c < 3; ++c) result.wb[nativePlanes ? planes[c] : c] = wb[c];

        const auto& matrix = converter.get_transform_matrix();
        Matrix transform{};
        if (matrix.size() != 3) throw std::runtime_error("rawtoaces did not return a 3 × 3 colour matrix");
        for (int r = 0; r < 3; ++r) {
            if (matrix[r].size() != 3) throw std::runtime_error("rawtoaces did not return a 3 × 3 colour matrix");
            for (int c = 0; c < 3; ++c) transform[r][nativePlanes ? planes[c] : c] = matrix[r][c];
        }
        const auto input = hints.find("raw:ColorSpace");
        if (input != hints.end() && input->get_string() == "XYZ") {
            // In this mode rawtoaces expects LibRaw's XYZ develop as its input.
            // Reproduce its linear camera-to-XYZ step without a uint16 develop.
            constexpr Matrix srgbToXYZ{{{0.4124564, 0.3575761, 0.1804375},
                                        {0.2126729, 0.7151522, 0.0721750},
                                        {0.0193339, 0.1191920, 0.9503041}}};
            transform = product(transform, product(srgbToXYZ, ref.cameraToSRGB));
        }
        for (const auto& row : transform)
            for (double value : row)
                if (!std::isfinite(value)) throw std::runtime_error("rawtoaces returned a non-finite colour matrix");
        result.matrix = transform;
    } catch (const std::exception& error) {
        if (o.space == "aces" || (o.matrix != "metadata" && o.matrix != "auto"))
            throw std::runtime_error("rawtoaces colour configuration: " + std::string(error.what()));
        warnings.push_back("rawtoaces matrix unavailable in camera-space output: " + std::string(error.what()));
    }
    const double minimum = *std::min_element(result.wb.begin(), result.wb.end());
    for (double& w : result.wb) w /= minimum;
    return result;
}
}

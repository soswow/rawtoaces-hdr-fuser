#include "pipeline.h"

#include <libraw/libraw.h>

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace rawfuse {
namespace {

void readLevels(RawFrame& frame, LibRaw& raw, Warnings& warnings) {
    const auto& meta = frame.metadata;
    const auto& levels = raw.imgdata.color;
    frame.white = float(levels.maximum);
    if (!meta.white.empty()) {
        if (meta.white.size() != 1) throw std::runtime_error("expected a single Bayer WhiteLevel");
        frame.white = float(meta.white.front());
        if (levels.maximum && std::abs(double(levels.maximum) - frame.white) > 1)
            warnings.push_back("DNG WhiteLevel differs from LibRaw; using the tag");
    }
    if (!meta.black.empty()) {
        frame.blackRows = meta.blackRepeatRows;
        frame.blackCols = meta.blackRepeatCols;
        if (meta.black.size() != size_t(frame.blackRows) * frame.blackCols)
            throw std::runtime_error("BlackLevel does not match BlackLevelRepeatDim");
        frame.blackPattern.assign(meta.black.begin(), meta.black.end());
        frame.blackDeltaRows.assign(meta.blackRows.begin(), meta.blackRows.end());
        frame.blackDeltaCols.assign(meta.blackCols.begin(), meta.blackCols.end());
        if (meta.black.size() == 1 && frame.blackDeltaCols.empty() && frame.blackDeltaRows.empty() &&
            std::abs(double(levels.black) + levels.cblack[0] + levels.cblack[6] - meta.black[0]) > 1)
            warnings.push_back("DNG BlackLevel differs from LibRaw; using the tag");
    } else {
        const unsigned rows = std::max(1u, levels.cblack[4]);
        const unsigned cols = std::max(1u, levels.cblack[5]);
        if (uint64_t(rows) * cols > LIBRAW_CBLACK_SIZE - 6)
            throw std::runtime_error("invalid LibRaw black-level pattern");
        auto period = [](unsigned n) { return n % 2 ? 2 * n : n; };
        frame.blackRows = int(period(rows));
        frame.blackCols = int(period(cols));
        for (int row = 0; row < frame.blackRows; ++row)
            for (int col = 0; col < frame.blackCols; ++col)
                frame.blackPattern.push_back(float(double(levels.black) + levels.cblack[raw.COLOR(row, col)] +
                                                  levels.cblack[6 + (row % rows) * cols + col % cols]));
        if (std::all_of(frame.blackPattern.begin(), frame.blackPattern.end(),
                        [&](float value) { return value == frame.blackPattern[0]; })) {
            frame.blackRows = frame.blackCols = 1;
            frame.blackPattern.resize(1);
        }
    }
    auto checkDelta = [](std::vector<float>& data, int count) {
        if (data.empty()) return 0.0f;
        if (data.size() != size_t(count)) throw std::runtime_error("black-level deltas do not match active area");
        return *std::max_element(data.begin(), data.end());
    };
    const float maxBlack = *std::max_element(frame.blackPattern.begin(), frame.blackPattern.end()) +
                           checkDelta(frame.blackDeltaRows, frame.layout.height) +
                           checkDelta(frame.blackDeltaCols, frame.layout.width);
    if (!std::isfinite(maxBlack) || !std::isfinite(frame.white) || frame.white <= maxBlack ||
        frame.white > 65535)
        throw std::runtime_error("invalid black or white level");
}

DarkModel readDarkModel(const fs::path& path, const Layout& size) {
    DarkModel result;
    if (path.empty()) return result;
    std::ifstream file(path);
    auto fields = [](std::string line) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::vector<std::string> columns;
        std::istringstream values(line);
        for (std::string value; std::getline(values, value, ',');) columns.push_back(value);
        return columns;
    };
    std::string line;
    if (!std::getline(file, line)) throw std::runtime_error("could not read dark model: " + path.string());
    const auto header = fields(line);
    auto column = [&](const std::string& name) {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end()) throw std::runtime_error("dark model needs column " + name);
        return size_t(it - header.begin());
    };
    const auto rowCol = column("row"), colCol = column("col");
    const auto offsetCol = column("offset_dn"), rateCol = column("rate_dn_per_s");
    const auto required = std::max({rowCol, colCol, offsetCol, rateCol});
    size_t lineNumber = 1;
    while (std::getline(file, line)) {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r") == line.npos) continue;
        try {
            const auto row = fields(line);
            if (row.size() <= required) throw std::runtime_error("missing field");
            const double y = number(row[rowCol]), x = number(row[colCol]);
            const float offset = float(number(row[offsetCol])), rate = float(number(row[rateCol]));
            if (y < 0 || y >= size.height || x < 0 || x >= size.width ||
                y != std::floor(y) || x != std::floor(x) || !std::isfinite(offset) || !std::isfinite(rate))
                throw std::runtime_error("out-of-range coordinate or coefficient");
            if (!result.emplace(size.index(int(y), int(x)), std::array<float, 2>{offset, rate}).second)
                throw std::runtime_error("duplicate coordinate");
        } catch (const std::exception& error) {
            throw std::runtime_error("dark model line " + std::to_string(lineNumber) + ": " + error.what());
        }
    }
    if (file.bad()) throw std::runtime_error("failed reading dark model");
    return result;
}

std::vector<BadPixel> readBadPixels(const std::string& path, const Layout& layout) {
    if (path == "dng" || path == "none") return {};
    std::ifstream file(path);
    if (!file) throw std::runtime_error("could not open bad-pixel list: " + path);
    std::vector<BadPixel> result;
    for (std::string line; std::getline(file, line);) {
        line.resize(std::min(line.size(), line.find('#')));
        if (line.find_first_not_of(" \t\r") == line.npos) continue;
        std::istringstream values(line);
        int col, row;
        long long timestamp = 0;
        if (!(values >> col >> row)) throw std::runtime_error("invalid bad-pixel list row");
        values >> std::ws;
        if (!values.eof() && !(values >> timestamp)) throw std::runtime_error("invalid bad-pixel timestamp");
        if (row < 0 || col < 0 || row >= layout.height || col >= layout.width)
            throw std::runtime_error("bad-pixel list coordinate is outside the active image");
        values >> std::ws;
        if (!values.eof() || timestamp < 0) throw std::runtime_error("invalid bad-pixel timestamp");
        result.push_back({{row, col}, timestamp});
    }
    if (file.bad()) throw std::runtime_error("failed reading bad-pixel list");
    return result;
}

void markClipping(RawFrame& frame, const Options& o, std::vector<uint8_t>& horizontal, const Mask& validity) {
    const auto& size = frame.layout;
    const int radius = o.clipDilate;
    frame.rejected.resize(size.size());
    if (radius) horizontal.resize(size.size());
    parallelRows(size.height, o.threads, [&](int row) {
        const size_t start = size.index(row, 0);
        auto* mask = frame.rejected.data() + start;
        const auto* data = frame.mosaic + size_t(row) * frame.stride;
        for (int col = 0; col < size.width; ++col) {
            const float black = frame.black(row, col);
            const bool eligible = validity.empty() || validity[start + col];
            mask[col] = eligible && data[col] >= frame.white - o.clipMargin * (frame.white - black) ? saturated : 0;
        }
        if (!radius) return;
        int hits = 0;
        for (int col = 0; col < std::min(radius, size.width); ++col) hits += mask[col];
        for (int col = 0; col < size.width; ++col) {
            if (col + radius < size.width) hits += mask[col + radius];
            horizontal[start + col] = hits != 0;
            if (col >= radius) hits -= mask[col - radius];
        }
    });
    if (!radius) return;
    const int bands = (size.height + 31) / 32;
    parallelRows(bands, o.threads, [&](int band) {
        const int first = band * 32, end = std::min(first + 32, size.height);
        std::vector<int> hits(size.width);
        for (int row = std::max(first - radius, 0); row < std::min(first + radius, size.height); ++row) {
            const auto* enter = horizontal.data() + size.index(row, 0);
            for (int col = 0; col < size.width; ++col) hits[col] += enter[col];
        }
        for (int row = first; row < end; ++row) {
            const auto* enter = row + radius < size.height ? horizontal.data() + size.index(row + radius, 0) : nullptr;
            const auto* leave = row >= radius ? horizontal.data() + size.index(row - radius, 0) : nullptr;
            auto* mask = frame.rejected.data() + size.index(row, 0);
            for (int col = 0; col < size.width; ++col) {
                if (enter) hits[col] += enter[col];
                if (hits[col] && !(mask[col] & saturated)) mask[col] = bloom;
                if (leave) hits[col] -= leave[col];
            }
        }
    });
}

void markDefects(RawFrame& frame, const Options& options, const std::vector<BadPixel>& sidecar,
                 const DarkModel& model) {
    const auto& layout = frame.layout;
    auto exclude = [&](size_t i) {
        if (!(frame.rejected[i] & excluded)) ++frame.hotMissing;
        frame.rejected[i] |= excluded;
    };
    // Without a capture time, a dated entry cannot be ruled out; excluding a
    // good site costs an interpolation, keeping a defect leaves an artefact.
    for (auto p : sidecar) {
        if (p.since && !frame.captureTimestamp) ++frame.undatedDefects;
        if (p.since <= frame.captureTimestamp || !frame.captureTimestamp)
            exclude(layout.index(p.site.row, p.site.col));
    }
    for (const auto& entry : model) exclude(entry.first);
    if (options.badPixels != "dng") return;
    for (auto p : frame.metadata.badPoints) exclude(layout.index(p.row, p.col));
    for (auto r : frame.metadata.badRects) {
        for (int y = r.top; y < r.bottom; ++y)
            for (int x = r.left; x < r.right; ++x) {
                const size_t i = layout.index(y, x);
                exclude(i);
                frame.rejected[i] |= badRect;
            }
    }
}

// Robust spread of LibRaw's masked (optical black) areas outside the active grid,
// measured around each Bayer phase's own median. Returns 0 without enough pixels.
double opticalBlackSigma(const LibRaw& raw, const RawFrame& frame) {
    const auto& sizes = raw.imgdata.sizes;
    const auto* data = raw.imgdata.rawdata.raw_image;
    std::array<std::vector<double>, 4> phases;
    constexpr size_t cap = 1'000'000;
    for (const auto& area : sizes.mask) {
        // Skip a few pixels at each edge, where masked areas often ramp.
        const int top = area[0] + 4, left = area[1] + 4, bottom = area[2] - 4, right = area[3] - 4;
        if (top < 0 || left < 0 || top >= bottom || left >= right ||
            bottom > frame.rawHeight || right > frame.rawWidth)
            continue;
        if (top < frame.top + frame.layout.height && bottom > frame.top &&
            left < frame.left + frame.layout.width && right > frame.left)
            continue;
        for (int y = top; y < bottom; ++y)
            for (int x = left; x < right; ++x) {
                auto& values = phases[(y & 1) * 2 + (x & 1)];
                if (values.size() < cap) values.push_back(data[size_t(y) * frame.stride + x]);
            }
    }
    std::vector<double> deviations;
    for (auto& values : phases) {
        if (values.empty()) continue;
        auto copy = values;
        const double centre = median(copy);
        for (double value : values) deviations.push_back(value - centre);
    }
    if (deviations.size() < 10000) return 0;
    return std::sqrt(robustVariance(std::move(deviations)));
}

} // namespace

std::string Layout::pattern() const {
    std::string result;
    for (const auto c : cfa) result += "RGB"[c];
    return result;
}

RawFrame readRaw(const fs::path& path, Warnings& warnings) {
    RawFrame frame;
    frame.path = path;
    auto raw = std::make_shared<LibRaw>();
    auto check = [&](int code) {
        if (code != LIBRAW_SUCCESS) throw std::runtime_error(libraw_strerror(code));
    };
    try {
        check(raw->open_file(path.string().c_str()));
        const auto& info = raw->imgdata.idata;
        if (info.is_foveon || info.colors != 3 || info.filters < 1000)
            throw std::runtime_error("expected a 2 × 2 Bayer mosaic");
        check(raw->unpack());
        const auto& sizes = raw->imgdata.sizes;
        if (!raw->imgdata.rawdata.raw_image || sizes.raw_pitch % 2 ||
            sizes.raw_pitch / 2 < sizes.raw_width)
            throw std::runtime_error("expected a uint16 Bayer mosaic with a valid row stride");
        frame.layout.width = sizes.width;
        frame.layout.height = sizes.height;
        frame.top = sizes.top_margin;
        frame.left = sizes.left_margin;
        frame.rawWidth = sizes.raw_width;
        frame.rawHeight = sizes.raw_height;
        frame.stride = sizes.raw_pitch / 2;
        frame.dng = info.dng_version != 0;
        Warnings local;
        if (frame.dng) {
            frame.metadata = readDngMetadata(path, local);
            // LibRaw may round odd DNG margins up for its develop path. We consume
            // raw_image directly, so retain the exact grid used by DNG calibration
            // tags and caller masks. CFAPattern is relative to ActiveArea.
            const auto area = frame.metadata.activeArea.value_or(Rect{0, 0, frame.rawHeight, frame.rawWidth});
            frame.top = area.top;
            frame.left = area.left;
            frame.layout.height = area.bottom - area.top;
            frame.layout.width = area.right - area.left;
        }
        if (frame.layout.width < 2 || frame.layout.height < 2 ||
            frame.left + frame.layout.width > frame.rawWidth || frame.top + frame.layout.height > frame.rawHeight)
            throw std::runtime_error("invalid Bayer active area");
        if (frame.metadata.cfa) {
            frame.layout.cfa = *frame.metadata.cfa;
            // Cross-check the exact-grid DNG parser against LibRaw's decoded pattern.
            // LibRaw's pattern is reliable only for RGB plane order: it assumes
            // green on plane 1 and relabels a duplicate plane-1 site as plane 3.
            if (frame.metadata.planeColours == std::array<int, 3>{0, 1, 2}) {
                const int dy = frame.top - sizes.top_margin, dx = frame.left - sizes.left_margin;
                bool agrees = true;
                for (int row = 0; row < 2; ++row)
                    for (int col = 0; col < 2; ++col) {
                        int plane = raw->COLOR(((row + dy) % 2 + 2) % 2, ((col + dx) % 2 + 2) % 2);
                        if (plane == 3) plane = 1;
                        agrees &= plane == frame.layout.cfa[row * 2 + col];
                    }
                if (!agrees)
                    local.push_back("DNG CFAPattern disagrees with LibRaw's decoded Bayer pattern; using the DNG tags");
            }
        } else for (int row = 0; row < 8; ++row) {
            for (int col = 0; col < 2; ++col) {
                const char channel = info.cdesc[raw->COLOR(row, col)];
                const int colour = channel == 'R' ? 0 : channel == 'G' ? 1 : channel == 'B' ? 2 : -1;
                if (colour < 0) throw std::runtime_error("unsupported Bayer colour");
                auto& entry = frame.layout.cfa[(row & 1) * 2 + col];
                if (row < 2) entry = colour;
                else if (entry != colour) throw std::runtime_error("CFA does not repeat every 2 rows");
            }
        }
        const auto pattern = frame.layout.pattern();
        if (pattern != "RGGB" && pattern != "BGGR" && pattern != "GRBG" && pattern != "GBRG")
            throw std::runtime_error("expected a standard RGGB, BGGR, GRBG or GBRG Bayer mosaic");

        frame.make = info.make;
        frame.model = info.model;
        auto& m = frame.metadata;
        const auto& other = raw->imgdata.other;
        frame.captureTimestamp = std::max<long long>(0, other.timestamp);
        if (m.seconds <= 0) m.seconds = other.shutter;
        if (m.seconds <= 0 || !std::isfinite(m.seconds)) throw std::runtime_error("missing exposure time");
        if (m.exposure.empty()) {
            std::ostringstream value;
            value << std::setprecision(9) << m.seconds;
            m.exposure = value.str();
        }
        if (m.iso <= 0 || m.iso == 65535) m.iso = other.iso_speed == 65535 ? 0 : other.iso_speed;
        if (m.aperture <= 0) m.aperture = other.aperture;
        if (m.description.empty()) m.description = other.desc;
        if (m.captureDate.empty() && other.timestamp > 0) {
            // LibRaw builds the timestamp from the camera's local wall-clock
            // date with mktime(), so localtime() reproduces the recorded date.
            std::tm date{};
            if (localtime_r(&other.timestamp, &date)) {
                std::ostringstream text;
                text << std::put_time(&date, "%Y:%m:%d %H:%M:%S");
                m.captureDate = text.str();
            }
        }
        const auto& c = raw->imgdata.color;
        for (int i = 0; i < 3; ++i) {
            const int colour = m.planeColours[i];
            frame.wb[colour] = m.neutral ? 1 / (*m.neutral)[colour] : c.cam_mul[i];
        }
        frame.hasWB = std::all_of(frame.wb.begin(), frame.wb.end(), [](double w) {
            return std::isfinite(w) && w > 0;
        });
        if (frame.hasWB) {
            const double green = frame.wb[1];
            for (auto& w : frame.wb) w /= green;
        } else {
            if (m.neutral) local.push_back("invalid AsShotNeutral");
            frame.wb = unitWB;
        }
        // List 1 DNG opcodes use raw sensor coordinates. Use active-area coordinates
        // internally for all defect and calibration inputs.
        for (auto& p : m.badPoints) {
            p.row -= frame.top;
            p.col -= frame.left;
        }
        m.badPoints.erase(std::remove_if(m.badPoints.begin(), m.badPoints.end(), [&](Point p) {
            return p.row < 0 || p.col < 0 || p.row >= frame.layout.height || p.col >= frame.layout.width;
        }), m.badPoints.end());
        for (auto& rect : m.badRects) {
            rect.top = std::max(0, rect.top - frame.top);
            rect.left = std::max(0, rect.left - frame.left);
            rect.bottom = std::min(frame.layout.height, rect.bottom - frame.top);
            rect.right = std::min(frame.layout.width, rect.right - frame.left);
        }
        readLevels(frame, *raw, local);
        // LibRaw's rgb_cam maps white-balanced camera values to linear sRGB.
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                frame.cameraToSRGB[i][m.planeColours[j]] = c.rgb_cam[i][j];

        frame.opticalBlackSigma = opticalBlackSigma(*raw, frame);
        frame.mosaic = raw->imgdata.rawdata.raw_image + size_t(frame.top) * frame.stride + frame.left;
        frame.mosaicOwner = std::move(raw);
        for (const auto& warning : local) warnings.push_back(path.filename().string() + ": " + warning);
        return frame;
    } catch (const std::exception& error) {
        throw std::runtime_error(path.string() + ": " + error.what());
    }
}

double validateFrames(const FrameGroups& groups, Warnings& warnings) {
    const auto& first = groups.front()->front();
    auto same = [&](const RawFrame& frame) {
        return frame.layout.width == first.layout.width && frame.layout.height == first.layout.height &&
               frame.top == first.top && frame.left == first.left && frame.rawWidth == first.rawWidth &&
               frame.rawHeight == first.rawHeight && frame.layout.cfa == first.layout.cfa &&
               frame.make == first.make && frame.model == first.model;
    };
    // Black and white levels may legitimately differ per frame: LibRaw measures
    // some cameras' black level from each frame's masked pixels. Every sample
    // uses its own frame's levels; only report differences beyond small jitter.
    auto meanBlack = [](const RawFrame& frame) {
        double sum = 0;
        for (float value : frame.blackPattern) sum += value;
        return sum / double(frame.blackPattern.size());
    };
    const double firstBlack = meanBlack(first);
    double blackSpread = 0, whiteSpread = 0;
    double aperture = 1;
    for (const auto& frame : *groups.front())
        if (std::isfinite(frame.metadata.aperture) && frame.metadata.aperture > 0) {
            aperture = frame.metadata.aperture;
            break;
        }
    const double firstISO = first.metadata.iso > 0 ? first.metadata.iso : 100;
    bool isoDiffers = false;
    size_t baselineFrames = 0;
    for (auto* group : groups) {
        for (auto& frame : *group) {
            auto& meta = frame.metadata;
            if (!same(frame))
                throw std::runtime_error(frame.path.string() + ": camera, CFA or active area differ");
            blackSpread = std::max(blackSpread, std::abs(meanBlack(frame) - firstBlack));
            whiteSpread = std::max(whiteSpread, double(std::abs(frame.white - first.white)));
            if (!std::isfinite(meta.iso) || meta.iso <= 0) {
                warnings.push_back(frame.path.filename().string() + ": no ISO; assuming 100");
                meta.iso = 100;
            }
            if (!std::isfinite(meta.aperture) || meta.aperture <= 0) {
                warnings.push_back(frame.path.filename().string() + ": no f-number; assuming the reference aperture");
                meta.aperture = aperture;
            }
            const double ratio = aperture / meta.aperture;
            frame.effective = meta.seconds * (meta.iso / 100) * ratio * ratio;
            if (!(frame.effective > 0) || !std::isfinite(frame.effective))
                throw std::runtime_error(frame.path.string() + ": invalid effective exposure");
            frame.exifEffective = frame.effective;
            isoDiffers |= meta.iso != firstISO;
            baselineFrames += meta.baseline != 0;
        }
    }
    const double tolerance = 0.001 * (first.white - firstBlack);
    if (blackSpread > tolerance || whiteSpread > tolerance) {
        std::ostringstream message;
        message << "black or white levels differ between frames (up to " << blackSpread << " / "
                << whiteSpread << " DN); using each frame's own levels";
        warnings.push_back(message.str());
    }
    if (baselineFrames)
        warnings.push_back("BaselineExposure is recorded without applying it (" +
                           std::to_string(baselineFrames) + " input frames)");
    if (isoDiffers) warnings.push_back("ISO differs between frames; using each frame's reported ISO and noise");
    return aperture;
}

Sidecars loadSidecars(const Options& o, const Layout& layout) {
    Sidecars result;
    result.darkModel = readDarkModel(o.darkModel, layout);
    result.badPixels = readBadPixels(o.badPixels, layout);
    for (const auto& group : o.groups) {
        result.valid.push_back(readMask(group.validMask, layout));
        result.contribution.push_back(readMask(group.contributionMask, layout));
        result.calibration.push_back(readMask(group.calibrationMask, layout));
        result.drift.push_back(readMask(group.driftMask, layout));
    }
    if (o.vignette != "dng" && o.vignette != "none") {
        const int channels = imageChannels(o.vignette, layout);
        if (channels != 1 && channels != 3)
            throw std::runtime_error("flat-field EXR must have one or three channels");
    }
    return result;
}

void releaseMosaics(std::vector<RawFrame>& frames) {
    for (auto& frame : frames) {
        frame.mosaic = nullptr;
        frame.mosaicOwner.reset();
        std::vector<uint8_t>().swap(frame.rejected);
        std::vector<std::pair<size_t, float>>().swap(frame.dark);
    }
}

Preparation prepareFrames(const FrameGroups& groups, const Sidecars& sidecars,
                          const Options& o, Warnings& warnings) {
    const auto& model = sidecars.darkModel;
    const auto& bad = sidecars.badPixels;
    std::vector<uint8_t> horizontal;
    size_t undated = 0;
    Preparation result;
    result.darkFits.resize(groups.size());
    result.exposureSteps.resize(groups.size());
    for (size_t g = 0; g < groups.size(); ++g) {
        const Mask empty;
        const auto& validity = sidecars.valid.empty() ? empty : sidecars.valid[g];
        for (auto& frame : *groups[g]) {
            markClipping(frame, o, horizontal, validity);
            markDefects(frame, o, bad, model);
            undated += frame.undatedDefects != 0;
            for (size_t i = 0; i < validity.size(); ++i)
                if (!validity[i]) frame.rejected[i] |= outsideValidity;
        }
    }
    // Exposure ratios first: the noise fit compares frames at their true ratios.
    auto byExposure = [](const RawFrame& a, const RawFrame& b) { return a.effective < b.effective; };
    result.exposureSteps = measureExposures(groups, o, warnings);
    for (auto* group : groups) std::stable_sort(group->begin(), group->end(), byExposure);
    result.noiseFits = assignNoise(groups, o, warnings);
    for (size_t g = 0; g < groups.size(); ++g) {
        if (o.darkModel.empty() || groups[g]->empty()) continue;
        result.darkFits[g] = fitDarkModel(*groups[g], model, o);
        if (result.darkFits[g].method == "unavailable")
            warnings.push_back(o.groups[g].name +
                " sweep has too few usable dark-model sites to fit; treating model sites as missing");
    }
    if (undated)
        warnings.push_back(std::to_string(undated) + " input frames have no capture time; "
                           "excluding every dated bad-pixel list entry");
    return result;
}

Sample sample(const RawFrame& frame, int row, int col, size_t index, const Options& o) {
    const auto mask = frame.rejected[index];
    if (mask & (excluded | outsideValidity)) return {SampleState::missing};
    if (mask & saturated) return {SampleState::clipped};
    const float black = frame.black(row, col);
    const double full = frame.white - black;
    double dark = 0, extra = 0;
    if (mask & modeled) {
        const auto it = std::lower_bound(frame.dark.begin(), frame.dark.end(), index,
            [](const auto& entry, size_t site) { return entry.first < site; });
        dark = it->second / full;
        const double sigma = frame.darkUncertainty * dark;
        extra = frame.noise[frame.layout.colour(row, col)].shot * std::max(0.0, dark) + sigma * sigma;
    }
    const float fraction = float((frame.raw(row, col) - black) / full - dark);
    if (mask & bloom) return {SampleState::adjacentClip, fraction, float(extra)};
    return {fraction <= o.high ? SampleState::usable : SampleState::aboveHigh, fraction, float(extra)};
}
}

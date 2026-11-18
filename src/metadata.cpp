#include "pipeline.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace rawfuse {
namespace {

using Bytes = std::vector<uint8_t>;

uint64_t unsignedValue(const Bytes& data, size_t offset, size_t length, bool little) {
    if (offset > data.size() || length > data.size() - offset)
        throw std::runtime_error("truncated TIFF or DNG opcode");
    uint64_t result = 0;
    for (size_t i = 0; i < length; ++i)
        result = (result << 8) | data[offset + (little ? length - 1 - i : i)];
    return result;
}

double floating(const Bytes& data, size_t offset, size_t length, bool little) {
    uint64_t bits = unsignedValue(data, offset, length, little);
    if (length == 4) {
        const uint32_t word = uint32_t(bits);
        float result;
        std::memcpy(&result, &word, sizeof(result));
        return result;
    }
    double result;
    std::memcpy(&result, &bits, sizeof(result));
    return result;
}

struct Entry { unsigned type; size_t count; uint64_t offset; };
using Directory = std::unordered_map<unsigned, Entry>;

class Tiff {
public:
    explicit Tiff(const fs::path& path) : in(path, std::ios::binary | std::ios::ate) {
        if (!in || in.tellg() < 0) throw std::runtime_error("could not open TIFF metadata");
        size = uint64_t(in.tellg());
        const auto header = read(0, 8);
        little = header[0] == 'I';
        if (!(little && header[1] == 'I') && !(header[0] == 'M' && header[1] == 'M'))
            throw std::runtime_error("invalid TIFF byte order");
        if (u(header, 2, 2) != 42) throw std::runtime_error("expected classic TIFF in DNG");
        first = u(header, 4, 4);
    }

    Directory directory(uint64_t offset, uint64_t& next) {
        const size_t n = u(read(offset, 2), 0, 2);
        const auto table = read(offset + 2, n * 12 + 4);
        Directory dir;
        for (size_t i = 0; i < n * 12; i += 12) {
            const unsigned type = unsigned(u(table, i + 2, 2));
            const size_t count = u(table, i + 4, 4);
            const size_t unit = typeSize(type);
            if (!unit) continue;
            const uint64_t location = count <= 4 / unit ? offset + 2 + i + 8 : u(table, i + 8, 4);
            dir.emplace(unsigned(u(table, i, 2)), Entry{type, count, location});
        }
        next = u(table, n * 12, 4);
        return dir;
    }

    Bytes bytes(const Directory& dir, unsigned tag) {
        const auto it = dir.find(tag);
        if (it == dir.end()) return {};
        return read(it->second.offset, it->second.count * typeSize(it->second.type));
    }

    std::string string(const Directory& dir, unsigned tag, bool trim = true) {
        auto data = bytes(dir, tag);
        if (trim) while (!data.empty() && !data.back()) data.pop_back();
        return {data.begin(), data.end()};
    }

    std::vector<double> numbers(const Directory& dir, unsigned tag) {
        const auto it = dir.find(tag);
        if (it == dir.end()) return {};
        const auto data = bytes(dir, tag);
        std::vector<double> values;
        values.reserve(it->second.count);
        const unsigned type = it->second.type;
        const size_t width = typeSize(type);
        for (size_t offset = 0; offset < data.size(); offset += width) {
            double value;
            if (type == 11 || type == 12) value = floating(data, offset, width, little);
            else if (type == 5 || type == 10) {
                const auto ratio = rational(data, offset, type);
                if (!ratio.second) throw std::runtime_error("zero denominator in TIFF tag " + std::to_string(tag));
                value = double(ratio.first) / double(ratio.second);
            } else {
                const auto bits = u(data, offset, width);
                if (type == 6) value = int8_t(bits);
                else if (type == 8) value = int16_t(bits);
                else if (type == 9) value = int32_t(bits);
                else value = double(bits);
            }
            if (!std::isfinite(value))
                throw std::runtime_error("non-finite TIFF tag " + std::to_string(tag));
            values.push_back(value);
        }
        return values;
    }

    double scalar(const Directory& dir, unsigned tag, double fallback = 0) {
        auto values = numbers(dir, tag);
        return values.empty() ? fallback : values.front();
    }

    std::string exposure(const Directory& dir) {
        const auto it = dir.find(33434);
        if (it == dir.end() || it->second.count != 1 ||
            (it->second.type != 5 && it->second.type != 10)) return {};
        const auto ratio = rational(bytes(dir, 33434), 0, it->second.type);
        return std::to_string(ratio.first) + '/' + std::to_string(ratio.second);
    }

    uint64_t first;

private:
    std::ifstream in;
    uint64_t size;
    bool little;

    static size_t typeSize(unsigned type) {
        constexpr std::array<size_t, 14> sizes{0, 1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8, 4};
        return type < sizes.size() ? sizes[type] : 0;
    }

    Bytes read(uint64_t offset, size_t amount) {
        constexpr size_t maxMetadataBytes = 64 * 1024 * 1024;
        if (offset > size || amount > size - offset || amount > maxMetadataBytes)
            throw std::runtime_error("TIFF metadata extends outside the file or exceeds 64 MiB");
        Bytes result(amount);
        if (amount) {
            in.seekg(std::streamoff(offset));
            if (!in.read(reinterpret_cast<char*>(result.data()), std::streamsize(amount)))
                throw std::runtime_error("could not read TIFF metadata");
        }
        return result;
    }

    uint64_t u(const Bytes& data, size_t offset, size_t width) const {
        return unsignedValue(data, offset, width, little);
    }

    std::pair<int64_t, int64_t> rational(const Bytes& data, size_t offset, unsigned type) const {
        const auto a = u(data, offset, 4);
        const auto b = u(data, offset + 4, 4);
        if (type == 10) return {int32_t(a), int32_t(b)};
        return {int64_t(a), int64_t(b)};
    }
};

class Opcode {
public:
    explicit Opcode(const Bytes& data) : data(data) {}
    size_t remaining() const { return data.size() - offset; }

    uint32_t u32() {
        const auto value = unsignedValue(data, offset, 4, false);
        offset += 4;
        return uint32_t(value);
    }

    int coordinate() {
        const auto value = u32();
        if (value > uint32_t(std::numeric_limits<int>::max()))
            throw std::runtime_error("out-of-range DNG opcode coordinate");
        return int(value);
    }

    double real(size_t width) {
        const double value = floating(data, offset, width, false);
        offset += width;
        if (!std::isfinite(value)) throw std::runtime_error("non-finite DNG opcode value");
        return value;
    }

    Bytes take(size_t length) {
        if (length > remaining()) throw std::runtime_error("truncated DNG opcode parameters");
        Bytes result(data.begin() + offset, data.begin() + offset + length);
        offset += length;
        return result;
    }

    Rect rect() {
        Rect result{coordinate(), coordinate(), coordinate(), coordinate()};
        if (result.top >= result.bottom || result.left >= result.right)
            throw std::runtime_error("empty DNG opcode region");
        return result;
    }

private:
    const Bytes& data;
    size_t offset = 0;
};

void parseOpcodes(Metadata& result, const Bytes& bytes, int list, Warnings& warnings) {
    if (bytes.empty()) return;
    Opcode in(bytes);
    const auto count = in.u32();
    for (uint32_t i = 0; i < count; ++i) {
        const auto id = in.u32();
        in.u32(); // Minimum DNG version
        const auto flags = in.u32();
        const auto params = in.take(in.u32());
        Opcode op(params);
        if (id == 5 && list == 1) {
            if (op.u32() > 3) throw std::runtime_error("invalid FixBadPixelsList Bayer phase");
            const auto points = op.u32();
            const auto rects = op.u32();
            if (uint64_t(points) * 8 + uint64_t(rects) * 16 != op.remaining())
                throw std::runtime_error("invalid FixBadPixelsList size");
            for (uint32_t j = 0; j < points; ++j)
                result.badPoints.push_back({op.coordinate(), op.coordinate()});
            for (uint32_t j = 0; j < rects; ++j)
                result.badRects.push_back(op.rect());
        } else if (id == 9 && list == 2) {
            GainMap map{};
            map.area = op.rect();
            const auto plane = op.u32();
            const auto planes = op.u32();
            map.rowPitch = op.coordinate();
            map.colPitch = op.coordinate();
            map.rows = op.coordinate();
            map.cols = op.coordinate();
            map.spacingV = op.real(8);
            map.spacingH = op.real(8);
            map.originV = op.real(8);
            map.originH = op.real(8);
            const auto mapPlanes = op.u32();
            if (plane != 0 || planes != 1 || mapPlanes != 1)
                throw std::runtime_error("GainMap must address the single plane of a Bayer mosaic");
            if (!map.rowPitch || !map.colPitch || !map.rows || !map.cols ||
                map.spacingV <= 0 || map.spacingH <= 0 ||
                uint64_t(map.rows) * map.cols * 4 != op.remaining())
                throw std::runtime_error("invalid GainMap grid");
            map.gains.reserve(size_t(map.rows) * map.cols);
            while (op.remaining()) {
                const float gain = float(op.real(4));
                if (gain <= 0) throw std::runtime_error("nonpositive GainMap value");
                map.gains.push_back(gain);
            }
            result.maps.push_back(std::move(map));
        } else {
            warnings.push_back("ignored DNG opcode " + std::to_string(id) + " in list " +
                               std::to_string(list) + ((flags & 1) ? " (optional)" : ""));
        }
    }
    if (in.remaining()) throw std::runtime_error("trailing bytes in DNG opcode list");
}

} // namespace

bool GainMap::operator==(const GainMap& other) const {
    return std::tie(area.top, area.left, area.bottom, area.right, rowPitch, colPitch, rows, cols,
                    spacingV, spacingH, originV, originH, gains) ==
           std::tie(other.area.top, other.area.left, other.area.bottom, other.area.right,
                    other.rowPitch, other.colPitch, other.rows, other.cols,
                    other.spacingV, other.spacingH, other.originV, other.originH, other.gains);
}

Metadata readDngMetadata(const fs::path& path, Warnings& warnings) {
    Tiff tiff(path);
    uint64_t next;
    const auto primary = tiff.directory(tiff.first, next);
    std::vector<uint64_t> candidates;
    if (next) candidates.push_back(next);
    std::set<uint64_t> visited{tiff.first};
    auto ifdOffset = [](double value) {
        if (value < 0 || value > std::numeric_limits<uint32_t>::max() || std::floor(value) != value)
            throw std::runtime_error("invalid TIFF directory offset");
        return uint64_t(value);
    };
    for (double offset : tiff.numbers(primary, 330)) candidates.push_back(ifdOffset(offset));
    auto raw = primary;
    auto score = [&](const Directory& d) {
        if (tiff.scalar(d, 262) != 32803 || tiff.scalar(d, 254) != 0) return 0.0;
        return tiff.scalar(d, 256) * tiff.scalar(d, 257);
    };
    double best = score(primary);
    for (size_t i = 0; i < candidates.size(); ++i) {
        if (i >= 1024 || candidates.size() > 1024) throw std::runtime_error("too many TIFF directories");
        const auto offset = candidates[i];
        if (!offset || !visited.insert(offset).second) continue;
        auto dir = tiff.directory(offset, next);
        if (next) candidates.push_back(next);
        for (double child : tiff.numbers(dir, 330)) candidates.push_back(ifdOffset(child));
        const double area = score(dir);
        if (area > best) {
            best = area;
            raw = std::move(dir);
        }
    }

    if (tiff.scalar(raw, 50711, 1) != 1)
        throw std::runtime_error("expected a rectangular DNG CFA layout");
    Metadata meta;
    auto imageTag = [&](unsigned tag) -> const Directory& { return raw.count(tag) ? raw : primary; };
    const auto active = tiff.numbers(imageTag(50829), 50829);
    if (!active.empty()) {
        if (active.size() != 4 || std::any_of(active.begin(), active.end(), [](double v) {
            return v < 0 || v > 65535 || v != std::floor(v);
        }) || active[0] >= active[2] || active[1] >= active[3])
            throw std::runtime_error("invalid DNG ActiveArea");
        meta.activeArea = Rect{int(active[0]), int(active[1]), int(active[2]), int(active[3])};
    }
    const uint64_t exifOffset = ifdOffset(tiff.scalar(primary, 34665, tiff.scalar(raw, 34665)));
    const auto exif = exifOffset ? tiff.directory(exifOffset, next) : Directory{};
    auto exposureTag = [&](unsigned tag) -> const Directory& { return exif.count(tag) ? exif : imageTag(tag); };
    meta.seconds = tiff.scalar(exposureTag(33434), 33434);
    meta.exposure = tiff.exposure(exposureTag(33434));
    meta.iso = tiff.scalar(exposureTag(34855), 34855);
    if (meta.iso == 65535 || meta.iso <= 0)
        meta.iso = tiff.scalar(exposureTag(34867), 34867, meta.iso);
    meta.aperture = tiff.scalar(exposureTag(33437), 33437);
    meta.captureDate = tiff.string(exposureTag(36867), 36867);
    if (meta.captureDate.empty()) meta.captureDate = tiff.string(imageTag(306), 306);
    meta.captureOffset = tiff.string(exposureTag(36881), 36881);
    if (meta.captureOffset.empty()) meta.captureOffset = tiff.string(exposureTag(36880), 36880);
    meta.description = tiff.string(imageTag(270), 270);
    meta.xmp = tiff.string(imageTag(700), 700, false);
    meta.baseline = tiff.scalar(imageTag(50730), 50730);
    meta.noise = tiff.numbers(imageTag(51041), 51041);
    auto planeColours = tiff.numbers(imageTag(50710), 50710);
    if (planeColours.empty()) planeColours = {0, 1, 2};
    {
        auto sorted = planeColours;
        std::sort(sorted.begin(), sorted.end());
        if (sorted != std::vector<double>{0, 1, 2})
            throw std::runtime_error("unsupported DNG plane colours");
    }
    for (size_t i = 0; i < 3; ++i) meta.planeColours[i] = int(planeColours[i]);
    const auto repeatCfa = tiff.numbers(imageTag(33421), 33421);
    const auto cfa = tiff.numbers(imageTag(33422), 33422);
    if (repeatCfa != std::vector<double>{2, 2} || cfa.size() != 4)
        throw std::runtime_error("expected a 2 × 2 DNG CFA pattern");
    // CFAPattern holds TIFF/EP colour codes (0 = red, 1 = green, 2 = blue), not
    // plane indices. CFAPlaneColor maps planes to those codes; it only reorders
    // per-plane tags such as AsShotNeutral, NoiseProfile and the colour matrices.
    meta.cfa.emplace();
    for (size_t i = 0; i < cfa.size(); ++i) {
        if (cfa[i] < 0 || cfa[i] >= 3 || cfa[i] != std::floor(cfa[i]))
            throw std::runtime_error("unsupported DNG CFA colour code");
        (*meta.cfa)[i] = int(cfa[i]);
    }
    if (meta.noise.size() == 6) {
        const auto noise = meta.noise;
        for (size_t i = 0; i < planeColours.size(); ++i)
            std::copy_n(noise.begin() + 2 * i, 2, meta.noise.begin() + 2 * int(planeColours[i]));
    }
    const auto neutral = tiff.numbers(imageTag(50728), 50728);
    if (neutral.size() == 3) {
        meta.neutral.emplace();
        for (size_t i = 0; i < 3; ++i) (*meta.neutral)[int(planeColours[i])] = neutral[i];
    }
    meta.black = tiff.numbers(imageTag(50714), 50714);
    meta.white = tiff.numbers(imageTag(50717), 50717);
    const auto repeat = tiff.numbers(imageTag(50713), 50713);
    if (!repeat.empty()) {
        if (repeat.size() != 2 || repeat[0] < 1 || repeat[1] < 1 ||
            repeat[0] > 65535 || repeat[1] > 65535 ||
            std::floor(repeat[0]) != repeat[0] || std::floor(repeat[1]) != repeat[1])
            throw std::runtime_error("invalid DNG BlackLevelRepeatDim");
        meta.blackRepeatRows = int(repeat[0]);
        meta.blackRepeatCols = int(repeat[1]);
    }
    meta.blackCols = tiff.numbers(imageTag(50715), 50715);
    meta.blackRows = tiff.numbers(imageTag(50716), 50716);
    constexpr unsigned opcodeTags[] = {51008, 51009, 51022};
    for (int i = 0; i < 3; ++i)
        parseOpcodes(meta, tiff.bytes(imageTag(opcodeTags[i]), opcodeTags[i]), i + 1, warnings);
    return meta;
}
}

#include "pipeline.h"

#include <Imath/half.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfChromaticitiesAttribute.h>
#include <OpenEXR/ImfFloatAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfMatrixAttribute.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfStandardAttributes.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfVecAttribute.h>
#include <OpenImageIO/imageio.h>

#include <stdexcept>

namespace rawfuse {
namespace exr = OPENEXR_IMF_NAMESPACE;
namespace imath = IMATH_NAMESPACE;
namespace {

void addChannel(exr::Header& header, exr::FrameBuffer& buffer, const char* name,
                exr::PixelType type, const void* data, size_t stride, int width) {
    header.channels().insert(name, exr::Channel(type));
    buffer.insert(name, exr::Slice(type, const_cast<char*>(static_cast<const char*>(data)),
                                  stride, stride * width));
}

void headerAttributes(exr::Header& header, const RawFrame& reference, const Colour& colour,
                      double scale, double referenceAperture, const Options& o) {
    header.insert("rawfuse:version", exr::StringAttribute(RAWFUSE_VERSION));
    header.insert("rawfuse:space", exr::StringAttribute(o.space == "aces" ? "aces-ap0" : "camera"));
    header.insert("rawfuse:exposureReference", exr::FloatAttribute(float(scale)));
    header.insert("rawfuse:referenceFNumber", exr::FloatAttribute(float(referenceAperture)));
    header.insert("rawfuse:whiteBalance", exr::V3fAttribute(
        imath::V3f(float(colour.wb[0]), float(colour.wb[1]), float(colour.wb[2]))));
    if (colour.matrix) {
        imath::M33f values;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) values[i][j] = float((*colour.matrix)[i][j]);
        header.insert("rawfuse:matrix", exr::M33fAttribute(values));
    }
    header.insert("rawfuse:report", exr::StringAttribute(o.report.filename().string()));
    if (o.space == "aces") {
        exr::addChromaticities(header, exr::Chromaticities(
            {0.7347f, 0.2653f}, {0, 1}, {0.0001f, -0.0770f}, {0.32168f, 0.33767f}));
    }
    const auto& meta = reference.metadata;
    if (!meta.captureDate.empty()) exr::addCapDate(header, meta.captureDate);
    if (o.normalize == "geomean") {
        exr::addExpTime(header, float(scale));
        exr::addIsoSpeed(header, 100);
        exr::addAperture(header, float(referenceAperture));
    }
}

} // namespace

namespace {

std::unique_ptr<OIIO::ImageInput> openImage(const fs::path& path, const Layout& layout) {
    OIIO::ImageSpec config;
    config.attribute("oiio:UnassociatedAlpha", 1);
    auto in = OIIO::ImageInput::open(path.string(), &config);
    if (!in) throw std::runtime_error("could not open " + path.string() + ": " + OIIO::geterror());
    const auto& spec = in->spec();
    if (spec.width != layout.width || spec.height != layout.height || spec.depth != 1 ||
        spec.deep || spec.nchannels < 1 || spec.nchannels > 4)
        throw std::runtime_error(path.string() + ": expected a full active-area image with 1 to 4 channels");
    return in;
}

} // namespace

int imageChannels(const fs::path& path, const Layout& layout) {
    return openImage(path, layout)->spec().nchannels;
}

std::vector<float> readImage(const fs::path& path, const Layout& layout, int& channels) {
    auto in = openImage(path, layout);
    channels = in->spec().nchannels;
    std::vector<float> result(layout.size() * channels);
    if (!in->read_image(0, 0, 0, channels, OIIO::TypeFloat, result.data()))
        throw std::runtime_error(path.string() + ": " + in->geterror());
    return result;
}

Mask readMask(const fs::path& path, const Layout& layout) {
    if (path.empty()) return {};
    int channels;
    const auto image = readImage(path, layout, channels);
    Mask mask(layout.size());
    for (size_t i = 0; i < layout.size(); ++i) {
        for (int c = 0; c < (channels >= 3 ? 3 : 1); ++c) {
            const float value = image[i * channels + c];
            if (!std::isfinite(value)) throw std::runtime_error(path.string() + ": non-finite mask value");
            mask[i] |= value != 0;
        }
    }
    return mask;
}

void writeImage(const fs::path& path, const Layout& layout, const Output& output,
                const RawFrame& ref, const Colour& colour,
                double scale, double referenceAperture, const Options& o) {
    exr::Header header(layout.width, layout.height);
    header.compression() = exr::ZIP_COMPRESSION;
    headerAttributes(header, ref, colour, scale, referenceAperture, o);
    exr::FrameBuffer buffer;
    constexpr const char* rgbNames[] = {"R", "G", "B"};
    constexpr const char* varNames[] = {"var.R", "var.G", "var.B"};
    std::vector<imath::half> frameCounts;
    for (int c = 0; c < 3; ++c) {
        addChannel(header, buffer, rgbNames[c], exr::FLOAT, output.rgb.data() + c, 3 * sizeof(float), layout.width);
        if (o.aux)
            addChannel(header, buffer, varNames[c], exr::FLOAT, output.variance.data() + c, 3 * sizeof(float), layout.width);
    }
    if (o.aux) {
        frameCounts.reserve(output.frames.size());
        for (const auto count : output.frames) frameCounts.emplace_back(float(std::min(int(count), 65504)));
        addChannel(header, buffer, "frames", exr::HALF, frameCounts.data(), sizeof(imath::half), layout.width);
        addChannel(header, buffer, "flags", exr::UINT, output.flags.data(), sizeof(uint32_t), layout.width);
        if (!output.sources.empty())
            addChannel(header, buffer, "source", exr::UINT, output.sources.data(), sizeof(uint32_t), layout.width);
    }
    exr::OutputFile file(path.string().c_str(), header, int(o.threads));
    file.setFrameBuffer(buffer);
    file.writePixels(layout.height);
}

void writeCfa(const fs::path& path, const Layout& layout, const std::vector<Site>& sites,
              const Options& o, const std::vector<uint32_t>& sources) {
    exr::Header header(layout.width, layout.height);
    header.compression() = exr::ZIP_COMPRESSION;
    header.insert("rawfuse:version", exr::StringAttribute(RAWFUSE_VERSION));
    header.insert("rawfuse:cfa", exr::StringAttribute(layout.pattern()));
    header.insert("rawfuse:space", exr::StringAttribute("camera-cfa-per-second"));
    std::vector<uint32_t> flags;
    flags.reserve(sites.size());
    for (const auto& site : sites) flags.push_back(site.flags);
    exr::FrameBuffer buffer;
    addChannel(header, buffer, "Y", exr::FLOAT, &sites.front().signal, sizeof(Site), layout.width);
    addChannel(header, buffer, "var", exr::FLOAT, &sites.front().variance, sizeof(Site), layout.width);
    addChannel(header, buffer, "flags", exr::UINT, flags.data(), sizeof(uint32_t), layout.width);
    if (!sources.empty())
        addChannel(header, buffer, "source", exr::UINT, sources.data(), sizeof(uint32_t), layout.width);
    exr::OutputFile file(path.string().c_str(), header, int(o.threads));
    file.setFrameBuffer(buffer);
    file.writePixels(layout.height);
}
}

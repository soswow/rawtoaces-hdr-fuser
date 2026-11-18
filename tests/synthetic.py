"""Generate small Bayer DNGs and check the public CLI with only Python's stdlib."""

import json
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import zlib

BLACK, WHITE = 4096, 64096
FULL = WHITE - BLACK
SIZE = 128
CFA = (0, 1, 1, 2)
S, O = 0.0002, 0.000001
XMP = b'<x:xmpmeta xmlns:x="adobe:ns:meta/">test</x:xmpmeta>\n'


def opcode(code, data):
    return struct.pack(">5I", 1, code, 0x01030000, 1, len(data)) + data


def dng(path, exposure, value, *, order="<", width=SIZE, height=SIZE,
        bad=b"", gains=b"", neutral=((1, 2), (1, 1), (5, 8)), subifd=False,
        active=None, final_opcode=b"", iso=100, aperture=(11, 5),
        cfa=CFA, black_pattern=None, black_rows=None, black_cols=None,
        plane_colours=(0, 1, 2), noise=(S, O), date=b"2026:01:02 03:04:05\0", offset=None):
    def pack(form, *values):
        return struct.pack(order + form, *values)

    def entry(tag, kind, form, *values):
        return tag, kind, len(values), pack(form * len(values), *values)

    def blob(tag, kind, value):
        return tag, kind, len(value), value

    def ratio(tag, *values, signed=False):
        return tag, 10 if signed else 5, len(values), b"".join(
            pack("ii" if signed else "II", a, b) for a, b in values
        )

    shared = [
        blob(271, 2, b"Rawfuse\0"), blob(272, 2, b"Test Bayer\0"),
        blob(270, 2, b"Synthetic bracket\0"),
        blob(50706, 1, bytes((1, 4, 0, 0))), blob(50707, 1, bytes((1, 1, 0, 0))),
        blob(50708, 2, b"Rawfuse Test Bayer\0"), blob(700, 1, XMP),
        blob(50710, 1, bytes(plane_colours)),
        ratio(50721, *((v, 1000) for c in plane_colours
                       for v in ((671, -237, -83), (-305, 1175, 142), (-51, 106, 507))[c]), signed=True),
        entry(50778, 3, "H", 21), ratio(50728, *(neutral[c] for c in plane_colours)),
        ratio(50730, (1, 4), signed=True),
    ]
    raw = [
        entry(254, 4, "I", 0), entry(256, 4, "I", width), entry(257, 4, "I", height),
        entry(258, 3, "H", 16), entry(259, 3, "H", 1), entry(262, 3, "H", 32803),
        entry(273, 4, "I", 8192), entry(274, 3, "H", 1), entry(277, 3, "H", 1),
        entry(278, 4, "I", height), entry(279, 4, "I", width * height * 2),
        entry(284, 3, "H", 1), entry(33421, 3, "H", 2, 2),
        # CFAPattern stores colour codes; CFAPlaneColor maps planes to the same codes.
        blob(33422, 1, bytes(cfa)),
        ratio(50714, (BLACK, 1)), entry(50717, 4, "I", WHITE),
        ratio(50718, (1, 1), (1, 1)), entry(50719, 4, "I", 0, 0),
        entry(50720, 4, "I", width, height),
    ]
    if noise:
        raw.append(entry(51041, 12, "d", *noise))
    if active:
        top, left, bottom, right = active
        raw = [tag for tag in raw if tag[0] != 50720]
        raw += [entry(50829, 4, "I", *active), entry(50720, 4, "I", right - left, bottom - top)]
    if black_pattern:
        raw = [tag for tag in raw if tag[0] != 50714]
        raw += [entry(50713, 3, "H", 2, 2), ratio(50714, *((v, 1) for v in black_pattern))]
    if black_rows:
        raw.append(ratio(50716, *((v, 1) for v in black_rows), signed=True))
    if black_cols:
        raw.append(ratio(50715, *((v, 1) for v in black_cols), signed=True))
    if final_opcode:
        raw.append(blob(51022, 7, final_opcode))
    if bad:
        raw.append(blob(51008, 7, bad))
    if gains:
        raw.append(blob(51009, 7, gains))
    exif = [
        ratio(33434, exposure), ratio(33437, aperture), entry(34855, 3, "H", iso),
    ]
    if date:
        exif.append(blob(36867, 2, date))
    if offset:
        exif.append(blob(36881, 2, offset))
    data = bytearray(8192)
    tail = 4096

    def directory(offset, entries):
        nonlocal tail
        table = pack("H", len(entries))
        for tag, kind, count, content in sorted(entries):
            table += pack("HHI", tag, kind, count)
            if len(content) > 4:
                table += pack("I", tail)
                data[tail:tail + len(content)] = content
                tail += len(content)
                tail += tail % 2
            else:
                table += content.ljust(4, b"\0")
        table += pack("I", 0)
        data[offset:offset + len(table)] = table
        assert tail <= 8192

    shared.append(entry(34665, 4, "I", 2048))
    data[:8] = (b"II" if order == "<" else b"MM") + pack("HI", 42, 8)
    if subifd:
        # IFD0 carries a tiny RGB thumbnail. CFA and raw-specific tags live in its SubIFD.
        preview = [
            entry(254, 4, "I", 1), entry(256, 4, "I", 1), entry(257, 4, "I", 1),
            entry(258, 3, "H", 8, 8, 8), entry(259, 3, "H", 1), entry(262, 3, "H", 2),
            entry(273, 4, "I", 8188), entry(277, 3, "H", 3), entry(278, 4, "I", 1),
            entry(279, 4, "I", 3), entry(330, 4, "I", 1024),
        ]
        directory(8, shared + preview)
        directory(1024, raw)
    else:
        directory(8, shared + raw)
    directory(2048, exif)
    for y in range(height):
        for x in range(width):
            ay, ax = y - (active[0] if active else 0), x - (active[1] if active else 0)
            black = black_pattern[(ay % 2) * 2 + ax % 2] if black_pattern else BLACK
            if black_rows and 0 <= ay < len(black_rows):
                black += black_rows[ay]
            if black_cols and 0 <= ax < len(black_cols):
                black += black_cols[ax]
            v = round(black + (WHITE - black) * value(y, x, cfa[(ay % 2) * 2 + ax % 2]))
            data += pack("H", max(0, min(WHITE, v)))
    path.write_bytes(data)


def read_exr(path):
    data = path.read_bytes()
    pos = 8

    def string():
        nonlocal pos
        end = data.index(0, pos)
        value = data[pos:end].decode()
        pos = end + 1
        return value

    attrs = {}
    while name := string():
        kind = string()
        size, = struct.unpack_from("<I", data, pos)
        pos += 4
        attrs[name] = (kind, data[pos:pos + size])
        pos += size
    xmin, ymin, xmax, ymax = struct.unpack("<4i", attrs["dataWindow"][1])
    assert (xmin, ymin) == (0, 0)
    width, height = xmax + 1, ymax + 1
    types = {}
    channel_data = attrs["channels"][1]
    offset = 0
    while channel_data[offset]:
        end = channel_data.index(0, offset)
        name = channel_data[offset:end].decode()
        types[name] = {0: "I", 1: "e", 2: "f"}[channel_data[end + 1]]
        offset = end + 17

    assert attrs["compression"][1] == b"\x03"  # ZIP, 16 scanlines per chunk
    channels = {name: [] for name in types}
    bytes_per_row = width * sum(struct.calcsize(fmt) for fmt in types.values())
    for first_row in range(0, height, 16):
        chunk, = struct.unpack_from("<Q", data, pos + (first_row // 16) * 8)
        y, length = struct.unpack_from("<2i", data, chunk)
        assert y == first_row
        raw = data[chunk + 8:chunk + 8 + length]
        rows = min(16, height - y)
        if length < rows * bytes_per_row:
            predicted = bytearray(zlib.decompress(raw))
            for i in range(1, len(predicted)):
                predicted[i] = (predicted[i - 1] + predicted[i] - 128) & 255
            split = (len(predicted) + 1) // 2
            flattened = bytearray(len(predicted))
            flattened[::2] = predicted[:split]
            flattened[1::2] = predicted[split:]
            raw = flattened
        offset = 0
        for _ in range(rows):
            for name, fmt in types.items():
                channels[name].extend(struct.unpack_from(f"<{width}{fmt}", raw, offset))
                offset += width * struct.calcsize(fmt)
    return attrs, channels


def close(actual, expected, tolerance=1e-4):
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{actual} differs from {expected} (tolerance {tolerance})")


def main(exe, folder):
    def fuse(name, *args, error=None):
        output = folder / f"{name}.exr"
        command = [exe, "--threads", "3", "--space", "camera", "--clip-dilate", "0",
                   "-o", str(output), *map(str, args)]
        result = subprocess.run(command, capture_output=True, text=True)
        if error:
            assert result.returncode == 2 and error in result.stderr, (command, result.stderr)
            return
        assert result.returncode == 0, (command, result.stderr)
        return json.loads(output.with_suffix(".json").read_text()), read_exr(output)

    dark_index = 10 * SIZE + 10
    point = opcode(5, struct.pack(">5I", 0, 1, 0, 10, 10))
    a, b = folder / "a.dng", folder / "b.dng"
    dng(a, (1, 10), lambda y, x, c: 0.7 if c == 0 else 0.4, bad=point)
    dng(b, (1, 20), lambda y, x, c: 0.35 if c == 0 else 0.2, order=">", subifd=True,
        neutral=((1, 4), (1, 1), (1, 2)))
    cfa = folder / "mosaic.exr"
    report, (attrs, rgb) = fuse("plain", b, a, "--normalize", "none", "--write-cfa", cfa, "--aux")
    assert attrs["rawfuse:space"][1] == b"camera" and "chromaticities" not in attrs
    assert not {"expTime", "isoSpeed", "aperture"} & attrs.keys()
    assert report["white_balance"]["method"] == "median"
    for w, expected in zip(report["white_balance"]["multipliers"], (3, 1, 1.8)):
        close(w, expected)
    assert sum("BaselineExposure" in w for w in report["warnings"]) == 1
    assert report["inputs"][0]["exposure_time"] == "1/20"
    assert report["inputs"][0]["xmp"].encode() == XMP
    assert report["inputs"][1]["hot_sites"]["missing"] == 1
    assert report["inputs"][1]["baseline_exposure"] == 0.25
    assert report["inputs"][1]["image_description"] == "Synthetic bracket"
    assert report["vignette"]["source"] == "absent"
    for c, expected in zip("RGB", (7, 4, 4)):
        close(rgb[c][1000], expected)
    close(rgb["frames"][1000], 2)
    _, mosaic = read_exr(cfa)
    full_weight = 0.1 ** 2 / (S * 0.7 + O) + 0.05 ** 2 / (S * 0.35 + O)
    close(mosaic["var"][0], 1 / full_weight, 4e-6)
    close(mosaic["Y"][dark_index], 7)

    different_settings = folder / "different-settings.dng"
    dng(different_settings, (1, 10), lambda y, x, c: 0.35 if c == 0 else 0.2,
        iso=200, aperture=(22, 5))
    report, (header, _) = fuse("exposure-tags", a, different_settings, "--reference", "2")
    def header_float(name):
        return struct.unpack("<f", header[name][1])[0]
    close(header_float("expTime"), math.sqrt(0.1 * 0.05))
    close(header_float("rawfuse:exposureReference"), report["exposure_reference"])
    close(header_float("isoSpeed"), 100)
    close(header_float("aperture"), 2.2)
    close(header_float("rawfuse:referenceFNumber"), 2.2)
    assert report["inputs"][1]["iso"] == 200 and report["reference_frame"] == 2

    shadow = folder / "shadow.dng"
    dng(shadow, (1, 10), lambda y, x, c: -0.002 if c == 0 else 0)
    report, (_, signed) = fuse("signed", shadow, "--aux")
    close(signed["R"][1000], -0.002)
    assert signed["G"][1000] == signed["B"][1000] == 0
    assert report["output_floor"] is None and report["merge"]["floored_fraction"] == 0
    assert not signed["flags"][1000] & 64
    report, (_, floor) = fuse("floor", shadow, "--floor", "0.0005", "--aux")
    for channel in "RGB":
        close(floor[channel][1000], 0.0005, 1e-7)
    assert report["merge"]["floored_fraction"] == 1
    assert floor["flags"][1000] & 64
    fuse("invalid-floor", a, "--floor", "0", error="--floor requires a positive float32 value")

    fuse("plain", a, error="output already exists")
    fuse("duplicate", a, a, error="duplicate input frame")
    fuse("alias", a, "--force", "--report", a, error="output aliases input")
    linked = folder / "hardlink.dng"
    os.link(a, linked)
    fuse("hardlinked", a, linked, error="duplicate input frame")

    _, (_, clipped) = fuse("clip", a, "--high", "0.65", "--highlights", "bound")
    close(clipped["R"][1000], 0.7)  # The sole near-white reading remains available.
    matrix = folder / "matrix.txt"
    matrix.write_text("0 1 0\n1 0 0\n0 0 1\n")
    report, (header, transformed) = fuse("colour", a, "--space", "aces", "--wb", "custom:2,1,3",
                                         "--matrix", f"custom:{matrix}", "--demosaic", "mhc")
    assert "chromaticities" in header and report["matrix"]["values"][0] == [0, 1, 0]
    for c, expected in zip("RGB", (0.4, 1.4, 1.2)):
        close(transformed[c][1000], expected)

    sat = folder / "saturated.dng"
    dng(sat, (1, 5), lambda y, x, c: 1 if x < 8 else 0.5)
    _, (_, highlight) = fuse("highlight", sat, "--space", "aces", "--wb", "custom:2,1,3",
                              "--matrix", f"custom:{matrix}", "--aux", "--highlights", "neutral")
    for c in "RGB":
        close(highlight[c][2 * SIZE + 2], 3 * 0.999)  # The clipping threshold, at the default margin
    assert highlight["flags"][2 * SIZE + 2] & 1
    # Reconstruction gives the fully clipped area the rim's colour (equal raw
    # values) at the lowest level its bounds allow.
    _, (_, highlight) = fuse("highlight-rim", sat, "--space", "aces", "--wb", "custom:2,1,3",
                              "--matrix", f"custom:{matrix}", "--aux")
    for c, expected in zip("RGB", (1, 2, 3)):
        close(highlight[c][2 * SIZE + 2], expected * 0.999)
    assert highlight["flags"][2 * SIZE + 2] & 1

    edge_short = folder / "edge-short.dng"
    edge_long = folder / "edge-long.dng"
    dng(edge_short, (1, 10), lambda y, x, c: 1 if x < 8 else 0.35)
    dng(edge_long, (1, 5), lambda y, x, c: 1 if x <= 8 else 0.7)
    edge_cfa = folder / "edge-cfa.exr"
    report, _ = fuse("edge", edge_short, edge_long, "--clip-dilate", "1",
                     "--write-cfa", edge_cfa, "--highlights", "bound")
    _, samples = read_exr(edge_cfa)
    edge = 10 * SIZE + 8
    close(samples["Y"][edge], 3.5)
    assert samples["flags"][edge] == 128
    assert samples["flags"][edge - 2] & 1
    assert report["merge"]["adjacent_clip_fallback_fraction"] > 0

    ramp = []
    rng = random.Random(55)
    for n in range(7):
        path = folder / f"ramp-{n}.dng"
        exposure = 4 ** n / 4096
        def noisy(y, x, c, exposure=exposure):
            signal = math.exp2(-9 + x / 127 * 21) * exposure
            return signal + rng.gauss(0, math.sqrt(S * signal + O)) if signal < 1 else signal
        dng(path, (4 ** n, 4096), noisy)
        ramp.append(path)
    _, (_, rampRGB) = fuse("ramp", *ramp, "--normalize", "none", "--no-deghost")
    errors = []
    for x in range(31, 111, 2):
        expected = math.exp2(-9 + x / 127 * 21)
        centre = sum(rampRGB["R"][y * SIZE + x - 1] for y in range(0, SIZE, 2)) / (SIZE / 2)
        target = expected / math.exp2(21 / 127)
        errors.append(centre / target - 1)
    assert abs(sum(errors) / len(errors)) < 0.005, errors

    # A fixed region differs in three frames. All Bayer phases come from frame 2.
    motion = []
    for n in range(3):
        path = folder / f"motion-{n}.dng"
        dng(path, (1, 10), lambda y, x, c, n=n: (0.6 if n == 1 else 0.2)
            if 38 <= x < 62 and 38 <= y < 62 else 0.2)
        motion.append(path)
    report, (_, moving) = fuse("motion", *motion, "--reference", "2", "--aux")
    assert report["reference_frame"] == 2
    close(moving["R"][50 * SIZE + 50], 0.6)
    assert moving["flags"][50 * SIZE + 50] & 4 and moving["frames"][50 * SIZE + 50] == 1
    close(moving["R"][25 * SIZE + 25], 0.2)

    sparse_model = folder / "sparse.csv"
    sparse_model.write_text("row,col,offset_dn,rate_dn_per_s\n10,10,300,6000\n")
    hot = folder / "hot.dng"
    dng(hot, (1, 10), lambda y, x, c: 0.4 + (900 / FULL if (y, x) == (10, 10) else 0),
        bad=point)
    report, (_, repaired) = fuse("repair", hot, "--aux", "--dark-model", sparse_model)
    close(repaired["R"][dark_index], 0.4)
    assert report["groups"][0]["dark_model"]["method"] == "unavailable"
    assert repaired["flags"][dark_index] & (8 | 16) == (8 | 16)

    offsets = {(y, x): 2400.0 for y in range(10, 70, 8) for x in (10, 18, 26, 34, 42)}
    rates = {(y, x): 160000.0 for y in range(10, 70, 8) for x in (58, 66, 74, 82, 90)}
    dead = (10, 10)
    dark_model = folder / "dark.csv"
    defects = offsets.keys() | rates.keys()
    dark_model.write_text("row,col,offset_dn,rate_dn_per_s\n" + "".join(
        f"{y},{x},{offsets.get((y, x), 0)},{rates.get((y, x), 0)}\n" for y, x in sorted(defects)))
    hot_frames = []
    for n, t in ((20, 0.05), (10, 0.1)):
        hot_frame = folder / f"hot-{n}.dng"
        def hot_value(y, x, c):
            location = (y, x)
            gain_error = 1 + 0.14 * (1 if (y // 8) % 2 else -1)
            dark = (0.65 * offsets.get(location, 0) + 0.25 * rates.get(location, 0) * t) * gain_error
            return 0.3 * t + (dark / FULL if location != dead else 0)
        dng(hot_frame, (1, n), hot_value)
        hot_frames.append(hot_frame)
    report, (_, corrected) = fuse("dark", *hot_frames, "--dark-model", dark_model, "--aux",
                                  "--normalize", "none", "--noise", "0.000001,0.000000001")
    fit = report["groups"][0]["dark_model"]
    assert fit["method"] == "separate_scales"
    close(fit["offset_scale"], 0.65, 0.017)
    close(fit["rate_scale"], 0.25, 0.007)
    assert fit["relative_scatter"] > 0.1
    for frame in report["inputs"]:
        assert frame["hot_sites"]["kept"] == len(defects) - 1
        assert frame["hot_sites"]["model_rejected"] == 1
    assert corrected["flags"][dark_index] & (8 | 16) == (8 | 16)
    close(sum(corrected["R"][y * SIZE + x] for y, x in defects) / len(defects), 0.3, 0.01)

    crop = folder / "active.dng"
    crop_bad = opcode(5, struct.pack(">9I", 0, 1, 1, 62, 64, 6, 8, 26, 28))
    crop_dark = folder / "active-dark.csv"
    crop_dark.write_text("row,col,offset_dn,rate_dn_per_s\n60,60,600,0\n")
    dng(crop, (1, 10), lambda y, x, c: 0.25 if 2 <= y < 126 and 4 <= x < 124 else 0,
        active=(2, 4, 126, 124), bad=crop_bad, final_opcode=opcode(222, b""))
    report, (_, active) = fuse("active", crop, "--aux", "--dark-model", crop_dark)
    assert (report["width"], report["height"]) == (120, 124)
    assert any("DNG opcode 222 in list 3" in warning for warning in report["warnings"])
    assert active["flags"][10 * 120 + 10] & (8 | 16) == (8 | 16)
    assert active["flags"][60 * 120 + 60] & (8 | 16) == (8 | 16)
    assert report["inputs"][0]["hot_sites"]["missing"] == 401
    for name in "RGB":
        close(active[name][10 * 120 + 10], 0.25)
    fuse("incompatible", a, crop, error="camera, CFA or active area differ")

    gain_data = struct.pack(">10I4dI4f", 0, 0, SIZE, SIZE, 0, 1, 1, 1, 2, 2,
                            1 - 1 / SIZE, 1 - 1 / SIZE, 0.5 / SIZE, 0.5 / SIZE, 1, 1, 2, 3, 4)
    flat = folder / "flat.dng"
    def flat_value(y, x, c):
        return 0.4 / (1 + 2 * y / (SIZE - 1) + x / (SIZE - 1))
    dng(flat, (1, 10), flat_value, gains=opcode(9, gain_data), order=">")
    report, (_, flattened) = fuse("gain", flat, "--aux")
    assert report["vignette"]["gain_maps"] == 1
    for name in "RGB":
        for i in (0, SIZE * SIZE - 1, 5048):
            close(flattened[name][i], 0.4)
    close(report["vignette"]["max_gain"]["B"], 4)

    base, filter_frame, no_overlap = (folder / f"nd-{s}.dng" for s in ("base", "filter", "no-overlap"))
    dims = {"width": 256, "height": 256}
    evs = (9.7, 9.8, 9.9)
    mask = folder / "mask.pgm"
    def plain_signal(y, x):
        bright = 192 <= x < 224 or (16 <= y < 24 and 16 <= x < 24)
        return (2 if bright else 0.65) if x < 224 else 0.25
    dng(base, (1, 10000), lambda y, x, c: plain_signal(y, x), **dims)
    dng(filter_frame, (1, 40), lambda y, x, c: 250 * plain_signal(y, x) /
        math.exp2(evs[c] + (0.55 if y < 64 and x < 64 else 0)) + 0.00015 if x < 224 else 1, **dims)
    uncovered = folder / "nd-uncovered.dng"
    drift = folder / "nd-drift.dng"
    dng(uncovered, (1, 10000), lambda y, x, c: plain_signal(y, x) /
        math.exp2(evs[c]) if x < 224 else 0.25, **dims)
    dng(drift, (1, 10000), lambda y, x, c: plain_signal(y, x) /
        math.exp2(evs[c]) if x < 224 else 0.125, **dims)
    dng(no_overlap, (1, 10000), lambda y, x, c: 0.005, **dims)
    mask.write_bytes(b"P5\n256 256\n255\n" + (bytes([255]) * 224 + bytes(32)) * 256)
    unfiltered = folder / "unfiltered.pgm"
    unfiltered.write_bytes(b"P5\n256 256\n255\n" + (bytes(224) + bytes([255]) * 32) * 256)
    nd_manifest = folder / "nd-groups.json"

    def nd_groups(frames, ev="auto", drift_mask=True, match=None):
        filtered = {"name": "filter", "frames": [f.name for f in frames], "ev": ev, "valid_mask": mask.name}
        if drift_mask:
            filtered["drift_mask"] = unfiltered.name
        if match:
            filtered["match"] = match
        nd_manifest.write_text(json.dumps({"groups": [{"name": "plain", "frames": [base.name]}, filtered]}))
        return nd_manifest

    report, (_, nd_image) = fuse("nd", "--groups", nd_groups([filter_frame, uncovered]), "--aux")
    calibration = report["groups"][1]["calibration"]
    for c, ev in zip("RGB", evs):
        close(report["groups"][1]["ev"][c], ev, 0.001)
        close(nd_image[c][130 * 256 + 208], 2, 0.005)
        close(nd_image[c][130 * 256 + 242], 0.25)
    # Without a contribution mask the filtered group only fills clipped sites.
    assert nd_image["source"][130 * 256 + 208] == 2 and nd_image["source"][130 * 256 + 130] == 1
    assert nd_image["flags"][130 * 256 + 208] & 2
    # A deviating calibration block is reported, but still fills its clipped sites.
    assert {"row": 0, "col": 0} in calibration["blocks_flagged"]
    assert any("calibration blocks differ" in w for w in report["warnings"])
    # Offset matching against the surrounding plain samples pulls that fill towards
    # the truth (2); without it the deviating block's own value is used.
    deviating = 2 / math.exp2(0.55)
    assert deviating + 0.1 < nd_image["G"][20 * 256 + 20] < 2
    assert report["groups"][1]["offset_match"]["corrected_sites"] > 0
    unmatched, (_, raw_fill) = fuse("nd-unmatched", "--groups", nd_groups([filter_frame, uncovered], match="none"),
                                    "--aux")
    close(raw_fill["G"][20 * 256 + 20], deviating, 0.005)
    assert unmatched["groups"][1]["offset_match"] is None
    close(calibration["scene_drift_ev"], 0)
    assert calibration["scene_drift_samples"] >= 100
    fuse("nd-fail", "--groups", nd_groups([no_overlap], drift_mask=False),
         error="filter: insufficient overlap with the calibration reference (R/G/B 0/0/0")
    fuse("nd-drift-fail", "--groups", nd_groups([filter_frame, drift], 9.8),
         error="scene drift of -1.00 EV exceeds 0.20 EV")
    assert not (folder / "nd-drift-fail.exr").exists()
    auto, _ = fuse("nd-auto-drift", "--groups", nd_groups([filter_frame, drift]))
    assert any("scene drift" in w for w in auto["warnings"])
    report, _ = fuse("nd-manual", "--groups", nd_groups([filter_frame, drift], 9.8), "--allow-scene-drift")
    close(report["groups"][1]["calibration"]["scene_drift_ev"], -1)
    assert report["groups"][1]["calibration"]["method"] == "manual"
    assert any("scene drift" in w for w in report["warnings"])
    fuse("drift-without-manual-ev", "--groups", nd_groups([filter_frame, drift]), "--allow-scene-drift",
         error="--allow-scene-drift requires a group with a manual ev and a drift_mask")


    def write_mask(name, width, selected):
        path = folder / (name + ".pgm")
        row = bytes(255 if selected(x) else 0 for x in range(width))
        path.write_bytes(f"P5\n{width} {width}\n255\n".encode() + row * width)
        return path.name

    # Independent masks for a blocked bracket, an open bracket and a filtered view.
    base_frame, open_frame, nd_frame = (folder / f"group-{g}.dng" for g in ("base", "open", "nd"))
    dng(base_frame, (1, 10), lambda y, x, c: 0.2 if x < 48 or 85 <= x < 100 else 1)
    dng(open_frame, (1, 10), lambda y, x, c: 1 if 65 <= x < 85 else 0.3)
    group_evs = (2, 3, 4)
    dng(nd_frame, (1, 10), lambda y, x, c: 0.6 / math.exp2(group_evs[c]) if 69 <= x < 78 else 0.1)
    blocked_valid = write_mask("base-valid", SIZE, lambda x: x < 48 or 85 <= x < 100)
    open_contribution = write_mask("open-contribution", SIZE, lambda x: 42 <= x < 96)
    open_valid = write_mask("open-valid", SIZE, lambda x: x < 100)
    nd_valid = write_mask("nd-valid", SIZE, lambda x: 69 <= x < 78)
    nd_patch = write_mask("nd-patch", SIZE, lambda x: 56 <= x < 78)
    manifest = folder / "groups.json"
    # The brackets read deliberately different levels to show which one wins, so
    # offset matching (tested separately) is off.
    entries = [
        {"name": "base", "frames": [base_frame.name], "valid_mask": blocked_valid},
        {"name": "clear", "frames": [open_frame.name], "valid_mask": open_valid,
         "contribution_mask": open_contribution, "match": "none"},
        {"name": "filter", "frames": [nd_frame.name], "valid_mask": nd_valid,
         "contribution_mask": nd_patch, "ev": list(group_evs), "match": "none"},
    ]
    manifest.write_text(json.dumps({"groups": entries}))
    mosaic = folder / "groups-cfa.exr"
    report, (_, composed) = fuse("group-manual", "--groups", manifest, "--clip-dilate", "1",
                                 "--aux", "--highlights", "bound", "--write-cfa", mosaic)
    assert [g["name"] for g in report["groups"]] == ["base", "clear", "filter"]
    assert report["groups"][2]["ev"] == dict(zip("RGB", group_evs))
    assert report["groups"][2]["measured_sites"] == SIZE * 9
    assert report["groups"][2]["calibration"] is None
    assert report["merge"]["uncovered_fraction"] == 28 / SIZE
    row = 40 * SIZE
    for c in "RGB":
        close(composed[c][row + 20], 0.2)
        close(composed[c][row + 46], 0.3)  # Later clean measurement wins in overlap.
        close(composed[c][row + 60], 0.3)  # ND contribution outside validity is ignored.
        close(composed[c][row + 73], 0.6, 0.001)  # ND beats the clipped clear bracket.
        close(composed[c][row + 82], 0.999)  # Clipped lower bound at the clipping threshold
        close(composed[c][row + 99], 0.2)  # Boundary interpolation skips uncovered taps.
        close(composed[c][row + 108], 0)
    assert composed["source"][row + 20] == 1
    assert composed["source"][row + 42] == 65535
    assert composed["source"][row + 73] == 3
    assert composed["source"][row + 108] == 0
    assert composed["flags"][row + 82] & 1
    assert composed["flags"][row + 73] & 32  # Replaces a higher clipping bound, flagged.
    assert composed["flags"][row + 108] == 16 | 256
    _, cfa_data = read_exr(mosaic)
    assert [cfa_data["source"][row + x] for x in (20, 42, 60, 73, 82, 108)] == [1, 2, 2, 3, 2, 0]
    fuse("group-conflict", "--groups", manifest, a, error="takes input frames")
    fuse("group-alias", "--groups", manifest, "--report", folder / nd_patch,
         error="output aliases input")
    entries[2]["valid_mask"] = "missing-mask.png"
    manifest.write_text(json.dumps({"groups": entries}))
    fuse("group-missing-mask", "--groups", manifest, error="missing-mask.png")

    # Calibration uses eligible sites outside both contribution masks and excludes an unfiltered region.
    fake_base = folder / "unrelated-base.dng"
    dng(fake_base, (1, 10000), lambda y, x, c: 0.33, **dims)
    calibration = write_mask("nd-cal", 256, lambda x: 64 <= x < 192 or 230 <= x)
    filtered = write_mask("nd-filter-valid", 256, lambda x: x < 224)
    replacement = write_mask("nd-only-highlight", 256, lambda x: 196 <= x < 224)
    clear_patch = write_mask("clear-only-highlight", 256, lambda x: 190 <= x < 224)
    # The first group is unrelated to the scene, so offset matching against it is off.
    entries = [
        {"name": "base", "frames": [fake_base.name]},
        {"name": "filter", "frames": [filter_frame.name], "contribution_mask": replacement,
         "ev": "auto", "calibration_reference": "clear", "calibration_mask": calibration,
         "valid_mask": filtered, "match": "none"},
        {"name": "clear", "frames": [base.name], "contribution_mask": clear_patch, "match": "none"},
    ]
    manifest.write_text(json.dumps({"groups": entries}))
    report, (_, composed) = fuse("group-auto", "--groups", manifest, "--aux")
    for c, ev in zip("RGB", evs):
        close(report["groups"][1]["ev"][c], ev, 0.001)
        close(composed[c][130 * 256 + 208], 2, 0.005)
        close(composed[c][130 * 256 + 130], 0.33)
    assert composed["source"][130 * 256 + 208] == 2
    assert report["groups"][1]["calibration"]["scene_drift_samples"] == 0
    entries[1]["calibration_mask"] = replacement
    manifest.write_text(json.dumps({"groups": entries}))
    fuse("group-auto-empty", "--groups", manifest, error="filter: insufficient overlap")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="rawfuse-tests-") as directory:
        main(str(Path(sys.argv[1]).resolve()), Path(directory))

"""Regressions found during the full pipeline review; no external Python packages."""

import calendar
import json
import itertools
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile
import unittest

from synthetic import SIZE, close, dng, opcode, read_exr


EXE = str(Path(sys.argv.pop(1)).resolve())


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix="rawfuse-review-")
        self.addCleanup(self.scratch.cleanup)
        self.folder = Path(self.scratch.name)

    def fuse(self, name, *args, error=None):
        output = self.folder / (name + ".exr")
        result = subprocess.run([EXE, "--threads", "3", "--space", "camera",
                                 "--clip-dilate", "0", "-o", str(output), *map(str, args)],
                                capture_output=True, text=True)
        if error:
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn(error, result.stderr)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix(".json").exists())
            return
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(output.with_suffix(".json").read_text()), read_exr(output)[1]

    def test_odd_active_area_preserves_grid_and_black_levels(self):
        # DNG CFA and black patterns start at ActiveArea, including odd margins.
        for cfa in ((0, 1, 1, 2), (2, 1, 1, 0), (1, 0, 2, 1), (1, 2, 0, 1)):
            with self.subTest(cfa=cfa):
                frame = self.folder / "crop.dng"
                dng(frame, (1, 10), lambda y, x, c: (0.2, 0.4, 0.6)[c],
                    active=(3, 5, 125, 123), cfa=cfa, black_pattern=(4000, 4100, 4200, 4300),
                    black_rows=[y % 3 for y in range(122)], black_cols=[x % 5 for x in range(118)])
                report, rgb = self.fuse("crop", frame, "--force")
                self.assertEqual((report["width"], report["height"]), (118, 122))
                self.assertEqual(report["cfa"], "".join("RGB"[c] for c in cfa))
                self.assertFalse(any("disagrees with LibRaw" in w for w in report["warnings"]))
                for c, expected in zip("RGB", (0.2, 0.4, 0.6)):
                    for i in (0, 1000, 118 * 122 - 1):
                        close(rgb[c][i], expected)

    def test_repaired_group_provenance(self):
        frame = self.folder / "filtered.dng"
        dng(frame, (1, 10), lambda y, x, c: 0.2,
            bad=opcode(5, struct.pack(">5I", 0, 1, 0, 10, 10)))
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [{"name": "filter", "frames": [frame.name], "ev": 1}]}))
        matrix = self.folder / "red.txt"
        matrix.write_text("1 0 0\n1 0 0\n1 0 0\n")
        _, rgb = self.fuse("repaired", "--groups", manifest, "--aux", "--space", "aces",
                           "--wb", "custom:1,1,1", "--matrix", "custom:" + str(matrix))
        i = 10 * SIZE + 10
        close(rgb["R"][i], 0.4)
        self.assertEqual(rgb["source"][i], 1)
        self.assertEqual(rgb["flags"][i] & (2 | 8 | 16), 2 | 8 | 16)

        # Donors on opposite sides of a group boundary produce mixed provenance.
        other = self.folder / "other.dng"
        dng(other, (1, 10), lambda y, x, c: 0.4,
            bad=opcode(5, struct.pack(">5I", 0, 1, 0, 10, 10)))
        mask = self.folder / "right.pgm"
        mask.write_bytes(f"P5\n{SIZE} {SIZE}\n255\n".encode() +
                         (bytes(11) + bytes([255]) * (SIZE - 11)) * SIZE)
        groups = json.loads(manifest.read_text())
        groups["groups"].append({"name": "other", "frames": [other.name], "contribution_mask": mask.name})
        manifest.write_text(json.dumps(groups))
        _, rgb = self.fuse("mixed", "--groups", manifest, "--aux", "--space", "aces",
                           "--wb", "custom:1,1,1", "--matrix", "custom:" + str(matrix))
        close(rgb["R"][i], 0.4)
        self.assertEqual(rgb["source"][i], 65535)

    def test_cropped_gain_map_and_defect_coordinates(self):
        frame = self.folder / "crop.dng"
        height, width = 122, 118
        gain = struct.pack(">10I4dI4f", 0, 0, height, width, 0, 1, 1, 1, 2, 2,
                           1 - 1 / height, 1 - 1 / width, 0.5 / height, 0.5 / width, 1, 1, 2, 3, 4)
        dng(frame, (1, 10), lambda y, x, c: 0.4 / (1 + 2 * (y - 3) / (height - 1) +
                                                                 (x - 5) / (width - 1)),
            active=(3, 5, 125, 123), gains=opcode(9, gain),
            bad=opcode(5, struct.pack(">5I", 0, 1, 0, 13, 15)))
        _, rgb = self.fuse("gain", frame, "--aux")
        self.assertEqual(rgb["flags"][10 * width + 10] & (8 | 16), 8 | 16)
        for c in "RGB":
            for i in (0, 10 * width + 10, height * width - 1):
                close(rgb[c][i], 0.4)

    def test_dng_plane_order(self):
        frame = self.folder / "bgr.dng"
        dng(frame, (1, 10), lambda y, x, c: (0.2, 0.4, 0.6)[c], plane_colours=(2, 1, 0),
            noise=(0.003, 0.00003, 0.002, 0.00002, 0.001, 0.00001))
        # Explicit frame WB exercises the metadata-to-RGB mapping independently
        # of rawtoaces' own metadata WB reader.
        report, rgb = self.fuse("bgr", frame, "--wb", "frame:1")
        # CFAPattern holds colour codes, independent of CFAPlaneColor.
        self.assertEqual(report["cfa"], "RGGB")
        for c, expected in zip("RGB", (0.2, 0.4, 0.6)):
            close(rgb[c][1000], expected)
        for actual, expected in zip(report["white_balance"]["multipliers"], (2, 1, 1.6)):
            close(actual, expected)
        for c, expected in zip("RGB", (0.001, 0.002, 0.003)):
            close(report["noise_profile"]["planes"][c]["S"], expected)

        # Encoding the same physical camera with different plane indices must
        # leave AP0, custom transforms and reported RGB WB unchanged.
        matrix = self.folder / "identity.txt"
        matrix.write_text("1 0 0\n0 1 0\n0 0 1\n")
        expected = {}
        for planes in itertools.permutations((0, 1, 2)):
            dng(frame, (1, 10), lambda y, x, c: (0.2, 0.4, 0.6)[c], plane_colours=planes,
                cfa=(1, 2, 0, 1))
            if planes != (0, 1, 2):
                self.fuse("unsupported-adobe", frame, "--space", "aces", "--matrix", "adobe",
                          error="Adobe matrix mode requires RGB DNG plane order")
            for method in ("metadata", "auto", "custom:" + str(matrix)):
                for wb in ("reference", "frame:1", "custom:3,1,2"):
                    with self.subTest(planes=planes, matrix=method, wb=wb):
                        report, rgb = self.fuse("plane-aces", frame, "--space", "aces", "--force",
                                               "--wb", wb, "--matrix", method)
                        self.assertEqual(report["cfa"], "GBRG")
                        values = report["white_balance"]["multipliers"] + [rgb[c][1000] for c in "RGB"]
                        key = method, wb
                        if key not in expected:
                            expected[key] = values
                        for actual, reference in zip(values, expected[key]):
                            close(actual, reference, 2e-5)

    def test_neutral_highlight_provenance(self):
        base, nd = self.folder / "base.dng", self.folder / "nd.dng"
        dng(base, (1, 10), lambda y, x, c: 1 if c == 0 else 0.2)
        dng(nd, (1, 10), lambda y, x, c: 0.4)
        mask = self.folder / "nonred.pgm"
        mask.write_bytes(f"P5\n{SIZE} {SIZE}\n255\n".encode() +
                         bytes(0 if y % 2 == x % 2 == 0 else 255
                               for y in range(SIZE) for x in range(SIZE)))
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [
            {"name": "base", "frames": [base.name]},
            {"name": "nd", "frames": [nd.name], "ev": 1, "valid_mask": mask.name,
             "contribution_mask": mask.name}]}))
        matrix = self.folder / "red.txt"
        matrix.write_text("1 0 0\n1 0 0\n1 0 0\n")
        _, rgb = self.fuse("highlight", "--groups", manifest, "--aux", "--space", "aces",
                           "--wb", "custom:1,2,1", "--matrix", "custom:" + str(matrix))
        i = 10 * SIZE + 10
        close(rgb["R"][i], 1.6)
        self.assertEqual(rgb["source"][i], 65535)
        self.assertEqual(rgb["flags"][i] & (1 | 2), 1 | 2)

    def test_reconstructed_highlights_take_the_rim_colour(self):
        # A disc clips R and G but not B. Its colour comes from the unclipped rim
        # and its level from B, so the bounds do not decide the colour.
        frame = self.folder / "disc.dng"
        colour = (0.6, 0.4, 0.1)
        # A lens-blurred profile: quads straddling a hard step would mix colours.
        dng(frame, (1, 10), lambda y, x, c: colour[c] * (1 + 3 / (1 + (math.hypot(y - 64, x - 64) / 10) ** 2)))
        report, rgb = self.fuse("disc", frame, "--aux", "--wb", "custom:1,1,1")
        i = 64 * SIZE + 64
        close(rgb["R"][i] / rgb["B"][i], 6, 0.03)
        close(rgb["G"][i] / rgb["B"][i], 4, 0.03)
        self.assertEqual(rgb["flags"][i] & (1 | 512), 1 | 512)
        self.assertGreater(report["merge"]["reconstructed_fraction"], 0)

        # The old neutral mode raises G to the R bound and stays available.
        _, rgb = self.fuse("disc-neutral", frame, "--aux", "--wb", "custom:1,1,1", "--highlights", "neutral")
        close(rgb["R"][i], rgb["G"][i])
        self.assertEqual(rgb["flags"][i] & 512, 0)

    def test_steep_edges_keep_their_hue(self):
        # A small bright disc whose brightness falls about tenfold per site at
        # the rim, like the sun. A straight average of the neighbours is
        # dominated by the bright side and gives each site the wrong hue.
        frame = self.folder / "sun.dng"
        colour = (1.0, 0.6, 0.25)
        def level(r):
            return 0.008 * 100 ** (1 / (1 + math.exp((r - 12) / 0.7)))
        # Low noise like a real sensor: log-space estimates need 8 σ everywhere in the window.
        dng(frame, (1, 10), lambda y, x, c: colour[c] * level(math.hypot(y - 64.3, x - 63.6)),
            noise=(1e-5, 6e-9))
        true = (math.log(colour[0] / colour[1]), math.log(colour[2] / colour[1]))

        def hue_errors(rgb, pixels):
            return [max(abs(math.log(rgb["R"][i] / rgb["G"][i]) - true[0]),
                        abs(math.log(rgb["B"][i] / rgb["G"][i]) - true[1])) for i in pixels]

        report, rgb = self.fuse("sun", frame, "--aux", "--wb", "custom:1,1,1")
        steep = [i for i in range(SIZE * SIZE) if rgb["flags"][i] & 1024]
        self.assertGreater(len(steep), 40)
        self.assertEqual(report["edges"]["pixels"], len(steep))
        ratio = hue_errors(rgb, steep)

        linear_report, linear = self.fuse("sun-linear", frame, "--aux", "--wb", "custom:1,1,1", "--edges", "linear")
        self.assertEqual(linear_report["edges"]["pixels"], 0)
        self.assertFalse(any(linear["flags"][i] & 1024 for i in range(SIZE * SIZE)))
        straight = hue_errors(linear, steep)
        # Hue error in log units: about 0.02 median and 0.25 at worst (inside the
        # blend ramp) against about 0.4 and 0.95 for the straight average.
        median = lambda values: sorted(values)[len(values) // 2]
        self.assertLess(median(ratio), 0.05)
        self.assertLess(max(ratio), 0.35)
        self.assertGreater(median(straight), 0.3)
        self.assertGreater(max(straight), 0.8)
        # Away from the rim both modes are the same demosaic.
        for i in (64 * SIZE + 64, 10 * SIZE + 10):
            self.assertEqual(rgb["flags"][i] & 1024, 0)
            for c in "RGB":
                self.assertEqual(rgb[c][i], linear[c][i])

    def seam(self, frame, name, selected, *args):
        # A second group with identical data that replaces only `selected` sites:
        # values stay the same, provenance must follow what each estimate read.
        copy = self.folder / ("copy-" + frame.name)
        copy.write_bytes(frame.read_bytes())
        mask = self.folder / (name + ".pgm")
        mask.write_bytes(f"P5\n{SIZE} {SIZE}\n255\n".encode() +
                         bytes(255 if selected(y, x) else 0 for y in range(SIZE) for x in range(SIZE)))
        manifest = self.folder / (name + "-groups.json")
        manifest.write_text(json.dumps({"groups": [
            {"name": "first", "frames": [frame.name]},
            {"name": "second", "frames": [copy.name], "contribution_mask": mask.name, "match": "none"}]}))
        return self.fuse(name, "--groups", manifest, "--aux", *args)

    def test_edge_estimate_sources_include_every_site_it_read(self):
        frame = self.folder / "sun.dng"
        colour = (1.0, 0.6, 0.25)
        level = lambda r: 0.008 * 100 ** (1 / (1 + math.exp((r - 12) / 0.7)))
        dng(frame, (1, 10), lambda y, x, c: colour[c] * level(math.hypot(y - 64.3, x - 63.6)),
            noise=(1e-5, 6e-9))
        column = 52  # The disc's left rim
        _, rgb = self.seam(frame, "edge-seam", lambda y, x: x == column, "--wb", "custom:1,1,1")
        # Bilinear taps reach one column; the log-space estimate reads two or three.
        beside = [y * SIZE + column + 2 for y in range(SIZE) if rgb["flags"][y * SIZE + column + 2] & 1024]
        self.assertGreater(len(beside), 3)
        for i in beside:
            self.assertEqual(rgb["source"][i], 65535)

    def test_reconstructed_highlight_sources_include_the_rim(self):
        frame = self.folder / "disc.dng"
        colour = (0.6, 0.4, 0.1)
        dng(frame, (1, 10), lambda y, x, c: colour[c] * (1 + 3 / (1 + (math.hypot(y - 64, x - 64) / 10) ** 2)))
        # Only the unclipped rim, which sets the disc's colour, comes from the second group.
        _, rgb = self.seam(frame, "rim-seam", lambda y, x: 20 <= math.hypot(y - 64, x - 64) <= 30,
                           "--wb", "custom:1,1,1")
        i = 64 * SIZE + 64
        self.assertEqual(rgb["flags"][i] & (1 | 512), 1 | 512)
        self.assertEqual(rgb["source"][i], 65535)

    def test_variance_overflow_is_rejected_before_publication(self):
        frame = self.folder / "frame.dng"
        dng(frame, (1, 10), lambda y, x, c: 0.2)
        matrix = self.folder / "large.txt"
        matrix.write_text("1e30 0 0\n0 1e30 0\n0 0 1e30\n")
        self.fuse("overflow", frame, "--aux", "--space", "aces", "--matrix", "custom:" + str(matrix),
                  error="output variance overflows float32")

        gain = struct.pack(">10I4dIf", 0, 0, SIZE, SIZE, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1e30)
        dng(frame, (1, 10), lambda y, x, c: 0.2, gains=opcode(9, gain))
        cfa = self.folder / "overflow-cfa.exr"
        self.fuse("gain-overflow", frame, "--write-cfa", cfa,
                  error="merged CFA signal or variance overflows float32")
        self.assertFalse(cfa.exists())

    def test_later_group_keeps_stronger_lower_bound(self):
        short, long = self.folder / "short.dng", self.folder / "long.dng"
        dng(short, (1, 1000), lambda y, x, c: 1)
        dng(long, (1, 10), lambda y, x, c: 1)
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [
            {"name": "short", "frames": [short.name]}, {"name": "long", "frames": [long.name]}]}))
        report, rgb = self.fuse("bounds", "--groups", manifest, "--aux", "--highlights", "bound")
        # Clipping at 1/1000 s bounds the rate at 999 per second (the default
        # 0.001 clip margin); the later 1/10 s bound is weaker and must not replace it.
        i = 10 * SIZE + 10
        close(rgb["R"][i], 0.999)
        self.assertEqual(rgb["source"][i], 1)
        self.assertEqual(rgb["flags"][i] & 1, 1)
        self.assertEqual(report["groups"][0]["lower_bound_sites"], SIZE * SIZE)

    def test_clipped_bound_is_the_clipping_threshold(self):
        # A reading at 95% of range clipped by a 10% margin is known only to be at
        # least 90%, and a reading at the white level likewise.
        for level in (0.95, 1):
            with self.subTest(level=level):
                frame = self.folder / "margin.dng"
                dng(frame, (1, 1), lambda y, x, c: level)
                _, rgb = self.fuse("margin", frame, "--force", "--aux", "--highlights", "bound",
                                   "--clip-margin", "0.10")
                i = 10 * SIZE + 10
                close(rgb["G"][i], 0.9)
                self.assertEqual(rgb["flags"][i] & 1, 1)

    def test_sidecars_fail_before_remaining_frames_with_warnings(self):
        good, broken = self.folder / "good.dng", self.folder / "broken.dng"
        dng(good, (1, 10), lambda y, x, c: 0.2, final_opcode=opcode(222, b""))
        broken.write_bytes(b"not a raw file")
        bad = self.folder / "bad.txt"
        bad.write_text(f"{SIZE} 0\n")
        result = subprocess.run([EXE, "-o", str(self.folder / "early.exr"), "--bad-pixels", str(bad),
                                 str(good), str(broken)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2, result.stderr)
        # The sidecar fails after the first decode, not after the broken second frame.
        self.assertIn("bad-pixel list coordinate is outside", result.stderr)
        self.assertNotIn("broken.dng", result.stderr)
        # Warnings gathered before a failure are still reported.
        self.assertIn("rawfuse: warning: good.dng: ignored DNG opcode", result.stderr)

        for args, error in ((["--noise", "-1,0.1"], "--noise coefficients"),
                            (["--noise", "0.1"], "comma-separated"),
                            (["--matrix", "custom:" + str(self.folder / "absent.txt")],
                             "could not open custom matrix")):
            with self.subTest(args=args):
                result = subprocess.run([EXE, "-o", str(self.folder / "x.exr"), *args, str(broken)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(error, result.stderr)

    def test_undated_frame_applies_dated_bad_pixels(self):
        frame = self.folder / "undated.dng"
        dng(frame, (1, 10), lambda y, x, c: 0.2, date=None)
        bad = self.folder / "bad.txt"
        bad.write_text("10 12 1700000000\n")
        report, rgb = self.fuse("undated", frame, "--bad-pixels", bad, "--aux")
        i = 12 * SIZE + 10
        self.assertEqual(rgb["flags"][i] & (8 | 16), 8 | 16)
        self.assertTrue(any("no capture time" in w for w in report["warnings"]))

    def test_dated_defects_use_the_capture_utc_offset(self):
        # Captured 2026-01-02 03:04:05 camera time; the defect dates from 00:00 UTC that day.
        bad = self.folder / "dated.txt"
        bad.write_text(f"10 12 {calendar.timegm((2026, 1, 2, 0, 0, 0))}\n")
        i = 12 * SIZE + 10

        def run(name, offset, *args, zone="UTC"):
            frame = self.folder / (name + ".dng")
            dng(frame, (1, 10), lambda y, x, c: 0.2, offset=offset)
            output = self.folder / (name + ".exr")
            result = subprocess.run([EXE, "--space", "camera", "--aux", "--bad-pixels", str(bad), *args,
                                     str(frame), "-o", str(output)], capture_output=True, text=True,
                                    env={**os.environ, "TZ": zone})
            self.assertEqual(result.returncode, 0, result.stderr)
            report = json.loads(output.with_suffix(".json").read_text())
            return report, bool(read_exr(output)[1]["flags"][i] & 16)

        # Without an offset the time is unknown wherever rawfuse runs: the entry applies.
        for zone in ("UTC", "Australia/Sydney"):
            report, excluded = run("plain-" + zone.replace("/", "-"), None, zone=zone)
            self.assertTrue(excluded)
            self.assertIsNone(report["inputs"][0]["capture_utc"])
            self.assertTrue(any("no capture time with a known UTC offset" in w for w in report["warnings"]))
        report, excluded = run("utc", None, "--camera-timezone", "+00:00", zone="Australia/Sydney")
        self.assertTrue(excluded)
        self.assertEqual(report["inputs"][0]["capture_time_source"], "camera_timezone")
        self.assertEqual(report["inputs"][0]["capture_utc"], calendar.timegm((2026, 1, 2, 3, 4, 5)))
        # At +10:00 the frame was taken the previous afternoon, before the defect.
        _, excluded = run("sydney", None, "--camera-timezone", "+10:00")
        self.assertFalse(excluded)
        # The file's own offset wins over the option.
        report, excluded = run("exif", b"+10:00\0", "--camera-timezone", "+00:00")
        self.assertFalse(excluded)
        self.assertEqual(report["inputs"][0]["capture_time_source"], "exif_offset")
        self.fuse("bad-zone", self.folder / "exif.dng", "--camera-timezone", "10", error="--camera-timezone expects")

    def test_frames_use_their_own_black_levels(self):
        # LibRaw measures some cameras' black level per frame (Canon CR3 varies by
        # 1 DN). Each frame must be fused against its own level.
        long, short = self.folder / "long.dng", self.folder / "short.dng"
        dng(long, (1, 10), lambda y, x, c: 0.2)
        dng(short, (1, 40), lambda y, x, c: 0.05, black_pattern=(4596,) * 4)
        report, rgb = self.fuse("own-black", long, short, "--no-deghost")
        # Both frames measure 2 per second; the geometric-mean scale is 0.05 s.
        for c in "RGB":
            close(rgb[c][1000], 0.1)
        self.assertTrue(any("levels differ between frames" in w for w in report["warnings"]))
        self.assertEqual([i["black_level"][0] for i in report["inputs"]], [4096, 4596])

    def test_later_group_without_contribution_mask_only_fills(self):
        base, nd = self.folder / "base.dng", self.folder / "nd.dng"
        dng(base, (1, 10), lambda y, x, c: 1 if x < 64 else 0.2)
        dng(nd, (1, 10), lambda y, x, c: 0.3 if x < 64 else 0.05)
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [
            {"name": "base", "frames": [base.name]}, {"name": "nd", "frames": [nd.name], "ev": 2}]}))
        report, rgb = self.fuse("fill", "--groups", manifest, "--aux")
        row = 40 * SIZE
        close(rgb["G"][row + 20], 1.2)  # The clipped base is filled from the filtered group.
        close(rgb["G"][row + 100], 0.2)  # Clean base sites are kept.
        self.assertEqual([rgb["source"][row + x] for x in (20, 100)], [2, 1])

        # A filtered value below the clipping bound contradicts it; the bound is kept.
        dng(nd, (1, 10), lambda y, x, c: 0.1 if x < 64 else 0.05)
        report, rgb = self.fuse("inconsistent", "--groups", manifest, "--aux", "--highlights", "bound",
                                "--force")
        self.assertEqual(rgb["flags"][row + 20] & (1 | 32), 1 | 32)
        self.assertEqual(rgb["source"][row + 20], 1)
        self.assertGreater(report["merge"]["inconsistent_fraction"], 0.4)

    def stick_brackets(self):
        # A blocked bracket whose blocker (a stick, columns 56-71) is excluded, and an
        # open bracket that sees the whole scene through a smooth additive glare veil.
        def scene(y, x):
            return 0.1 + 0.08 * (((x // 16) + (y // 16)) % 2)
        def veil(y, x):
            return 0.05 + 0.1 * x / SIZE
        blocked, open_ = self.folder / "blocked.dng", self.folder / "open.dng"
        dng(blocked, (1, 10), lambda y, x, c: 0 if 56 <= x < 72 else scene(y, x))
        dng(open_, (1, 10), lambda y, x, c: scene(y, x) + veil(y, x))
        def mask(name, selected):
            path = self.folder / name
            row = bytes(255 if selected(x) else 0 for x in range(SIZE))
            path.write_bytes(f"P5\n{SIZE} {SIZE}\n255\n".encode() + row * SIZE)
            return path.name
        return scene, veil, blocked, open_, mask

    def test_offset_match_removes_the_open_brackets_glare_veil(self):
        scene, veil, blocked, open_, mask = self.stick_brackets()
        manifest = self.folder / "groups.json"
        def write(match):
            group = {"name": "open", "frames": [open_.name],
                     "contribution_mask": mask("stick.pgm", lambda x: 56 <= x < 72)}
            if match:
                group["match"] = match
            manifest.write_text(json.dumps({"groups": [
                {"name": "blocked", "frames": [blocked.name],
                 "valid_mask": mask("outside.pgm", lambda x: not 56 <= x < 72)}, group]}))
            return manifest
        report, rgb = self.fuse("veil", "--groups", write(None), "--aux")
        match = report["groups"][1]["offset_match"]
        self.assertGreater(match["overlap_samples"], 10000)
        self.assertEqual(match["corrected_sites"], 16 * SIZE)
        close(match["median_offset"]["G"], veil(0, 64), 0.01)
        # Sites away from the texture's edges, where demosaicing reproduces it exactly.
        for y in (24, 40, 104):
            for x in (58, 61, 68):
                i = y * SIZE + x
                self.assertEqual(rgb["source"][i], 2)
                close(rgb["G"][i], scene(y, x), 0.01)
        self.assertEqual(match["status"], "corrected")
        report, rgb = self.fuse("veil-kept", "--groups", write("none"), "--aux")
        self.assertIsNone(report["groups"][1]["offset_match"])
        i = 40 * SIZE + 61
        close(rgb["G"][i], scene(40, 61) + veil(40, 61), 0.005)

    def test_offset_match_without_overlap_reports_uncorrected_sites(self):
        _, _, blocked, open_, mask = self.stick_brackets()
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [
            {"name": "left", "frames": [blocked.name], "valid_mask": mask("left.pgm", lambda x: x < SIZE // 2)},
            {"name": "right", "frames": [open_.name], "valid_mask": mask("right.pgm", lambda x: x >= SIZE // 2)}]}))
        output = self.folder / "disjoint.exr"
        result = subprocess.run([EXE, "--space", "camera", "--groups", str(manifest), "-o", str(output)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        match = json.loads(output.with_suffix(".json").read_text())["groups"][1]["offset_match"]
        self.assertEqual(match["status"], "no_overlap")
        self.assertEqual(match["corrected_sites"], 0)
        self.assertEqual(match["unresolved_sites"], SIZE * SIZE // 2)
        self.assertEqual(match["anchored"], {"R": False, "G": False, "B": False})
        self.assertIn("no precise overlap with earlier groups in R, G, B; 8,192 inserted sites not corrected",
                      result.stderr)

    def test_contribution_mask_group_also_fills_gaps_outside_the_mask(self):
        scene, veil, blocked, open_, mask = self.stick_brackets()
        manifest = self.folder / "groups.json"
        # The blocked bracket is also invalid in columns 100-109, outside the open
        # bracket's contribution mask; the open bracket fills them rather than leaving
        # them missing, and its glare veil is removed there too.
        manifest.write_text(json.dumps({"groups": [
            {"name": "blocked", "frames": [blocked.name],
             "valid_mask": mask("outside.pgm", lambda x: not (56 <= x < 72 or 100 <= x < 110))},
            {"name": "open", "frames": [open_.name], "contribution_mask": mask("stick.pgm", lambda x: 56 <= x < 72)}]}))
        report, rgb = self.fuse("gaps", "--groups", manifest, "--aux")
        i = 40 * SIZE + 104
        self.assertEqual(rgb["source"][i], 2)
        self.assertEqual(rgb["flags"][i] & 16, 0)
        close(rgb["G"][i], scene(40, 104), 0.01)
        self.assertEqual(report["merge"]["missing_fraction"], 0)
        # Clean blocked samples outside the mask are kept.
        self.assertEqual(rgb["source"][64 * SIZE + 20], 1)

    def test_alignment_measures_sub_pixel_shifts_within_and_between_groups(self):
        size = 256
        rng = random.Random(7)
        # Texture down to a 12-site wavelength, as lens blur leaves it. Finer detail,
        # near the quad grid's limit, biases a shift by a few hundredths of a site.
        waves = [(rng.uniform(0, 2 * math.pi), rng.uniform(12, 48), rng.uniform(0, math.pi), rng.uniform(0.3, 1))
                 for _ in range(12)]
        norm = sum(w[3] for w in waves)
        def scene(y, x):  # Rate at the (1, 10) exposure
            texture = sum(a * math.sin(2 * math.pi * (x * math.cos(t) + y * math.sin(t)) / l + p)
                          for p, l, t, a in waves) / norm
            return 0.1 * (1 + 0.8 * texture)
        def frame(name, exposure, scale, shift):
            path = self.folder / name
            dy, dx = shift
            dng(path, exposure, lambda y, x, c: scale * scene(y - dy, x - dx), width=size, height=size)
            return path.name
        # Group a: the (1, 10) frame is the reference; the long frame moved.
        a = [frame("a1.dng", (1, 40), 0.25, (0, 0)), frame("a2.dng", (1, 10), 1, (0, 0)),
             frame("a3.dng", (1, 3), 10 / 3, (0.3, -0.2))]
        # Group b, a second bracket after the camera was bumped.
        b = [frame("b1.dng", (1, 10), 1, (0.7, -1.3)), frame("b2.dng", (1, 3), 10 / 3, (0.7, -1.3))]
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [{"name": "a", "frames": a}, {"name": "b", "frames": b}]}))
        report, _ = self.fuse("aligned", "--groups", manifest, "--no-deghost", "--align", "report")
        frames = {(f["group"], f["input"]): f for f in report["alignment"]["frames"]}
        self.assertTrue(frames[("a", 2)]["reference"])
        expected = {("a", 1): (0, 0), ("a", 3): (0.3, -0.2), ("b", 1): (0.7, -1.3), ("b", 2): (0.7, -1.3)}
        for key, (dy, dx) in expected.items():
            f = frames[key]
            self.assertTrue(f["aligned"], key)
            self.assertIn("aligned_to", f, key)
            close(f["shift"]["rows"], dy, 0.02)
            close(f["shift"]["cols"], dx, 0.02)
            self.assertLess(f["max_shift"] - math.hypot(dy, dx), 0.05, key)
            close(f["gain_ev"], 0, 0.02)
        self.assertFalse(report["alignment"]["applied"])
        self.assertIn("the merge does not compensate motion yet", " ".join(report["warnings"]))
        report, _ = self.fuse("unaligned", "--groups", manifest)
        self.assertIsNone(report["alignment"])

    def test_nd_gate_counts_distinct_sites(self):
        # Repeated filtered frames of one small patch add observations, not coverage.
        base = self.folder / "base.dng"
        dng(base, (1, 10), lambda y, x, c: 0.2)

        def square(name, side):
            path = self.folder / name
            path.write_bytes(f"P5\n{SIZE} {SIZE}\n255\n".encode() + bytes(
                255 if y < side and x < side else 0 for y in range(SIZE) for x in range(SIZE)))
            return path.name

        def manifest(count, side):
            frames = []
            for k in range(count):
                frame = self.folder / f"nd{k}.dng"
                dng(frame, (1, 10), lambda y, x, c: 0.2 / 8)
                frames.append(frame.name)
            path = self.folder / "groups.json"
            path.write_text(json.dumps({"groups": [
                {"name": "clear", "frames": [base.name]},
                {"name": "nd", "frames": frames, "ev": "auto", "valid_mask": square(f"sq{side}.pgm", side),
                 "match": "none"}]}))
            return path

        self.fuse("one", "--groups", manifest(1, 32), error="256/512/256 distinct usable sites")
        self.fuse("nine", "--groups", manifest(9, 32), error="256/512/256 distinct usable sites")
        report, _ = self.fuse("wide", "--groups", manifest(2, 96))
        calibration = report["groups"][1]["calibration"]
        self.assertEqual(calibration["sites"], {"R": 2304, "G": 4608, "B": 2304})
        self.assertEqual(calibration["samples"], {"R": 4608, "G": 9216, "B": 4608})
        close(report["groups"][1]["ev"]["G"], 3, 0.01)

    def test_status_lines_and_quiet(self):
        frame = self.folder / "status.dng"
        dng(frame, (1, 10), lambda y, x, c: 0.2)
        output = self.folder / "status.exr"
        result = subprocess.run([EXE, str(frame), "-o", str(output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for line in ("rawfuse: read 1 frames in 1 group", "rawfuse: plain: 1 frame",
                     "rawfuse: noise: DNG NoiseProfile", "rawfuse: plain: supplies",
                     "rawfuse: result: 0% of sites missing", "rawfuse: wrote " + str(output)):
            self.assertIn(line, result.stderr)
        self.assertEqual(result.stdout, "")
        # The synthetic camera has no spectral data, so auto falls back to metadata;
        # that is reported, not printed under -q, and -q wins over -v.
        self.assertIn("rawfuse: colour: --matrix auto uses the metadata solver", result.stderr)
        report = json.loads(output.with_suffix(".json").read_text())
        self.assertEqual(report["matrix"]["method"], "auto")
        self.assertEqual(report["matrix"]["resolved"], "metadata")
        self.assertTrue(report["matrix"]["fallback_reason"])
        for extra in ([], ["-v"], ["--matrix", "metadata", "-v"]):
            result = subprocess.run([EXE, "-q", *extra, "--force", str(frame), "-o", str(output)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(all(line.startswith("rawfuse: warning:") for line in result.stderr.splitlines()),
                            result.stderr)

    def test_measured_exposure_ratio_replaces_rounded_exif(self):
        # The long frame's real exposure is 0.1 EV longer than its EXIF label.
        long, short = self.folder / "long.dng", self.folder / "short.dng"
        true_ratio = 4 * 2 ** 0.1
        dng(long, (1, 10), lambda y, x, c: 0.05 * true_ratio)
        dng(short, (1, 40), lambda y, x, c: 0.05)
        report, rgb = self.fuse("measured", long, short, "--no-deghost")
        step, = report["groups"][0]["exposure_steps"]
        self.assertTrue(step["applied"])
        close(step["exif_ratio"], 4)
        close(step["measured_ratio"], true_ratio, 1e-3)
        corrections = [i["exposure_correction_ev"] for i in report["inputs"]]
        close(corrections[0], 0.05, 1e-3)
        close(corrections[1], -0.05, 1e-3)
        # Both frames now agree: 0.05 at sqrt(0.025 · 0.1 / ratio) seconds, scaled to the 0.05 s mean.
        expected = 0.05 / math.sqrt(0.0025 / true_ratio) * 0.05
        for c in "RGB":
            close(rgb[c][1000], expected, 1e-3)

        report, _ = self.fuse("exif", long, short, "--no-deghost", "--exposures", "exif")
        self.assertEqual(report["groups"][0]["exposure_steps"], [])
        self.assertEqual(report["inputs"][0]["exposure_correction_ev"], 0)

        # A ratio far from EXIF is not trusted.
        dng(long, (1, 10), lambda y, x, c: 0.05 * 4 * 2 ** 0.5)
        report, _ = self.fuse("far", long, short, "--no-deghost")
        self.assertFalse(report["groups"][0]["exposure_steps"][0]["applied"])
        self.assertTrue(any("differs from EXIF by 0.50 EV" in w for w in report["warnings"]))

    def test_unmeasurable_step_borrows_the_camera_shutter_correction(self):
        true_ratio = 4 * 2 ** 0.1
        frames = {}
        for name, scale in (("plain", 1), ("filtered", 0.004)):
            long, short = self.folder / f"{name}-long.dng", self.folder / f"{name}-short.dng"
            dng(long, (1, 10), lambda y, x, c: 0.05 * true_ratio * scale)
            dng(short, (1, 40), lambda y, x, c: 0.05 * scale)
            frames[name] = [long.name, short.name]
        manifest = self.folder / "groups.json"
        manifest.write_text(json.dumps({"groups": [
            {"name": "plain", "frames": frames["plain"]},
            {"name": "filtered", "frames": frames["filtered"], "ev": 8}]}))
        report, _ = self.fuse("borrowed", "--groups", manifest, "--no-deghost")
        step, = report["groups"][1]["exposure_steps"]
        self.assertFalse(step["applied"])
        self.assertEqual(step["borrowed_from"], "plain")
        close(step["ratio"], true_ratio, 1e-3)
        self.assertFalse(any("could not be measured" in w for w in report["warnings"]))

    def test_noise_is_measured_without_noise_profile(self):
        shot, read = 2e-4, 1e-6
        generator = random.Random(7)
        width = 256
        def scene(scale):
            def value(y, x, c):
                signal = 0.002 * (0.85 / 0.002) ** (x / (width - 1)) * scale
                return signal + generator.gauss(0, math.sqrt(shot * signal + read))
            return value
        long, short = self.folder / "noisy-long.dng", self.folder / "noisy-short.dng"
        dng(long, (1, 10), scene(1), width=width, height=width, noise=None)
        dng(short, (1, 40), scene(0.25), width=width, height=width, noise=None)
        report, _ = self.fuse("noisy", long, short, "--no-deghost")
        fit, = report["noise_fits"]
        self.assertEqual(fit["read_noise"], "bracket")
        close(fit["S"], shot, 0.1)
        close(fit["O"], read, 0.3)
        self.assertEqual({i["noise_profile"]["source"] for i in report["inputs"]}, {"measured"})
        self.assertFalse(any("generic" in w for w in report["warnings"]))

        # A single frame cannot be measured and gets a generic profile with a warning.
        report, _ = self.fuse("generic", long)
        self.assertEqual(report["noise_profile"]["source"], "generic")
        self.assertTrue(any("generic profile" in w for w in report["warnings"]))


if __name__ == "__main__":
    unittest.main()

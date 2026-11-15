"""Regressions found during the full pipeline review; no external Python packages."""

import json
import itertools
from pathlib import Path
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
            {"name": "nd", "frames": [nd.name], "ev": 1, "valid_mask": mask.name}]}))
        matrix = self.folder / "red.txt"
        matrix.write_text("1 0 0\n1 0 0\n1 0 0\n")
        _, rgb = self.fuse("highlight", "--groups", manifest, "--aux", "--space", "aces",
                           "--wb", "custom:1,2,1", "--matrix", "custom:" + str(matrix))
        i = 10 * SIZE + 10
        close(rgb["R"][i], 1.6)
        self.assertEqual(rgb["source"][i], 65535)
        self.assertEqual(rgb["flags"][i] & (1 | 2), 1 | 2)

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
        # Clipping at 1/1000 s bounds the rate at 1000 per second; the later
        # 1/10 s bound (10 per second) is weaker and must not replace it.
        i = 10 * SIZE + 10
        close(rgb["R"][i], 1)
        self.assertEqual(rgb["source"][i], 1)
        self.assertEqual(rgb["flags"][i] & 1, 1)
        self.assertEqual(report["groups"][0]["lower_bound_sites"], SIZE * SIZE)

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

    def test_nd_all_updates_shadow_statistics(self):
        plain, nd = self.folder / "plain.dng", self.folder / "nd.dng"
        dng(plain, (1, 1000), lambda y, x, c: 0.001)
        dng(nd, (1, 10), lambda y, x, c: 0.05)
        report, _ = self.fuse("nd", plain, "--nd", nd, "--nd-ev", "1", "--nd-use", "all")
        self.assertEqual(report["merge"]["shadow_fraction"], 0)


if __name__ == "__main__":
    unittest.main()

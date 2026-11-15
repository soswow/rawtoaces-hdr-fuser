# rawfuse

`rawfuse` merges fixed-camera Bayer RAW exposures into a single linear, float32 OpenEXR. It merges the sensor samples before demosaicing and uses rawtoaces for the ACES AP0 colour transform. The default JSON sidecar records exposure metadata, the merge statistics and each input’s embedded DNG XMP.

## Build and install

Use a C++17 compiler, CMake 3.18+, pkg-config, LibRaw 0.21+, OpenImageIO, OpenEXR 3, nlohmann/json, and the rawtoaces 2.2 C++ utility library. Use the same OpenImageIO installation for rawtoaces and rawfuse. CTest uses Python 3.8+ and generates its own DNGs.

When rawtoaces is already installed, configure with its installation prefix on `CMAKE_PREFIX_PATH`. CMake also accepts `-DRAWTOACES_DIR=/path/to/the/package/config/directory`.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
cmake --install build --prefix "$PWD/build/local"
build/local/bin/rawfuse --version
```

### Build rawtoaces in a local prefix

The example assumes a neighbouring rawtoaces source checkout. On macOS with Homebrew, install the development dependencies and configure both projects against the resulting prefixes. The dependency build keeps its optional Python bindings and tests off. For spectral mode, install the rawtoaces spectral data with `RTA_INSTALL_DATABASE=ON` or point `RAWTOACES_DATA_PATH` at an existing rawtoaces database. The metadata path uses the camera's DNG tags.

```sh
brew install cmake ninja pkgconf libraw openimageio openexr nlohmann-json \
    eigen ceres-solver lensfun exiftool

prefix="$PWD/build/local"
packages="$(brew --prefix)"
cmake -S ../rawtoaces -B build/rawtoaces -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DCMAKE_PREFIX_PATH="$packages" \
    -DEigen3_DIR="$(brew --prefix eigen)/share/eigen3/cmake" \
    -DRTA_BUILD_PYTHON_BINDINGS=OFF -DRTA_BUILD_TESTS=OFF \
    -DRTA_INSTALL_DATABASE=OFF
cmake --build build/rawtoaces --parallel
cmake --install build/rawtoaces

cmake -S . -B build/fuser -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$prefix;$packages" -DCMAKE_INSTALL_PREFIX="$prefix"
cmake --build build/fuser --parallel
ctest --test-dir build/fuser --output-on-failure
cmake --install build/fuser
"$prefix/bin/rawfuse" --version
```

On other platforms, use the corresponding package prefixes for these dependencies. Pass `-DBUILD_TESTING=OFF` for a build without the Python test dependency. rawtoaces may invoke ExifTool for supplementary metadata.

## Usage

```sh
rawfuse shot-*.dng -o fused.exr
rawfuse shot-*.dng --nd nd-*.dng --nd-mask filter.png --nd-ev auto -o fused-nd.exr
rawfuse --space camera --normalize none --aux --write-cfa mosaic.exr shot-*.dng -o camera.exr
```

- Use a fixed camera and matching Bayer layouts, active areas, and black and white levels. DNG carries the richest metadata. Its exact `ActiveArea` defines the output and mask dimensions, including odd margins; `DefaultCrop` and orientation do not crop or rotate the sensor grid. Other LibRaw Bayer formats use the active area, levels and exposure metadata exposed by LibRaw.
- Plain filenames can arrive in any order. `--reference N` and `--wb frame:N` use the one-based order of the plain arguments. Auto reference selects the frame with the widest usable Bayer-quad coverage for motion handling. White balance defaults to the plain-bracket median. Use `--wb reference` for the reference frame or `--wb custom:R,G,B` for caller-specified multipliers.
- The default output is AP0 at the bracket’s geometric mean exposure. Its OpenEXR exposure tags describe that scale at ISO 100 and the normalization f-number. Camera RGB omits the final white balance and matrix and records them as metadata when available. `--normalize none` gives a rate per second at ISO 100 and the first supplied plain frame’s f-number. Linear RGB retains zero and negative values by default. Pass `--floor F` to set a positive minimum in output units for log-domain consumers.
- A full-frame ND mask uses nonzero luminance or RGB samples to mark the covered area. Auto calibration needs at least 2,000 usable overlapping samples per colour. Supply a measured neutral filter strength with `--nd-ev EV` if the captures lack overlap. Measured drift above 0.2 EV in the uncovered region warns with automatic EV and fails with manual EV. Use `--allow-nd-scene-drift` to proceed with manual EV and retain that warning in the report.
- An optional full-frame flat field has one gain per Bayer site or three RGB gains. A per-site dark model CSV takes zero-based active-area `row,col` and `offset_dn,rate_dn_per_s` in stored sensor units. The fuser fits offset and rate scales independently for the plain and ND sweeps, estimates residual uncertainty, and marks inconsistent sites missing. A text bad-pixel list uses active-area `col row [timestamp]` with a Unix timestamp in seconds. A dated entry applies to frames captured at or after that time; frames without a capture time receive every entry, with a warning. Sidecars and masks are checked after the first frame is decoded, before the others. Embedded DNG defect coordinates are translated to the active area.
- `--aux` adds per-channel propagated variance, contribution counts and flags. `--write-cfa` records the unscaled flat-fielded mosaic before missing-site and highlight repair. EXR variance models sensor noise and interpolation weights; it does not include uncertainty in ND strength, exposure metadata, white balance or flat-field calibration. Repairs introduce correlations that the variance planes do not fully model. Repaired or saturated values carry flags for interpretation. Non-finite CFA values and overflowing output variance cause an error before outputs are published.
- The output EXR, JSON report and optional mosaic use distinct paths. Use `--force` to replace existing files. Install and run the program on a filesystem that supports hard links for atomic publication of new files.

### Brackets with independent masks

Use `--groups` with a JSON manifest when an occluder, filter, or other change makes only part of a bracket useful. All frames share the camera's active Bayer grid. Mask images match that grid, in sensor orientation. A nonzero grayscale or RGB sample selects the site. Missing masks select the whole grid. Paths in the manifest are relative to its directory.

```json
{
  "groups": [
    {
      "name": "blocked",
      "frames": ["blocked-1.dng", "blocked-2.dng"],
      "valid_mask": "outside-blocker.png"
    },
    {
      "name": "open",
      "frames": ["open-1.dng", "open-2.dng"],
      "valid_mask": "reliable-open-samples.png",
      "contribution_mask": "blocker-and-support.png"
    },
    {
      "name": "filtered",
      "frames": ["nd-1.dng", "nd-2.dng"],
      "valid_mask": "filter-coverage.png",
      "contribution_mask": "highlight-patch.png",
      "ev": [9.7, 9.8, 9.9]
    }
  ]
}
```

```sh
rawfuse --groups brackets.json --aux --write-cfa fused-cfa.exr -o fused.exr
```

- `valid_mask` excludes occlusion, flare, unfiltered lenses and other unreliable samples from measurement and calibration. In the example, `outside-blocker.png` selects the usable area and `blocker-and-support.png` selects the area to replace. An area can be invalid in every group.
- `contribution_mask` limits where a bracket supplies the output. Samples outside it remain available for calibration. At each Bayer site, a clean measurement wins over a sample adjacent to saturation or a clipping lower bound. Within the same quality class, the later group wins, except that between two lower bounds the higher bound is kept. These are binary selections, with ordinary demosaic interpolation at boundaries.
- `ev` is optional attenuation in stops. Give a neutral number, an RGB array, or `"auto"`. Auto compares valid samples against the first group or the unattenuated group named by `calibration_reference`. An optional `calibration_mask` limits the overlap independently of the insertion patch. Supplying a calibration region or reference with manual EV checks that value against the captures. Group calibration samples only eligible filtered sites. A masked lens stays excluded.
- The first group supplies the reference metadata, white balance and geometric-mean exposure scale. The JSON report lists each group's source ID and selected-site counts. With `--aux`, the EXR stores a UINT `source` channel. IDs start at 1 in manifest order, 0 means no group, and 65535 means repair, highlight neutralization or demosaicing combined groups. Repairs retain donor provenance and quality flags. The optional CFA EXR records source IDs before any repair.
- Where no group allows a contribution, the output has `missing` and `uncovered` (bit 8) flags and a value of zero, or the requested `--floor`. An eligible clipped site is reported as a lower bound. To leave a flare-contaminated clipped halo unresolved, exclude it with the open bracket's validity mask.

The same manifest applies to a single fisheye or a packed multi-view frame. For example, make the filtered group's validity mask zero on the unfiltered half of a dual-fisheye image. The traditional positional and `--nd` syntax remains available for two-bracket jobs.

Read `flags` and `source` as UINT channels. OpenImageIO normalizes integer values
when converting them to float, so a whole-image float read must not be cast back
to integers to interpret these fields. Request the native `uint` type instead, for
example by passing `"uint"` as the format to OpenImageIO's `ImageInput.read_image`.

Run `rawfuse --help` for all switches. The [algorithm and output reference](docs/ALGORITHM.md) describes the processing steps, output schema and known limits.

DNG `CFAPattern` entries are colour codes. `CFAPlaneColor` assigns those colours to
planes, and the per-plane colour metadata is mapped to RGB, including white balance,
noise and the rawtoaces metadata matrix. For DNGs with a non-RGB plane order, use
`--matrix metadata` (the default), `auto`, or a custom RGB-to-AP0 matrix. Forced
`adobe` mode rejects that layout because LibRaw assumes green on plane 1 and can
supply a fallback identity matrix for those files.

## Licence

Apache License 2.0. See [LICENSE](LICENSE).

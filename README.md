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
rawfuse shot-*.CR3 -o fused.exr
rawfuse --groups brackets.json -o fused-nd.exr
rawfuse --space camera --normalize none --aux --write-cfa mosaic.exr shot-*.dng -o camera.exr
```

- Use a fixed camera and matching Bayer layouts and active areas. Each frame is fused against its own black and white levels: LibRaw measures some cameras' black level per frame (Canon CR3 varies by 1 DN), and larger differences are reported as a warning. DNG carries the richest metadata. Its exact `ActiveArea` defines the output and mask dimensions, including odd margins; `DefaultCrop` and orientation do not crop or rotate the sensor grid. Other LibRaw Bayer formats use the active area, levels and exposure metadata exposed by LibRaw.
- Relative exposures are measured from the frames themselves. Cameras record rounded shutter labels (a "1/125" exposure may really be 1/131 s), which can leave visible steps where the merge moves from one frame to the next. Each pair of neighbouring frames is compared on sites well exposed in both, and the bracket keeps the EXIF geometric mean as its absolute scale. A step with too little overlap borrows the correction measured for the same shutter settings in another group, or uses EXIF. Corrections beyond 0.25 EV are not trusted and are reported. `--exposures exif` uses the metadata as recorded.
- The noise model comes from the DNG `NoiseProfile`. Without one (CR3, NEF, ARW and most other formats), rawfuse measures it: read noise from the sensor's masked optical-black pixels where LibRaw exposes them, and shot noise from temporal differences between neighbouring frames. A single frame without a profile gets a generic profile and a warning. `--noise S,O` overrides both.
- Progress goes to stderr: what was read, how each bracket's exposures and noise were measured, each group's attenuation and contribution, the merge result and the files written. `-q` prints only warnings and errors.
- Positional frames can arrive in any order and form one unattenuated bracket. `--reference N` and `--wb frame:N` use the one-based order of those arguments, or of the first group's frames. Auto reference selects the frame with the widest usable Bayer-quad coverage for motion handling. White balance defaults to the first bracket's median. Use `--wb reference` for the reference frame or `--wb custom:R,G,B` for caller-specified multipliers.
- The default output is AP0 at the bracket’s geometric mean exposure. Its OpenEXR exposure tags describe that scale at ISO 100 and the normalization f-number. Camera RGB omits the final white balance and matrix and records them as metadata when available. `--normalize none` gives a rate per second at ISO 100 and the first group's first f-number. Linear RGB retains zero and negative values by default. Pass `--floor F` to set a positive minimum in output units for log-domain consumers.
- An optional full-frame flat field has one gain per Bayer site or three RGB gains. A per-site dark model CSV takes zero-based active-area `row,col` and `offset_dn,rate_dn_per_s` in stored sensor units. The fuser fits offset and rate scales independently for each group, estimates residual uncertainty, and marks inconsistent sites missing. A text bad-pixel list uses active-area `col row [timestamp]` with a Unix timestamp in seconds. A dated entry applies to frames captured at or after that time; frames without a capture time receive every entry, with a warning. Sidecars and masks are checked after the first frame is decoded, before the others. Embedded DNG defect coordinates are translated to the active area.
- `--aux` adds per-channel propagated variance, contribution counts, flags and group source IDs. `--write-cfa` records the unscaled flat-fielded mosaic before missing-site and highlight repair. EXR variance models sensor noise and interpolation weights; it does not include uncertainty in ND strength, exposure metadata, white balance or flat-field calibration. Repairs introduce correlations that the variance planes do not fully model. Repaired or saturated values carry flags for interpretation. Non-finite CFA values and overflowing output variance cause an error before outputs are published.
- The output EXR, JSON report and optional mosaic use distinct paths. Use `--force` to replace existing files. Install and run the program on a filesystem that supports hard links for atomic publication of new files.

### Several brackets: ND filters, occluders and masks

Use `--groups` with a JSON manifest for more than one bracket: an ND-filtered bracket for the sun, a bracket with the sun blocked by a hand, or any change that makes only part of a bracket useful. The simplest ND job needs no masks:

```json
{
  "groups": [
    {"name": "plain", "frames": ["plain-1.CR3", "plain-2.CR3", "plain-3.CR3"]},
    {"name": "nd", "frames": ["nd-1.CR3", "nd-2.CR3", "nd-3.CR3"], "ev": "auto"}
  ]
}
```

The ND bracket's strength is measured against the plain bracket, and its frames fill only what the plain bracket clipped. The fuller example below combines a blocked bracket, an open bracket and a filtered view. All frames share the camera's active Bayer grid. Mask images match that grid, in sensor orientation. A nonzero grayscale or RGB sample selects the site. Missing masks select the whole grid. Paths in the manifest are relative to its directory.

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
- At each Bayer site, a clean measurement wins over a sample adjacent to saturation, which wins over a clipping lower bound; between two lower bounds the higher bound is kept. A later group fills every valid site that earlier groups measured worse. With a `contribution_mask` it also replaces equal-quality samples inside the mask. These are binary selections, with ordinary demosaic interpolation at boundaries.
- Later groups are offset-matched to what earlier groups recorded. With the sun in view, an open bracket carries a veil of glare that a blocked bracket does not; its samples would paste a bright patch in place of the blocker. rawfuse measures the difference where both brackets are precise, interpolates it smoothly into the replaced area and subtracts it. Set `"match": "none"` on a group to keep its values as recorded. Glare structure that exists only inside the replaced area (a diffraction spike running along a stick that points at the sun) cannot be measured; move the blocker between two blocked brackets to avoid it.
- A measurement more than 10% below an earlier group's clipping bound contradicts it and is flagged `INCONSISTENT`. A filling group keeps the bound there; inside a contribution mask the later measurement is used, because blocked and open brackets legitimately see different glare.
- `ev` is optional attenuation in stops. Give a neutral number, an RGB array, or `"auto"`. Auto compares valid samples against the first group or the unattenuated group named by `calibration_reference`; it needs 2,000 sites per colour that the reference measures with an SNR of at least 10 and an attenuated frame records with a predicted SNR of at least 5, and otherwise fails with the exposure gap between the brackets. With the sun in view, calibrate against the open bracket rather than a blocked one: the filtered bracket's glare is attenuated with the sun, so it matches the open bracket's. Give the filtered group a `valid_mask` covering only what the filter covers; a filter smaller than a fisheye's field of view, its mount, or a cloth around it otherwise bias the measurement. An optional `calibration_mask` limits the overlap independently of the insertion patch. Supplying a calibration region or reference with manual EV checks that value against the captures. A masked lens stays excluded.
- `drift_mask` marks sites that the group's frames record without the attenuation, such as the unfiltered lens of a dual-fisheye camera. They are compared against the reference to measure scene change between the brackets. Drift above 0.2 EV warns with automatic `ev` and fails with a manual `ev`, unless `--allow-scene-drift` is given.
- The first group supplies the reference metadata, white balance and geometric-mean exposure scale. The JSON report lists each group's source ID and selected-site counts. With `--aux`, the EXR stores a UINT `source` channel. IDs start at 1 in manifest order, 0 means no group, and 65535 means repair, highlight neutralization or demosaicing combined groups. Repairs retain donor provenance and quality flags. The optional CFA EXR records source IDs before any repair.
- Where no group is valid, the output has `missing` and `uncovered` (bit 8) flags and a value of zero, or the requested `--floor`. An eligible clipped site is reported as a lower bound. A halo the open bracket clips is filled by a later filtered group that is valid there; to leave it unresolved, exclude it from every group's validity mask.

The same manifest applies to a single fisheye or a packed multi-view frame. For example, make the filtered group's validity mask zero on the unfiltered half of a dual-fisheye image, and give that half as its `drift_mask`.

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

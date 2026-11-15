# rawfuse algorithm and output reference

This document describes what `rawfuse` does to its inputs, what it writes, and where its limits are. The [README](../README.md) covers building and the command line.

`rawfuse` reads unprocessed Bayer mosaics and exposure metadata. It merges the frames in camera space at each sensor site, applies the sensor gains once, demosaics the float mosaic once, and optionally converts to linear ACES AP0 with rawtoaces. The same sensor site must view the same scene point in every frame. The pipeline contains no camera-model table, filename pattern or lens geometry. Alignment, geometry, flare removal and tonemapping belong to the caller.

## Inputs

- One or more plain Bayer RAWs, and optionally a second bracket taken through a neutral density (ND) filter, or ordered `--groups` brackets (see below).
- All frames come from the same camera and match in size, CFA pattern, black and white level, make and model.
- Optional caller files: an ND or validity mask, a flat-field image, a per-site dark model CSV and a bad-pixel list. All use the active-area grid.

### DNG metadata

`rawfuse` reads BlackLevel, WhiteLevel, CFAPattern, CFAPlaneColor, exact Exif exposure time, ISO, f-number, AsShotNeutral, NoiseProfile, ImageDescription, XMP, BaselineExposure and DNG opcodes with its own TIFF/DNG tag walker, because LibRaw does not expose all of them reliably. ImageDescription and XMP are preserved in the JSON report. Opcode headers and parameters are big-endian regardless of TIFF byte order. OpcodeList1 `FixBadPixelsList` supplies defective points and rectangles, and OpcodeList2 `GainMap` supplies phase- and region-specific gains in normalized whole-image coordinates. Other opcodes produce a warning and are not applied. `BaselineExposure` is reported, not applied.

For DNG the exact `ActiveArea` (or the whole raw image when absent) defines the output grid, even when LibRaw rounds odd margins. `CFAPattern` entries are TIFF/EP colour codes (0 red, 1 green, 2 blue), not plane indices. `CFAPlaneColor` assigns those codes to planes and reorders per-plane tags. White balance, noise and matrix columns are converted between native plane order and RGB at the rawtoaces seam. For RGB plane order the pattern is cross-checked against LibRaw. LibRaw assumes green on plane 1, so it is not consulted for other plane orders.

### Dark model CSV

Columns are zero-based active-area `row,col` and `offset_dn,rate_dn_per_s` in stored DN above black. The nominal dark signal of a frame is `offset_dn + rate_dn_per_s × t`, with `t` the exposure in seconds.

## Notation

- Site `s` is one Bayer sample (row, col). Its colour `c(s)` is R, G or B (both greens are G).
- Frame `k` has raw sample `raw_k(s)`, black `b`, white `w` and full scale `F = w − b`.
- The net fraction is `x_k(s) = (raw_k(s) − b − d_k(s)) / F`, with `d_k` the dark signal (0 unless a dark model lists `s`). Negative values are kept because they are unbiased noise.
- The effective exposure is `E_k,c = t_k · (ISO_k / 100) · (N_ref / N_k)² · τ_k,c`, with `t` the exposure time, `N` the f-number, `N_ref` the first plain frame's f-number and `τ` the ND transmission (1 for plain frames).
- The rate is `L_k(s) = x_k(s) / E_k,c(s)`: what a 1 s, ISO 100 exposure would have recorded, as a fraction of full scale. The merge estimates `L(s)`.
- Noise follows `var(x) = S·x + O` from the DNG NoiseProfile, so `var(L_k) = (S·x + O + extra_k(s)) / E_k²`.

## Processing

### 1. Read and validate

LibRaw opens each frame and supplies the unprocessed mosaic. LibRaw's own processing (`dcraw_process`) is not used. Black comes from the DNG `BlackLevel` (per CFA position if present) or LibRaw, and white from the DNG `WhiteLevel` or LibRaw. For DNGs the tags are preferred. All frames must agree on size, CFA, levels, make and model. Plain frames are sorted by effective exposure. Grid-dependent caller files are loaded after the first frame is decoded and before the others, so a bad sidecar fails early.

### 2. Clipped, usable and defective sites

- A site is **clipped** if `raw ≥ w − clip_margin·F`, then dilated by `--clip-dilate` sites in every direction to cover blooming.
- A site is **usable** if it is clear of direct or dilated clipping and `x ≤ high`. If only near-white unsaturated readings exist, the shortest is a fallback. If nearby clipping excludes every usable sample, the shortest non-saturated reading is retained and flagged `ADJACENT_CLIP_FALLBACK`.
- **Defective sites** from the bad-pixel list, embedded DNG defects and the dark model start as missing. Dated bad-pixel entries apply to frames captured at or after that time. Frames without a capture time receive every entry, with a warning.
- **Dark model fit.** The CSV is fitted against the median black-subtracted DN of neighbouring same-phase sites, separately for plain and ND sweeps. The fit needs at least eight informative sites. Offset and rate scales are fitted independently when the data constrain both, or as one shared scale when their effects are indistinguishable. Robust weights account for scene texture, sensor noise and model outliers. At an accepted site `d_k = a·offset_dn + b·rate_dn_per_s·t_k` and `extra_k = S·max(0, d_k/F) + (σ_m·d_k/F)²`. The scatter `σ_m` is estimated from the residuals, with `--dark-model-error` as its lower limit. Sites with strong local or repeated disagreement stay missing and are repaired from same-colour neighbours.
- Rectangles in `FixBadPixelsList` stay missing.

### 3. Merge the plain frames

1. **Seed.** Triangular weights inside `[low, high]`: `w1 = max(0, min(x − low, high − x))`, and `seed = Σ w1·L_k / Σ w1`. If no frame lies inside the interval, take the longest usable frame's `L`.
2. **Inverse-variance merge.** The expected fraction is `x̂_k = max(seed, 0)·E_k`, the weight is `w_k = E_k² / (S·x̂_k + O + extra_k)` for usable frames, and `L = Σ w_k L_k / Σ w_k` with `var(L) = 1 / Σ w_k`. The weight uses the predicted signal rather than the noisy sample, which avoids bias.
3. **Motion rejection on 2 × 2 quads**, so all four colours of a quad come from the same frames and motion does not turn into colour fringes. The quad value `Q_k` is the mean of its four sites' `L_k`, defined where all four are usable and at least `low`. The reference is `--reference`, or the frame with the most quads defined. A quad is moving if any frame has `|ln Q_k − ln Q_ref| > ln 2 · deghost_ev`. The moving mask is dilated by one quad, and moving quads take the reference frame alone where it is usable. This needs at least three plain frames. No registration or deghosting between independent sweeps is attempted.
4. **Nothing usable at a site.** Use the shortest near-white fallback if one exists. Otherwise keep the shortest unsaturated reading from the clip dilation zone, flagged `ADJACENT_CLIP_FALLBACK`. For direct saturation set `L = F / E_min` with `E_min` the shortest directly saturated exposure, flagged `LOWER_BOUND`. Mark `MISSING` when the site has neither a measurement nor a bound.
5. **Statistics** go into the report: the fractions of unsupported, shadow-only, lower-bound, missing and moving sites.

### 4. ND frames (optional)

ND frames go through the same steps with their own transmission `τ`. Only sites inside `--nd-mask` count.

**Measuring the strength (`--nd-ev auto`).** Per colour, samples are ND-frame sites inside the mask that are usable, at least `nd_min_signal` (bright, so glare in dark areas does not count) and at most `high`, where the plain merge is a real measurement with relative σ below 5%. Each sample gives `log2 τ = log2((x_k / (t_k·ISO_k/100·(N_ref/N_k)²)) / L_plain)`. The estimate is the median, with samples further than `max(4 MAD, 0.1 EV)` from it dropped and the median taken again. At least 2,000 samples per colour are required, otherwise the run fails with a suggestion. Large images are sampled on an odd-stride grid across all Bayer phases, with sparsely represented colours retried at full resolution.

- **Spatial check.** In 64 × 64 blocks with at least 50 samples, a block whose median differs from the global value by more than 0.25 EV (filter edge, glare, a moved cloud) is removed from the ND mask and listed in the report.
- **Scene-drift guard.** The median exposure-normalized change outside the mask is estimated when at least 100 samples support it. Drift above 0.2 EV warns with automatic EV and fails with manual EV unless `--allow-nd-scene-drift` is given.

**Given strength (`--nd-ev 9.8`).** `τ = 2^−9.8` in every colour. Measurement still runs when it can, and values differing by more than 0.2 EV warn.

**Merging.** With `--nd-use fill` (default), inside the mask where the plain result is a lower bound the usable ND frames are merged and take that value, but never below the plain lower bound (if lower by more than 10%, the bound is kept and flagged `ND_INCONSISTENT`). Sites are flagged `FROM_ND`, and sites still clipped in every ND frame keep a lower bound from the shortest ND frame. With `--nd-use all`, ND frames join the weighted merge at every site inside the mask.

### 5. Ordered input groups (`--groups`)

The manifest format is in the [README](../README.md). Validity masks apply before clipping dilation, dark-model fitting, motion rejection and ND measurement. Contribution masks select output sites without removing a group's calibration samples. Each group is merged in camera space and multiplied by its per-colour attenuation correction. At each site a clean sample is preferred over an adjacent-clipping fallback, which is preferred over a lower bound. The later group wins within the same quality class, except between two lower bounds, where the higher (more informative) bound is kept. The mosaics of groups no later calibration needs are released after each merge. Missing eligible defects may be repaired, but sites with no eligible contribution stay uncovered.

Automatic attenuation compares valid overlap against an unattenuated reference group, independently of output insertion. Calibration reports suspect 64 × 64 blocks as `blocks_flagged`. It does not remove them from contribution masks, and it makes no uncovered-area drift measurement, because invalid parts of a group need not be comparable. The caller must provide an appropriate calibration region and consistent captures.

The first group determines the reference metadata, white balance and normalization.

### 6. Flat field

`--vignette dng` (default) applies every `GainMap` in OpcodeList2 to the merged sites it covers. The Bayer sample uses the pixel centre, `v = (row + 0.5) / image_height` and `u = (col + 0.5) / image_width`, matching [Adobe's GainMap reference implementation](https://android.googlesource.com/platform/external/dng_sdk/+/refs/heads/android14-prebuilt-test/source/dng_gain_map.cpp). Gains are interpolated bilinearly with edge clamping, and `L`, lower bounds and variance are scaled by the gain and its square. A DNG without OpcodeList2 gets no flat field, and the report says so. `--vignette FILE` takes an EXR of the full frame with three channels (a site takes its colour's channel) or one channel per site. `none` skips it.

### 7. Highlights and missing sites

These run after the flat field so colour comparisons use flat-fielded values.

- **Lower-bound sites** (`--highlights neutral`): in white-balanced camera space, a quad with any lower-bound site gets each of those sites raised to at least the quad's brightest white-balanced value, so fully clipped areas come out neutral instead of magenta. `bound` leaves the raw bounds.
- **Missing sites** are filled with the variance-weighted mean of the nearest same-colour sites with data (R and B: the four sites two away horizontally and vertically, then the four diagonals. G: the four diagonal G neighbours) and flagged `FILLED`.

Repairs keep `MISSING` as a record of absent measurements and inherit quality flags and group source IDs from the donors actually used. Sites excluded by every contribution mask remain `MISSING | UNCOVERED` and provide no interpolation taps to their neighbours.

### 8. Demosaic, once

Demosaicing runs on the float mosaic, which can span more than 20 EV and is never converted to integers. `bilinear` (default) takes each missing colour as the mean of its nearest same-colour neighbours. It is linear, so it commutes with exposure scaling. `mhc` uses Malvar–He–Cutler 5 × 5 gradient-corrected filters, falling back to bilinear where the correction exceeds half the bilinear value or any site in the window is a lower bound, uncovered or an unfilled hole. Uncovered and unfilled taps are omitted. Image borders are mirrored.

### 9. Colour

With `--space aces`, white-balanced camera RGB is multiplied by the rawtoaces matrix to give AP0. `rawfuse` configures an `rta::util::ImageConverter` from the reference DNG to obtain the white-balance multipliers and the 3 × 3 camera-to-AP0 matrix, then applies `out_i = Σ_j matrix[i][j] · in_j` itself, so variance propagates through the combined demosaic and colour weights. The multipliers are normalized so the smallest of R, G, B is 1. For a white balance other than the reference frame's (`median`, `custom`), rawtoaces is configured with `WBMethod::Custom`, because it derives the illuminant, and so the matrix, from the white balance. `--matrix` is passed through to rawtoaces' `matrix_method`, and `custom:FILE` holds nine numbers (camera RGB to AP0). With `--space camera`, camera RGB is written without white balance or matrix, and both are recorded in the report and EXR header. Forced `--matrix adobe` rejects DNGs with a non-RGB plane order, because LibRaw may supply an identity fallback matrix for them. rawtoaces' lens correction is not used.

### 10. Scale and write

- `--normalize geomean` multiplies `L` by `E_ref`, the geometric mean of the plain frames' `t·ISO/100·(N_ref/N)²`, so values resemble a middle frame's. `none` writes `L` itself, a rate per second at ISO 100.
- Signed linear output is preserved. `--floor F` applies a positive minimum to every channel after the flat field, colour and normalization, and counts affected pixels.
- The EXR is float32 (half overflows for bright light sources), ZIP-compressed, with channels R, G, B at the full active-area size and `chromaticities` set to AP0 for `--space aces`.
- Non-finite merged CFA signal or variance, and non-finite output RGB or variance, abort before any output is published. Outputs are staged and published atomically with hard links.

## Outputs

### EXR attributes

| Attribute | Content |
|---|---|
| `rawfuse:version` | program version; its presence marks a fused file |
| `rawfuse:space` | `aces-ap0` or `camera` |
| `rawfuse:exposureReference` | `E_ref` (seconds at ISO 100, f/`N_ref`); 1 with `--normalize none` |
| `rawfuse:whiteBalance` | normalized multipliers applied (R G B) |
| `rawfuse:matrix` | the 3 × 3 matrix applied (or to apply) |
| `rawfuse:ndEV` | R G B stops, or absent |
| `rawfuse:report` | the report's file name |
| `rawfuse:referenceFNumber` | aperture used to normalize effective exposure |
| `capDate` | capture date of the reference frame |
| `expTime`, `isoSpeed`, `aperture` | geometric-mean equivalent exposure time, ISO 100, and normalization aperture; omitted with `--normalize none` |

A non-DNG capture date is LibRaw's local wall-clock timestamp, formatted as local time.

### Auxiliary channels (`--aux`)

- `var.R`, `var.G`, `var.B`: variance of the output values.
- `frames`: number of frames that contributed (half float).
- `flags` (UINT bit field): bit 0 `LOWER_BOUND`, 1 `FROM_ND`, 2 `MOVING`, 3 `FILLED`, 4 `MISSING`, 5 `ND_INCONSISTENT`, 6 `FLOORED`, 7 `ADJACENT_CLIP_FALLBACK`, 8 `UNCOVERED`. Repair, highlight and demosaic donors are included.
- `source` (UINT, with `--groups`): 1-based group ID in manifest order, 0 for no source, 65535 for mixed sources. RGB source IDs include repair, highlight and demosaic dependencies.

`--write-cfa PATH` writes a one-channel EXR of the merged mosaic after the flat field and before repair and highlight handling, with `var` and `flags` channels and the CFA pattern as an attribute. Its source IDs describe selection before repair.

OpenImageIO normalizes integer channels to [0, 1] when asked for floats, so read `flags` and `source` in their native UINT type.

### JSON report

The report (`OUT.json` by default) records, among other things: the version, each input's path, role, exposure, ISO, f-number, AsShotNeutral, ImageDescription, XMP and defect counts; the reference frame; the white balance and matrix used; the exposure reference; the dark model fit; the noise profile; the flat field; ND calibration results; group source IDs, masks and selected-site counts; merge statistics; and warnings.

## Errors and warnings

Exit code 2 with a one-line reason: frames of different size, CFA pattern, levels or camera; a non-Bayer or linear raw; a frame without an exposure time; the same file given twice, including through a hard link; an existing output without `--force`; `--nd-ev auto` without enough overlap; manual ND EV with measured drift over 0.2 EV.

Warnings go to stderr and the report: every exposure being an exact 1/N (nominal shutter names), differing ISO between frames, a missing NoiseProfile (defaults are used), an unimplemented opcode, AsShotNeutral spread across frames, nonzero BaselineExposure and differing embedded GainMaps (the reference frame's maps are used). Warnings collected before a failure are printed with the error.

## Program structure

`rawfuse_pipeline` is a static library with a thin command-line front end in `src/main.cpp`.

| Module | Job |
|---|---|
| `options` | command-line parsing and validation |
| `metadata` | TIFF/DNG tag walker and opcode parsing |
| `raw_frame` | LibRaw reading: mosaic, geometry, levels, CFA, exposure |
| `dark_model` | dark-model CSV, fitting, per-site dark signal and extra variance |
| `merge` | exposure merge and motion handling |
| `nd` | ND calibration and merge |
| `groups` | ordered groups and masks |
| `finish` | flat field, highlights, repair and demosaic |
| `colour` | white balance and matrix through rawtoaces |
| `writer` | EXR and JSON output |
| `pipeline` | orchestration |

## Known limits

- Variance covers the sensor noise model and interpolation weights. It excludes uncertainty in shutter calibration, ND strength, colour transforms and flat-field fits. Repaired sites share donors, so their covariance is not fully modelled. Interpret lower bounds and repairs through the flags.
- A shared flat field is applied after fusion. Captures whose physical vignetting changes between frames need a deliberate per-frame calibration, not an average of maps.
- Group calibration's `blocks_flagged` are diagnostic and do not erase contribution regions, unlike the legacy ND path.
- Automatic ND strength needs enough overlap between the plain and ND exposures. Without it, supply a measured `--nd-ev`.
- There is no registration or deghosting between independent sweeps. The caller must supply a fixed camera and a static scene.
- The CLI is the supported interface. `rawfuse_pipeline` is a build target, but there is no installed, versioned C++ package or stable ABI. Its low-level stages assume validated input, and `run()` sets process-wide OpenEXR/OIIO thread limits (restored on return).
- Peak memory is dominated by decoding all frames of a group at once, roughly 50 MB per full-resolution frame.
- Development and testing so far use macOS with synthetic Bayer DNGs. Linux and Windows CI and native CR3 coverage are not yet in place.

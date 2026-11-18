# rawfuse algorithm and output reference

This document describes what `rawfuse` does to its inputs, what it writes, and where its limits are. The [README](../README.md) covers building and the command line.

`rawfuse` reads unprocessed Bayer mosaics and exposure metadata. It merges the frames in camera space at each sensor site, applies the sensor gains once, demosaics the float mosaic once, and optionally converts to linear ACES AP0 with rawtoaces. The same sensor site must view the same scene point in every frame. The pipeline contains no camera-model table, filename pattern or lens geometry. Alignment, geometry, flare removal and tonemapping belong to the caller.

## Inputs

- One bracket of Bayer RAWs on the command line, or ordered `--groups` brackets (see below), for example a plain bracket and one taken through a neutral density (ND) filter. Positional frames form a single unattenuated group named `plain`; everything below applies to it as a one-group manifest.
- All frames come from the same camera and match in size, CFA pattern, make and model. Black and white levels may differ per frame.
- Optional caller files: group masks, a flat-field image, a per-site dark model CSV and a bad-pixel list. All use the active-area grid.

### DNG metadata

`rawfuse` reads BlackLevel, WhiteLevel, CFAPattern, CFAPlaneColor, exact Exif exposure time, ISO, f-number, AsShotNeutral, NoiseProfile, ImageDescription, XMP, BaselineExposure and DNG opcodes with its own TIFF/DNG tag walker, because LibRaw does not expose all of them reliably. ImageDescription and XMP are preserved in the JSON report. Opcode headers and parameters are big-endian regardless of TIFF byte order. OpcodeList1 `FixBadPixelsList` supplies defective points and rectangles, and OpcodeList2 `GainMap` supplies phase- and region-specific gains in normalized whole-image coordinates. Other opcodes produce a warning and are not applied. `BaselineExposure` is reported, not applied.

For DNG the exact `ActiveArea` (or the whole raw image when absent) defines the output grid, even when LibRaw rounds odd margins. `CFAPattern` entries are TIFF/EP colour codes (0 red, 1 green, 2 blue), not plane indices. `CFAPlaneColor` assigns those codes to planes and reorders per-plane tags. White balance, noise and matrix columns are converted between native plane order and RGB at the rawtoaces seam. For RGB plane order the pattern is cross-checked against LibRaw. LibRaw assumes green on plane 1, so it is not consulted for other plane orders.

### Dark model CSV

Columns are zero-based active-area `row,col` and `offset_dn,rate_dn_per_s` in stored DN above black. The nominal dark signal of a frame is `offset_dn + rate_dn_per_s × t`, with `t` the exposure in seconds.

## Notation

- Site `s` is one Bayer sample (row, col). Its colour `c(s)` is R, G or B (both greens are G).
- Frame `k` has raw sample `raw_k(s)`, its own black `b_k`, white `w_k` and full scale `F_k = w_k − b_k`.
- The net fraction is `x_k(s) = (raw_k(s) − b_k − d_k(s)) / F_k`, with `d_k` the dark signal (0 unless a dark model lists `s`). Negative values are kept because they are unbiased noise.
- The EXIF exposure is `t_k · (ISO_k / 100) · (N_ref / N_k)²`, with `t` the exposure time, `N` the f-number and `N_ref` the first group's f-number. The effective exposure `E_k` is the measured version of it (step 3). A group's attenuation `τ_c = 2^−ev_c` multiplies `E_k` for its frames.
- The rate is `L_k(s) = x_k(s) / E_k`: what a 1 s, ISO 100 exposure would have recorded, as a fraction of full scale. The merge estimates `L(s)`.
- Noise follows `var(x) = S·x + O` (step 4), so `var(L_k) = (S·x + O + extra_k(s)) / E_k²`.

## Processing

### 1. Read and validate

LibRaw opens each frame and supplies the unprocessed mosaic. LibRaw's own processing (`dcraw_process`) is not used. Black comes from the DNG `BlackLevel` (per CFA position if present) or LibRaw, and white from the DNG `WhiteLevel` or LibRaw. For DNGs the tags are preferred. All frames must agree on size, CFA, make and model. Levels may differ: LibRaw measures some cameras' black level from each frame's masked pixels (Canon CR3 varies by 1 DN), so every sample uses its own frame's levels, and differences above 0.1% of the signal range produce a warning. Each group's frames are sorted by exposure. Where LibRaw reports masked (optical-black) areas outside the active grid, their spread is measured per frame: a MAD start, then a 4σ-clipped standard deviation around each Bayer phase's median, which stays accurate for integer data. Grid-dependent caller files are loaded after the first frame is decoded and before the others, so a bad sidecar fails early.

### 2. Clipped, usable and defective sites

- A site is **clipped** if `raw ≥ w − clip_margin·F`, then dilated by `--clip-dilate` sites in every direction to cover blooming.
- A site is **usable** if it is clear of direct or dilated clipping and `x ≤ high`. If only near-white unsaturated readings exist, the shortest is a fallback. If nearby clipping excludes every usable sample, the shortest non-saturated reading is retained and flagged `ADJACENT_CLIP_FALLBACK`.
- **Defective sites** from the bad-pixel list, embedded DNG defects and the dark model start as missing. Dated bad-pixel entries apply to frames captured at or after that time. Frames without a capture time receive every entry, with a warning.
- **Dark model fit.** The CSV is fitted against the median black-subtracted DN of neighbouring same-phase sites, separately for each group, after steps 3 and 4. The fit needs at least eight informative sites. Offset and rate scales are fitted independently when the data constrain both, or as one shared scale when their effects are indistinguishable. Robust weights account for scene texture, sensor noise and model outliers. At an accepted site `d_k = a·offset_dn + b·rate_dn_per_s·t_k` and `extra_k = S·max(0, d_k/F) + (σ_m·d_k/F)²`. The scatter `σ_m` is estimated from the residuals, with `--dark-model-error` as its lower limit. Sites with strong local or repeated disagreement stay missing and are repaired from same-colour neighbours.
- Rectangles in `FixBadPixelsList` stay missing.

### 3. Measured exposures

Shutter times in EXIF are rounded labels, often by 0.05–0.1 EV per frame, and the error repeats at each setting. For each pair of neighbouring frames in a group, sites usable in both, without a dark-model correction and at least `max(0.03, low)` in the shorter frame give `x_long / x_short`. The median of at least 2,000 such ratios replaces the EXIF ratio of that step, unless it differs by more than 0.25 EV (reported as a warning). Large frames are sampled on an odd-stride grid, retried at full resolution when that leaves too few ratios. A step without a usable measurement takes the correction (measured / EXIF ratio) of a measured step with the same two exposure times and EXIF ratio in another group, typically the plain bracket for an ND bracket whose short frames are nearly black; otherwise it keeps the EXIF ratio. The chained ratios are scaled so the group's geometric-mean exposure equals the EXIF one, so absolute levels still follow the metadata. `--exposures exif` skips this step. The report lists each step's EXIF, measured and used ratio and each frame's correction.

### 4. Noise model

A DNG `NoiseProfile` (or `--noise S,O`) is used as given. Otherwise the model is measured from the bracket. For neighbouring frames at one ISO with `r = E_b / E_a`, the temporal difference `x_b − r·x_a` cancels the scene and has variance `S·x_b·(1 + r) + O_b + r²·O_a`. The difference is evaluated in 8 × 8 tiles of the sampling grid, each giving a 4σ-clipped variance and mean level. Motion, misregistration and texture only add variance, so estimates use the 10th percentile of the tile values, divided by its expected bias for Gaussian noise (the Wilson–Hilferty chi-squared quantile).

- With optical black in every frame, `O` per frame is the masked-pixel variance in normalized units, and `S` is the corrected 10th percentile of per-tile `(V − O_b − r²·O_a) / (x_b·(1 + r))` over tiles where shot noise is at least twice the read term.
- Otherwise the corrected 10th-percentile variance of each pair's tiles in 24 logarithmic level bins (at least 20 tiles each) feeds a relative-weight least-squares fit of `S` and `O`.
- A frame at an ISO without a fit scales the nearest ISO's fit in proportion. With no fit at all, a generic `S = 3·10⁻⁵ · ISO/100`, `O = 2·10⁻⁸ · ISO/100` (or the optical-black `O`) is used with a warning.

The noise model sets the merge weights, the variance output and which samples are reliable enough for attenuation calibration.

### 5. Merge each group

1. **Seed.** Triangular weights inside `[low, high]`: `w1 = max(0, min(x − low, high − x))`, and `seed = Σ w1·L_k / Σ w1`. If no frame lies inside the interval, take the longest usable frame's `L`.
2. **Inverse-variance merge.** The expected fraction is `x̂_k = max(seed, 0)·E_k`, the weight is `w_k = E_k² / (S·x̂_k + O + extra_k)` for usable frames, and `L = Σ w_k L_k / Σ w_k` with `var(L) = 1 / Σ w_k`. The weight uses the predicted signal rather than the noisy sample, which avoids bias.
3. **Motion rejection on 2 × 2 quads**, so all four colours of a quad come from the same frames and motion does not turn into colour fringes. The quad value `Q_k` is the mean of its four sites' `L_k`, defined where all four are usable and at least `low`. The reference is `--reference` (first group only), or the frame with the most quads defined. A quad is moving if any frame has `|ln Q_k − ln Q_ref| > ln 2 · deghost_ev`. The moving mask is dilated by one quad, and moving quads take the reference frame alone where it is usable. This needs at least three frames in the group. No registration or deghosting between independent sweeps is attempted.
4. **Nothing usable at a site.** Use the shortest near-white fallback if one exists. Otherwise keep the shortest unsaturated reading from the clip dilation zone, flagged `ADJACENT_CLIP_FALLBACK`. For direct saturation set `L = F / E_min` with `E_min` the shortest directly saturated exposure, flagged `LOWER_BOUND`. Mark `MISSING` when the site has neither a measurement nor a bound.
5. **Statistics** go into the report: the fractions of unsupported, shadow-only, lower-bound, missing and moving sites.

### 6. Attenuation calibration

A group with `ev: "auto"`, a `calibration_mask`, a `drift_mask` or a `calibration_reference` is compared with its unattenuated reference group (the first group by default), merged on the overlap.

**Measuring the strength.** Per colour, samples are sites inside both groups' validity and the calibration mask where the reference merge is a real measurement (not a bound, fallback, missing or moving) with SNR of at least 10, and where an attenuated frame reads the site unclipped. Each sample gives `ev = log2(E_k · L_ref / x_k)`. Selecting samples by their own noisy reading would favour those that happened to read high, so selection runs twice: a first estimate from samples whose reading has SNR of at least 5, then the final one from samples whose *predicted* reading `E_k · L_ref · 2^−ev` lies in `[low, 0.8·high]` with SNR of at least 5. A strong filter can leave only low-signal overlap (a 14 EV filter's longest frame may be darker than the reference's shortest), and thousands of such samples still give the median to about 0.01 EV. The estimate is the median, with samples further than `max(4 MAD, 0.1 EV)` from it dropped and the median taken again. At least 2,000 samples per colour are required. Otherwise the run fails, giving the counts and, when one colour was measured, how far the longest attenuated frame (at that strength) is from the shortest reference frame. Large images are sampled on an odd-stride grid across all Bayer phases, with sparsely represented colours retried at full resolution.

- **Spatial check.** In 64 × 64 blocks with at least 50 samples, blocks whose median differs from the group's ev by more than 0.25 EV (filter edge, glare, a moved cloud) are listed as `blocks_flagged` with a warning. They are diagnostic and do not change the composition.
- **Scene drift.** On `drift_mask` sites, which the group's frames record without attenuation, unclipped readings with SNR of at least 5 are compared with the reference merge. The median change is reported when at least 100 samples support it. Drift above 0.2 EV warns with automatic `ev` and fails with a manual `ev` unless `--allow-scene-drift` is given.
- **Given strength.** A manual `ev` is used as given. Measurement still runs when a calibration setting asks for it, and values differing by more than 0.2 EV warn.

### 7. Composition

The manifest format is in the [README](../README.md). Validity masks apply before clipping dilation, exposure and noise measurement, dark-model fitting, motion rejection and calibration. Each group is merged in camera space over its validity mask and multiplied by its per-colour attenuation correction `2^ev`; attenuated contributions are flagged `FROM_ND`. At each site a clean sample is preferred over an adjacent-clipping fallback, which is preferred over a lower bound, and between two lower bounds the higher (more informative) bound is kept. A later group replaces an equal-quality sample only inside its `contribution_mask`; everywhere else it fills only what earlier groups measured worse, so a site stays missing only when no group is valid there. A measurement more than 10% below an earlier lower bound is flagged `INCONSISTENT`: a filling group keeps the bound, a contribution mask takes the measurement. The mosaics of groups no later calibration needs are released after each merge. Missing eligible defects may be repaired, but sites with no eligible contribution stay uncovered.

**Offset matching.** Groups after the first are matched to the composite of earlier groups with an additive offset per colour, unless the manifest sets `"match": "none"`. The typical cause is the sun's veiling glare: an open bracket sees the scene plus a smooth veil that a blocked bracket does not, and pasting its samples in place of the blocker leaves a bright patch. The offset `D = L_group·2^ev − L_composite` is measured on sites where both are precise (clean, SNR of at least 10), outside the group's contribution mask, as the median in 32 × 32-site cells with at least 16 samples. Cells without a measurement take the membrane (harmonic) interpolation of the measured ones, and each inserted site subtracts the bilinear offset at its position. An offset larger than the site's own value is capped at that value, since the veil is part of what the group measured. Measured gain between real glare-affected brackets stays within a few percent of 1 inside a cell, so the additive model holds; it cannot recover structure that exists only inside the replaced area, such as a diffraction spike running along a blocker that points at the sun. The report gives the overlap samples, measured cells, corrected, capped and unresolved sites, and the median offset (camera RGB at output scale) and its median share of the group's signal.

### 8. Flat field

`--vignette dng` (default) applies every `GainMap` in OpcodeList2 to the merged sites it covers. The Bayer sample uses the pixel centre, `v = (row + 0.5) / image_height` and `u = (col + 0.5) / image_width`, matching [Adobe's GainMap reference implementation](https://android.googlesource.com/platform/external/dng_sdk/+/refs/heads/android14-prebuilt-test/source/dng_gain_map.cpp). Gains are interpolated bilinearly with edge clamping, and `L`, lower bounds and variance are scaled by the gain and its square. A DNG without OpcodeList2 gets no flat field, and the report says so. `--vignette FILE` takes an EXR of the full frame with three channels (a site takes its colour's channel) or one channel per site. `none` skips it.

### 9. Highlights and missing sites

These run after the flat field so colour comparisons use flat-fielded values.

- **Lower-bound sites** (`--highlights reconstruct`, the default) are estimated from the channels that did not clip, in white-balanced camera space, before the colour matrix mixes channels:
  - Quads with a lower-bound site form 8-connected regions. The unclipped quads within two quads of a region give its colour as ratios of channel sums, so bright rim quads outweigh dim, aliased edge quads. A region needs at least four rim quads.
  - Each clipped site is the mean of its 3 × 3 neighbourhood's unclipped other-colour sites, each multiplied by the region's ratio to the clipped colour. It never drops below its own bound, nor below a neighbouring bound carried through the ratio. Where all channels clip, the bounds set the level and the rim sets the colour.
  - Without a usable rim the colour is unknown and the estimate is neutral: the brightest carried value wins.
  - Raised sites keep `LOWER_BOUND`, gain `RECONSTRUCTED`, and inherit the donors' flags and sources. A value from a measurement carries that measurement's variance through the ratio; a value from a bound has zero variance.

  The estimate assumes the clipped area has the colour of its rim, which suits small sources such as the sun's disc and a lamp. A clipped area that differs in colour from everything around it cannot be recovered. `neutral`, the previous default, raises each lower-bound site to the quad's brightest white-balanced value: fully clipped areas come out neutral instead of magenta, but a site raised past its own bound (green to red's level, say) makes a partly clipped source yellow. `bound` leaves the raw bounds.
- **Missing sites** are filled with the variance-weighted mean of the nearest same-colour sites with data (R and B: the four sites two away horizontally and vertically, then the four diagonals. G: the four diagonal G neighbours) and flagged `FILLED`.

Repairs keep `MISSING` as a record of absent measurements and inherit quality flags and group source IDs from the donors actually used. Sites excluded by every contribution mask remain `MISSING | UNCOVERED` and provide no interpolation taps to their neighbours.

### 10. Demosaic, once

Demosaicing runs on the float mosaic, which can span more than 20 EV and is never converted to integers. `bilinear` (default) takes each missing colour as the mean of its nearest same-colour neighbours. It is linear, so it commutes with exposure scaling. `mhc` uses Malvar–He–Cutler 5 × 5 gradient-corrected filters, falling back to bilinear where the correction exceeds half the bilinear value or any site in the window is a lower bound, uncovered or an unfilled hole. Uncovered and unfilled taps are omitted. Image borders are mirrored.

`--edges ratio` (default) replaces the straight average on steep edges. Beside a source far brighter than its surroundings, such as the sun's disc, brightness can change tenfold per site. A mean of neighbours that differ 15–75× is dominated by the bright side, so the colours a site did not measure come out several times too high for the one it did. At a blue site the result is B/G far below the true ratio, which the AP0 matrix turns negative. There the missing colours are estimated in log space, which is exact for constant hue:

- Green at a red or blue site is taken along the direction with the smaller log gradient, Hamilton–Adams style. The site's own colour corrects the curvature (`½(g₋₁ + g₊₁) + ⅛(2c₀ − c₋₂ − c₊₂)`, exact to second order in log), and the result is limited to the range of the green neighbours.
- Red and blue are that green times the geometric mean of the neighbours' colour-to-green ratios, each neighbour using its own green estimate.

The log estimate blends in linearly in log contrast, from 2× to 4× between the largest and smallest taps a straight average would use. Every site in the 7 × 7 window must be measured or filled, positive and at least 8 σ above its noise, since logs of noisy values are biased; otherwise the site keeps the straight average. The kernel weights are rescaled to the new value, so variance keeps its relative size and flags and sources still come from the same taps. Both estimates scale with exposure. Pixels interpolated this way get `EDGE_RATIO`, and the report gives their count under `edges`. `--edges linear` keeps the straight average everywhere.

### 11. Colour

With `--space aces`, white-balanced camera RGB is multiplied by the rawtoaces matrix to give AP0. `rawfuse` configures an `rta::util::ImageConverter` from the reference DNG to obtain the white-balance multipliers and the 3 × 3 camera-to-AP0 matrix, then applies `out_i = Σ_j matrix[i][j] · in_j` itself, so variance propagates through the combined demosaic and colour weights. The multipliers are normalized so the smallest of R, G, B is 1. For a white balance other than the reference frame's (`median`, `custom`), rawtoaces is configured with `WBMethod::Custom`, because it derives the illuminant, and so the matrix, from the white balance. `--matrix` is passed through to rawtoaces' `matrix_method`. The default, `auto`, uses the spectral solver when the rawtoaces database has the camera's sensitivity curves and falls back to the metadata matrix otherwise, and `custom:FILE` holds nine numbers (camera RGB to AP0). With `--space camera`, camera RGB is written without white balance or matrix, and both are recorded in the report and EXR header. Forced `--matrix adobe` rejects DNGs with a non-RGB plane order, because LibRaw may supply an identity fallback matrix for them. rawtoaces' lens correction is not used.

### 12. Scale and write

- `--normalize geomean` multiplies `L` by `E_ref`, the geometric mean of the first group's effective exposures (equal to its EXIF geometric mean, see step 3), so values resemble a middle frame's. `none` writes `L` itself, a rate per second at ISO 100.
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
| `rawfuse:report` | the report's file name |
| `rawfuse:referenceFNumber` | aperture used to normalize effective exposure |
| `capDate` | capture date of the reference frame |
| `expTime`, `isoSpeed`, `aperture` | geometric-mean equivalent exposure time, ISO 100, and normalization aperture; omitted with `--normalize none` |

A non-DNG capture date is LibRaw's local wall-clock timestamp, formatted as local time.

### Auxiliary channels (`--aux`)

- `var.R`, `var.G`, `var.B`: variance of the output values.
- `frames`: number of frames that contributed (half float).
- `flags` (UINT bit field): bit 0 `LOWER_BOUND`, 1 `FROM_ND`, 2 `MOVING`, 3 `FILLED`, 4 `MISSING`, 5 `INCONSISTENT`, 6 `FLOORED`, 7 `ADJACENT_CLIP_FALLBACK`, 8 `UNCOVERED`, 9 `RECONSTRUCTED`, 10 `EDGE_RATIO`. `FROM_ND` marks any attenuated group's contribution. Repair, highlight and demosaic donors are included.
- `source` (UINT): 1-based group ID in manifest order (1 for positional frames), 0 for no source, 65535 for mixed sources. RGB source IDs include repair, highlight and demosaic dependencies.

`--write-cfa PATH` writes a one-channel EXR of the merged mosaic after the flat field and before repair and highlight handling, with `var` and `flags` channels and the CFA pattern as an attribute. Its source IDs describe selection before repair.

OpenImageIO normalizes integer channels to [0, 1] when asked for floats, so read `flags` and `source` in their native UINT type.

### JSON report

The report (`OUT.json` by default) records, among other things: the version, each input's path, group, role, exposure, ISO, f-number, EXIF and measured effective exposure, AsShotNeutral, ImageDescription, XMP, levels, optical-black spread, noise model and its source (`dng`, `custom`, `measured` or `generic`), and defect counts; the noise fits; the reference frame; the white balance and matrix used; the exposure reference; the flat field; per group its masks, ev, calibration, exposure steps, offset match, dark model fit and selected-site counts; merge statistics; and warnings.

## Errors and warnings

Exit code 2 with a one-line reason: frames of different size, CFA pattern or camera; a non-Bayer or linear raw; a frame without an exposure time; the same file given twice, including through a hard link; an existing output without `--force`; automatic `ev` without enough overlap; a manual `ev` with measured drift over 0.2 EV.

Progress and outcome lines go to stderr as each stage finishes: the frames read, each bracket's exposure range and how its steps were measured, the noise model, each group's attenuation, the sites it supplies and its offset match, the merge result and the files written. `-q` keeps only warnings and errors; `-v` also lists each frame as it is read.

Warnings go to stderr and the report: exposure steps that could not be measured or were not trusted, differing ISO between frames, a generic noise profile, flagged calibration blocks, scene drift, an unimplemented opcode, AsShotNeutral spread across frames, nonzero BaselineExposure and differing embedded GainMaps (the reference frame's maps are used). Warnings collected before a failure are printed with the error.

## Program structure

`rawfuse_pipeline` is a static library with a thin command-line front end in `src/main.cpp`.

| Module | Job |
|---|---|
| `options` | command-line parsing and validation |
| `metadata` | TIFF/DNG tag walker and opcode parsing |
| `raw_frame` | LibRaw reading: mosaic, geometry, levels, CFA, exposure |
| `dark_model` | dark-model CSV, fitting, per-site dark signal and extra variance |
| `bracket` | measured exposure ratios and noise model |
| `merge` | exposure merge and motion handling |
| `nd` | attenuation calibration and scene drift |
| `groups` | manifest, group masks and composition |
| `finish` | flat field, highlights, repair and demosaic |
| `colour` | white balance and matrix through rawtoaces |
| `writer` | EXR and JSON output |
| `pipeline` | orchestration |

## Known limits

- Variance covers the sensor noise model and interpolation weights. It excludes uncertainty in measured exposures, attenuation, colour transforms and flat-field fits. Repaired sites share donors, so their covariance is not fully modelled. Interpret lower bounds and repairs through the flags.
- Measured exposure ratios include any change of the light between neighbouring frames. A steady change during a bracket is absorbed into the exposures; a large one is rejected by the 0.25 EV limit.
- A measured noise model assumes a mostly static bracket at one ISO; at least a tenth of the shot-noise tiles must be free of motion and texture.
- A shared flat field is applied after fusion. Captures whose physical vignetting changes between frames need a deliberate per-frame calibration, not an average of maps.
- `blocks_flagged` are diagnostic and do not change the composition.
- Offset matching assumes earlier groups see the same scene; it only corrects smooth additive differences and reaches into a replaced area from its border.
- Automatic attenuation needs enough overlap between the reference and attenuated exposures. Without it, supply a measured `ev`.
- There is no registration or deghosting between independent sweeps. The caller must supply a fixed camera and a static scene.
- The CLI is the supported interface. `rawfuse_pipeline` is a build target, but there is no installed, versioned C++ package or stable ABI. Its low-level stages assume validated input, and `run()` sets process-wide OpenEXR/OIIO thread limits (restored on return).
- Peak memory is dominated by decoding all frames of a group at once, roughly 50 MB per full-resolution frame.
- Development and testing so far use macOS with synthetic Bayer DNGs. Linux and Windows CI and native CR3 coverage are not yet in place.

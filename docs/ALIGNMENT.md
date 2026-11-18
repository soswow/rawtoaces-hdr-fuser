# Alignment and splatting merge (design, in progress)

Status: stages 1 (motion estimation), 2 (motion-compensated merge on the
reference mosaic, `--align apply`) and 3 (full-colour merge, `--demosaic
splat`) are implemented on the `splat-merge` branch. Stages 4–5 are planned.
This note records the design so the stages can land one at a time.

## Goal

Fuse brackets and groups taken with small camera motion (a bumped tripod, a
hand-held bracket, a blocker that nudged the camera) without resampling the
Bayer mosaic and without a separate demosaic: each raw sample is placed at its
own motion-compensated sub-pixel position and accumulated into the output RGB
grid with a kernel, as in Wronski et al. 2019 (Handheld Multi-Frame
Super-Resolution) and Bayer drizzle in astrophotography. With no motion the
result must match today's merge plus demosaic within noise.

## Why not resample

Shifting a mosaic by a non-even amount moves red samples onto green positions.
Resampling each colour plane first (2 px spacing for R and B) blurs, correlates
neighbouring noise and breaks the per-site variance model. Even-pixel shifts
(HDRMerge) avoid that but leave up to 1 px of error. Splatting avoids both: a
sample is never interpolated, only weighted by its distance to the output pixel.

## Stages

1. **Motion estimation (done).** A motion field per frame, relative to the
   reference frame of the first group, in the JSON and status line
   (`--align report`).
2. **Motion-compensated merge on the reference mosaic (done).** `--align apply`
   reads every aligned frame through its field: the accumulation evaluated at
   the reference's own colour sites. Groups, calibration, highlights and the
   demosaic run unchanged on that mosaic.
3. **Full-colour accumulation (done).** `--demosaic splat` evaluates the same
   accumulation for all three colours at every pixel, so each frame's samples
   of every colour land where they fell and the separate demosaic goes. Repair,
   highlight reconstruction and offset matching work on per-colour sites with
   provenance.
4. **Per-pixel robustness.** Wronski's local-statistics test replaces quad
   motion rejection and decides where scene motion (clouds, foliage) is
   followed or rejected.
5. **Anisotropic kernels.** Structure-tensor kernels per output pixel (edge
   direction and strength from the merged luminance), which also replace the
   steep-edge ratio interpolation.

## Stage 1: motion estimation

Displacements are in site units: reference site `x` is seen at `x + d_k(x)` in
frame `k`.

**Images.** Each frame is reduced to a quad image: per 2 × 2 quad, the mean of
the four sites' black-subtracted fractions divided by the effective exposure,
defined where all four sites are usable (not clipped, defective or outside the
group's validity) and the mean has an SNR of at least 10. Matching uses the
logarithm, so a gain between frames (an ND group, a residual exposure error) is
a constant offset. Each quad carries an inverse-variance weight from the noise
model plus a 2 DN additive error (a frame's black level or dark signal), capped
at an SNR of 100. The additive term matters: in log it is an offset that grows
toward dark quads, and on a sky that brightens down the frame it reads as a
vertical shift. Without it, pairs of Katoomba frames reported about +0.3 sites
of row shift whichever frame was the moving one, and chains did not close.

**Matching.** Translation and log-gain are fitted by robust Gauss–Newton
(Cauchy weights on residuals normalized by their predicted noise). A global fit
on a pyramid (exhaustive search at the coarsest level, then every level,
subsampled to about 65,000 points) seeds a fit per tile of 64 × 64 quads
(128 × 128 sites) at the two finest levels. The finest level interpolates with a
Catmull–Rom cubic: bilinear interpolation pulled shifts toward whole pixels by
a few percent of their fractional part. A tile is measured when its fit
converges on at least a fifth of the tile with a translation uncertainty of at
most 0.2 sites. On synthetic frames with texture down to a 12-site wavelength
the error is about 0.01 sites; texture near the quad grid's limit (8 sites)
biases it by 0.02–0.03.

**Camera model.** Camera motion is rigid: a translation and a roll about the
frame centre. It is fitted to the measured tiles, translation first from their
median so the majority decides, then roll, with a fixed robust scale (tile
uncertainty plus 0.1 sites). A scale estimated from all tiles grows with the
share of moving ones; it let a cluster beside the sun tilt the model. Tiles
within max(0.5, 3σ) sites of the model are inliers. The frame's reported shift
is the model at the centre. The field keeps inlier tiles' own values and takes
the model elsewhere. Tiles that depart coherently with their neighbours are
labelled scene motion and reported but not followed (see open questions).

**Confidence.** An estimate is used only when the inliers number at least
max(16, 4% of tiles) and span at least 40% of the frame's height and width
(frames of fewer than 12 tiles: most tiles agree). A cluster of tiles (the sky
around the sun in an ND frame) usually measures the sun's or a cloud's motion.
Frames that fail, or that share too little with any aligned frame, take the
motion of the nearest exposure in their own bracket when one was measured, and
are otherwise reported unaligned.

**Which frame aligns to which.** Pairs inside a bracket see the same glare and
scene and match far better (residual scale 0.6–1.0 on Katoomba) than pairs
across groups (2–6, from the glare veil). Group 0 grows a tree from the
reference, each step pairing the unaligned and aligned frames that share the
most valid quads. Each later group is tied to the earlier ones by its single
best-overlapping pair (up to four candidates are tried), then grows inside
itself. Fields compose by addition, `d_k ≈ d_a + d_ka`.

**Output.** `--align off` (the default) skips the estimate, so output stays as
before; `report` estimates and reports it; `apply` also merges through it
(stage 2). The JSON `alignment` block lists per frame its parent, global shift
and uncertainty, roll, log-gain (`gain_ev`), residual scale, model RMS, the tile
grid (`dy`, `dx`, `measured`, `scene`) and whether the merge used it
(`applied`). Under `report`, a warning names how many frames are shifted by more
than 0.5 sites.

## Stage 2: motion-compensated merge

`sample()` is the one place frames are read for the merge, deghosting, ND
calibration and exposure ratios. A frame with a motion field is read at the
reference site's position in it, bilinear over the four nearest sites of the
same colour: red and blue sites form a square lattice of pitch 2, green sites a
quincunx, which is a square lattice of pitch 2 in u = row + col, v = row − col.
No colour is mixed into another, and a position on a site reads that site
alone, so zero motion reads the frame exactly as before. The sample's noise
terms are scaled by Σw², which the merge weights use. A position outside the
frame is missing; inside it, sites just beyond the edge mirror back, as in the
demosaic. Each site is scaled by its flat-field gain relative to the reference
site's: the frame saw the scene point through other sites, with their own
vignetting. (Applying the reference site's gain after the merge was off by the
vignetting gradient times the shift. On Katoomba the first version of stage 2
left 100 pixels with a negative channel at the image circle's rim; with this
and the interpolation weighting below, 2 remain, as without alignment. Which
of the two changes removed them was not separated.)

The merge also weighs interpolation error. Σw² alone favours the most averaged
sample, the best for noise but the most blurred. The detail around the pixel
is measured from the least noisy interpolated frame, as the spread of its sites
beyond ten times what their noise explains, and each interpolated frame's
error is (1 − Σw²) times it. The noise weight is multiplied by 1 / (1 + error /
σ²ₘ), σ²ₘ being the noise of the merge. Frames whose sites land alike get the
same factor, so a shared error (which no weighting removes) does not move
weight toward noisier frames; a frame on a site has no error and dominates
where detail is well above the noise. The gate keeps noise from reading as
detail: without it, frames with 1.5% noise came out 4% noisier than with plain
noise weights.

States combine conservatively. Any missing site (outside the frame or its
validity, a defect) makes the sample missing, so masks erode rather than blend;
any clipped site makes it a lower bound, the weighted sum of the clipped sites'
thresholds and the others' readings; a bloom site makes it an adjacent-clip
fallback.

After the estimate, exposure ratios are measured again through the fields
(starting from EXIF, so the first unaligned measurement does not compound). The
noise fit keeps its unaligned pairs: interpolation would shrink the pair
differences it measures. Frames whose motion is unconstrained or unknown are
read directly and named in the estimate's warnings.

This is interpolation of each colour plane, which stage 3 replaces: a moved
frame's samples are blurred by up to half a lattice pitch where they dominate
the merge, and the demosaic then interpolates the merged mosaic again. With
four 1–3-site shifts of a smooth texture, the merged green sites' RMS error
against the scene drops from 2.0% to 0.2% of the signal.

## Stage 3: full-colour merge

`--demosaic splat` runs the stage 2 accumulation for all three colours at every
pixel. A frame is read for colour c at the pixel's position in it, bilinear
over the four nearest sites of c; the frames are merged with the noise and
interpolation-error weights above. Without motion every frame interpolates
from the same sites, so the merge is the bilinear demosaic of the merged mosaic
(to 1e-5 in the tests, borders included). With motion each frame's samples land
where they fell: at a pixel where some frame measured a colour, that frame
dominates, so frames one site apart give every colour at every pixel.

After the merge the pipeline works on three sites per pixel (pixel-major):

- **Composition** ranks and replaces per colour sample. Masks are per pixel: a
  pixel's own site decides whether a group covers all three of its colours, and
  a pixel next to a validity edge has no reading of the colours whose sites lie
  beyond it, so a later group fills them. A mask that selects some Bayer
  colours and not others (a test checkerboard) therefore means something else
  here. Offset matching samples every colour at every pixel; a cell needs 16
  samples per site of its colour (64 for red and blue, 32 for green), so the
  field is anchored where the mosaic's is.
- **Flat field**: a full-colour site carries its colour's gains interpolated at
  the pixel with the demosaic's taps; the merge read sites relative to it.
- **Repair** fills a colour from the same colour of the axial, then diagonal,
  neighbours.
- **Highlights**: a clipped colour comes from the pixel's own unclipped colours
  times the clipped region's rim ratios, and never falls below its own bound or
  another colour's bound carried through them. `neutral` raises clipped colours
  to the pixel's brightest.
- **Render** applies the colour matrix to the pixel's three colours; there is no
  demosaic and no steep-edge estimate (`edges.method` is `none`).

Synthetic frames of an independent texture per colour, 192 × 192, RMS error
against the scene in camera RGB:

| frames | noise | stage 2 + bilinear demosaic | stage 3 |
|---|---|---|---|
| shifted (0,0), (0,1), (1,0), (1,1) | none | R 0.00029, G 0.00029, B 0.00022 | R 0.00003, G 0.00004, B 0.00003 |
| shifted 0.6–2.3 sites | none | R 0.00029, G 0.00029, B 0.00022 | R 0.00026, G 0.00025, B 0.00019 |
| shifted 0.6–2.3 sites | 1.5% | R 0.00069, G 0.00070, B 0.00054 | R 0.00055, G 0.00064, B 0.00045 |
| shifted 0.6–2.3 sites, flat scene | 1.5% | R 0.00031, G 0.00061, B 0.00028 | R 0.00044, G 0.00061, B 0.00039 |
| still | 1.5% | equal | equal |

On a flat scene stage 3 has the still frames' noise in red and blue; stage 2 is
quieter there only because it interpolates twice, which blurs detail as much.

## Katoomba side A (Canon R6, 8 mm fisheye, 21 frames)

- Open bracket against the blocked reference: **+2.4 to +3.0 sites in columns**,
  roll under 0.1 mrad, matching the bump noted when the captures were measured.
- Blocked frames 4–7 within 0.07 sites of frame 5; model RMS 0.08–0.13 sites.
- The three shortest frames of each bracket move 0.3–0.7 sites. Their only
  matchable content is sky, so the estimate follows the clouds as much as the
  camera: at 8 mm a site is about 0.6 mrad, and wind moves clouds several sites
  a second.
- ND frames: unconstrained. Only the sky beside the sun is matchable, in 11–18
  tiles. Their log-gain against the open bracket, −10.17 to −10.24 EV, agrees
  with the measured ND strength.
- Cost: 5.2 s of an 11 s run, mostly the tile fits and quad images.
- `--align apply` reads 13 frames through their fields (the ND frames
  directly). Where a dark silhouette crosses the open bracket's contribution
  boundary, its edge jumps 2.2 sites without alignment and 0.6 with it (within
  the silhouette's own raggedness). Open and ND attenuation change by under
  0.005 EV: both are medians over smooth sky. The run takes 18.8 s against
  8.2 s; the merge's four-site reads double its time.
- `--align apply --demosaic splat`: the same edge jumps 0.4 sites. Calibration
  is unchanged (it merges the mosaic), and the level step across the open
  bracket's mask boundary is within 0.006 EV in every colour, as for stage 2.
  Against stage 2, 2% of pixels differ by more than 7%, almost all at
  high-contrast detail. Without motion, against the bilinear demosaic with
  linear edges, 0.8% do: 0.2% at highlights and mask edges, and the rest in the
  open bracket's replacement area by the sun, where the offset field is
  extrapolated hundreds of sites from its anchors and differs by 1–4% between
  the two ways of sampling it. 102 pixels have a negative channel, in dark
  areas where red and blue noise is no longer smoothed by a demosaic. The run
  takes 22.7 s (merge 11.1 s) and peaks at 4.4 GB against 2.3 GB.

## Open questions

- **Scene motion.** Clouds and foliage are labelled per tile but not followed.
  Following them aligns their content but also follows bad fits (the dim rim
  of the image circle, validity-mask edges); stage 4's per-pixel robustness test
  is the place to decide, with the photometric residual as evidence.
- **Short frames.** Their camera motion should probably be shrunk toward the
  bracket's long frames unless the evidence is strong.
- **Unconstrained groups.** The ND group cannot be measured and is read
  unshifted, although it was shot between the open and blocked brackets. A
  manifest `shift` (TODO.md) would let the caller supply it.
- **Noise measurement** still pairs frames site by site, so handheld brackets
  inflate it; it needs aligned pairs with the Σw² factor.
- **Scene drift** sites in ND calibration read frames directly.
- **Steep edges under splat.** The full-colour merge has no log-ratio estimate,
  so the sun's rim fringes as with `--edges linear`. Stage 5's kernels are the
  place for it; until then a log-domain fallback per pixel could reuse the
  mosaic estimator's test.
- **Interpolation error** is measured from tap spread, which a linear gradient
  also produces, though bilinear interpolation is exact there. Curvature would
  be the better measure but needs more sites per frame.
- **Memory and time.** Three sites per pixel double the merge's memory; the
  per-frame reads (12 sites per pixel for a moved frame) dominate its time.

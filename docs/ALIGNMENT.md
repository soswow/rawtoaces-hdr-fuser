# Alignment and splatting merge (design, in progress)

Status: stage 1 (motion estimation, reported only) is implemented on the
`splat-merge` branch. Stages 2–4 are planned. This note records the design so
the stages can land one at a time.

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
   reference frame of the first group, reported in the JSON and status line.
   Nothing is warped yet.
2. **Sample cloud and splat accumulation for one group.** Each usable raw
   sample becomes (position in reference coordinates, colour, rate, variance,
   state). Output pixels accumulate kernel × inverse-variance weights per
   colour. Clipped samples contribute lower bounds; bloom and dark-model states
   keep their current meaning. Robustness (Wronski's local-statistics test)
   replaces quad motion rejection. Behind `--merge splat`; the mosaic path stays
   the default.
3. **Groups, calibration and highlights on RGB.** ND calibration, scene drift,
   offset matching and highlight reconstruction move from mosaic sites to RGB
   pixels with per-channel provenance. Exposure-ratio and noise measurement
   sample frame pairs through the motion fields.
4. **Anisotropic kernels.** Structure-tensor kernels per output pixel (edge
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

**Output.** `--align off` (the default) skips the estimate, so output stays as before; `report` estimates and reports
it. The JSON `alignment` block lists per frame its parent, global shift and
uncertainty, roll, log-gain (`gain_ev`), residual scale, model RMS and the tile
grid (`dy`, `dx`, `measured`, `scene`). `applied` is false: nothing is warped
yet, and a warning names how many frames are shifted by more than 0.5 sites.

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

## Open questions

- **Scene motion.** Clouds and foliage are labelled per tile but not followed.
  Following them aligns their content but also follows bad fits (the dim rim
  of the image circle, validity-mask edges); stage 2's per-pixel robustness test
  is the place to decide, with the photometric residual as evidence.
- **Short frames.** Their camera motion should probably be shrunk toward the
  bracket's long frames unless the evidence is strong.
- **Exposure and noise measurement** still pair frames site by site; stage 3
  samples them through the fields.

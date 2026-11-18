# TODO

## 1. Cross-group alignment

**Problem.** `rawfuse` assumes every frame sees the same scene point at the same
sensor site. Between `--groups` brackets (plain, sun blocked, ND, ...) the camera
can be nudged by a pixel or two. Measured on Katoomba: hand-blocked bracket is
about 2.3 px off the plain bracket in columns, ND showed no measurable shift.
Effects: doubled edges at group seams, biased auto attenuation calibration.

**Feature.** Estimate a small per-group geometric offset relative to a reference
group, report it, and optionally resample the group onto the reference grid before
calibration and merging. Off by default, so the "no resampling" guarantee holds
unless asked.

### Scope

- Between groups only. Frames inside one bracket are not aligned (tripod drift is
  ~0.2 px over a bracket). Revisit if the report shows otherwise.
- Reference group = the first group, or `calibration_reference` when set. It is
  never resampled.
- Transform: translation + roll, chosen over plain translation only when the fit
  improves. Anything richer (affine, homography, displacement field) is out of scope
  for v1.

### CLI and manifest

- `--align off|report|auto` (default `off`).
  - `report`: estimate and write to the JSON report, change nothing.
  - `auto`: estimate and apply when confident, otherwise warn and leave the group
    untouched.
- Per-group manifest key `shift: [dx, dy]` (sensor pixels, optional `roll` in
  degrees) overrides the estimate and is applied by the same resampler.
- Report per group: shift, roll, uncertainty, overlap sample count, applied yes/no
  and why not.

### Algorithm

1. Merge each group as today (`L`, variance, quality flags).
2. For each non-reference group, build a quad-averaged green image of `L` in the
   log domain, high-passed. This is gain-invariant, so ND strength is not needed
   first.
3. Mask: both groups' `valid_mask`, usable and unclipped in both, weighted by
   gradient and variance. Static foreground drives the fit, a robust loss rejects
   moving sky.
4. Coarse-to-fine masked sub-pixel registration. Accept the result only if the
   uncertainty is below ~0.2 px and there is enough overlap texture.
5. Apply (auto/manual): demosaic the group's merged `L`, warp with Lanczos3, sample
   at the reference CFA positions, hand it back as that group's mosaic. Erode clean
   masks and dilate clipped/adjacent-clip masks so lower bounds never blend into
   measurements. Variance scales by the squared kernel weights.
6. Run calibration, dark model, contribution masks and the rest unchanged.

### Acceptance

- Synthetic: known sub-pixel shift and roll recovered within 0.1 px, zero shift
  recovered as zero, sparse ND-like overlap reports "uncertain" instead of a wrong
  number.
- Katoomba: blocked vs plain recovers ~2.3 px in columns. Auto ND EV with and
  without alignment compared and recorded.
- `--align off` output is bit-identical to today.
- Edge rise-time of the warped patch vs reference measured and recorded.

### Build order

1. `report` mode (estimator only). Run on Katoomba before going further.
2. Manifest `shift` + resampler.
3. `auto` wiring, thresholds, docs ([ALGORITHM.md](docs/ALGORITHM.md) Known limits).

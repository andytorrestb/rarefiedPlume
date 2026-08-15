# Statistical health of the Cai Kn = 100 case

A 3 × 3 sweep over the two knobs that set statistical quality — **parcels per
cell** and **averaging duration** — at fixed `Kn = 100`, judged primarily
through imagery.

[`cases/cai2012`](../cai2012/) measures how well `dsmcFoam` reproduces Cai's
collisionless solution. It never measures how much *statistical noise* its own
answer carries.
[`plumetools/cai2012/checks.py`](../../plumetools/cai2012/checks.py)
`::_check_occupancy` prints an **estimated** exit-cell occupancy and says
outright that only a post-run audit of the sampled `dsmcRhoN` field can say what
it actually was. This family is that audit.

The occasion for it is that the rendered time series for
`cases/cai2012/Cases/Kn100` is visibly speckled. That is shot noise, not a bug —
those frames are instantaneous `rhoN`, one timestep's parcels, and `Kn = 100` is
the most rarefied case in the matrix. The question this study answers is whether
it is *only* that.

## Quick start

```bash
cd cases/cai2012-health

./generate_cases.py --dry-run    # the matrix, and what it will cost
./generate_cases.py              # nine cases, Cases/<weight>/<sampling>
./AllmeshCases                   # blockMesh + topoSet + createPatch, nine times
./AllrunCases                    # solve, cheapest first, timing each
./AllpostCases                   # validate, audit, render, dissect
```

Get the machinery working before committing the compute:

```bash
./AllrunCases --only Cases/ppc005/s0p5 --only Cases/ppc005/s1p5
```

is the two cheapest cases — about 40 minutes — and exercises the whole
generate → mesh → run → sample → render → dissect chain, including the
determinism check on its first real pair.

## The sweep

| | | |
|---|---|---|
| `resolution.target_particles_per_cell` | **5 / 20 / 40** | rows, `Cases/ppc005/` … |
| `dsmc.sampling_domain_transits` | **0.5 / 1.5 / 4.5** | columns, `.../s0p5/` … |

5 is `checks.min_particles_per_cell`, the documented floor. 20 and 1.5 are what
`cases/cai2012` runs every one of its three cases at, so **`ppc020/s1p5` is that
family's `Kn100` case re-run** and is the control for the whole matrix.

**Everything else is held fixed** — the same mesh (1 310 720 cells, 10 mm core),
the same Knudsen number, the same 1.373 × 10⁻⁶ s time step, the same 2 domain
transits of transient, the same write schedule — so every case is comparable
image to image. That is enforced rather than intended:

* `study.yaml` has **no override mechanism**, unlike `cases/cai2012`'s. The two
  axes are the entire per-case input; `generate_cases.py` writes exactly two
  keys into each `case.yaml` and rejects an `overrides:` block by name. Anything
  else worth changing goes in `baseCase/case.yaml`, where it applies to all nine
  at once.
* `manifest.yaml` records the mesh size, time step, transient and write interval
  **per case**, and says `NOT SHARED`, loudly, if any of them ever differs. Two
  cases on two meshes would still generate, still run, and still draw nine
  plausible pictures.
* `AllmeshCases` meshes all nine independently and then checks the cell counts
  agree. They *are* all the same mesh — the weight axis changes only
  `dsmc.n_equivalent_particles` — and meshing once and copying would be cheaper,
  but a study whose comparability rests on the meshes being identical should
  demonstrate that rather than arrange it.

### Why this one uses CaseFoam and `cases/cai2012` does not

[`cases/cai2012/README.md`](../cai2012/README.md) records that family's decision
not to use CaseFoam, and the reason holds: it is a **flat list** of three
Knudsen numbers, and requiring a dependency that would do a `copytree` and
nothing else would make `./generate_cases.py` fail on a machine where the core
install works.

This family is a 3 × 3 **hierarchy** — weight over sampling duration — and
hierarchies are what CaseFoam is for. It earns the dependency exactly the way
`cases/markelov1999`'s `<gap>/<pressure>` does, and `Cases/ppc020/s1p5/` is
CaseFoam's own `tree` layout. The two decisions are the same rule applied to two
different shapes, not a contradiction.

[`plumetools.markelov1999.study.clone_cases`](../../plumetools/markelov1999/study.py)
could not be reused directly: its first level is a single `gap_dir` and it
`rmtree`s the write directory on entry, so three calls to build three weight
rows would each destroy the last. Its `require_casefoam` and its
`mkCases(..., hierarchy="tree")` call are what get reused, in
[`plumetools/cai2012/health.py`](../../plumetools/cai2012/health.py).

### Many frames per run, and why the interval is what it is

`dsmc.write_interval_s` is pinned at **3.536 × 10⁻⁴ s — one eighth of a domain
transit, 258 timesteps**. `cases/cai2012` leaves it `null`, which derives about
two writes inside the sampling window: enough to post-process, not enough to
watch anything converge.

| sampling | frames with a `*Mean` field |
|---|---|
| 0.5 transits | 5 |
| 1.5 transits | 13 |
| 4.5 transits | 37 |

The 1:3:9 ratio of the axis itself, so **frame *k* is the same elapsed averaging
time in every case of a row**. Finer (1/16 transit) puts the matrix past 300 GB
and makes the imagery step outlast the solve; coarser (1/4) leaves the shortest
case with two frames, and two points are not a curve.

It is fixed in **seconds, not frames**. Every case therefore writes at the same
solver times, which is what the overlap check and the common-time contact sheet
line up on. A per-case frame count would put the frames at different times in
different cases and neither would have anything to compare.

Fifteen further frames are written during the transient, before `fieldAverage`
starts. They carry no `*Mean` field and are excluded from the imagery by
`viz.yaml`'s sampling window — see below.

## What it costs

`./generate_cases.py --dry-run` prints this before writing anything, because
run time goes as parcels × total transits and the top weight row is most of the
bill however the matrix is shaped:

```
  case                ppc  transits    parcels    steps   hours  writes      GB
  ppc005/s0p5           5       0.5   2.89e+05    5,160    0.28      20     4.8
  ppc005/s1p5           5       1.5   2.89e+05    7,224    0.39      28     6.7
  ppc005/s4p5           5       4.5   2.89e+05   13,416    0.73      52    12.4
  ppc020/s0p5          20       0.5   1.16e+06    5,160    0.69      20     4.8
  ppc020/s1p5          20       1.5   1.16e+06    7,224    0.97      28     6.7
  ppc020/s4p5          20       4.5   1.16e+06   13,416    1.80      52    12.4
  ppc040/s0p5          40       0.5   2.31e+06    5,160    1.25      20     4.8
  ppc040/s1p5          40       1.5   2.31e+06    7,224    1.74      28     6.7
  ppc040/s4p5          40       4.5   2.31e+06   13,416    3.24      52    12.4

  TOTAL      11.1 h of solver, 71 GB of time directories
```

The cost model has two terms — per parcel and per cell — because the weight axis
multiplies the first by eight and leaves the second alone. A single-term model
fitted on the cheap row understates the expensive one, which is the direction
that costs half a day to discover. `AllrunCases` records each case's wall clock
in `results/cost-model.yaml` and the estimate is refitted from it as the matrix
runs; until then it is labelled `prior` wherever it is printed.

**Those hours are not all solver, and the prior assumed they were.** Measured on
`ppc005/s0p5`: `dsmcFoam` took 646 s of the run, and `reconstructPar` took about
40 s per written time directory — 13 minutes for that case's twenty writes, and
proportionally more for the 52-write cases. Reconstruction is single-threaded
and scales with the *number of frames*, which is the one thing this family
deliberately maximises, so it adds roughly 3.4 h across the matrix.

It cancels out: the solver runs about 1.6× faster than the prior expected and
reconstruction makes up the difference, so the ~11 h total stands. It is
recorded because the *composition* is wrong in the prior, and because
`AllrunCases` times the whole of `Allrun` — solve, reconstruct and all — the
fitted model absorbs it automatically. A model fitted on solver time alone would
under-predict every long case in the matrix.

**Rendering is a third cost again.** 165 sampled frames across the nine cases ×
4 fields ≈ 660 images at ~20 s each ≈ **3.7 h**, in one `vifpara` invocation.
That is the price of the write interval this family chose, and it is why
`docs/viz-slices.md`'s warning about the imagery step outlasting the solve is
taken seriously here rather than quoted.

### What was trimmed, and what that costs

**The weight axis was cut from 80 to 40.** At 80 that row alone is ~14 h and the
matrix is 18.8 h. At 40 it is 11.1 h.

What is lost is the top half of a decade on the parcels axis. The remaining span
is 8× (5 → 40), which still brackets the floor and the baseline, and the
sampling axis spans 9× on its own — so the combined statistical budget still
covers a factor of 72 and the `1/√N` test survives. It is weaker than it would
have been, not absent.

**Disk, not CPU, is the tighter constraint.** A reconstructed time directory of
this mesh is ~255 MB of fields plus a parcel cloud that scales with the weight
(~200 MB at `ppc020`). `AllrunCases` prunes the reconstructed clouds from every
time but the newest as each case finishes — nothing here reads parcel positions,
and `processor*/`, which is what a resume actually reads, is untouched. Without
that the matrix is ~150 GB instead of ~71.

## The imagery, and dissecting it

This is the heart of the study, not a postscript.
[`viz.yaml`](viz.yaml) is the authoritative list of what is drawn, and it
differs from both existing study specs in three ways.

### 1. `prefer_mean: true`, which is the opposite of everywhere else

[`docs/viz-slices.md`](../../docs/viz-slices.md) argues that a *series* should
use the instantaneous field, because `rhoNMean` is a running average and "a
series of it shows the average converging, not the flow evolving — by
construction it barely moves." Every word of that is true.

**Here the convergence is the measurement.** This family is not watching a plume
fill a vacuum; it is watching a running mean settle, and how fast it settles is
precisely what the two axes are supposed to change. The property that makes an
averaged series useless for that study is the property this one exists to look
at. Both settings are kept and both documents say why.

### 2. A sampling window, without which the first setting is a lie

`fieldAverage` writes no `*Mean` field before its `timeStart`, and these runs
write 15 frames before then. With `prefer_mean` on and no window, those frames
have nothing to reach for and fall back to the **instantaneous** field: the
series starts as shot noise and becomes a running average part way through, with
nothing in the pictures to say so. Both halves look like a plume.

`sampling.time_min: 5.658e-3` — `dsmc.average_start_s` — keeps the series to the
frames where the average exists. The mechanism is new
(`plumetools/viz/resolve.py`); the trap it closes is the one
`docs/viz-slices.md` already warned about.

### 3. Pinned ranges on every field

Colour ranges are auto-pinned across the **frames of one series** but **not
across cases** — each case gets its own union. Nine cases would therefore be
drawn on nine different colour scales, and a contact sheet of nine
differently-scaled images compares nothing while looking entirely fine: the
cells differ visibly, and the differences are the colour maps.

Every field carries an explicit `range:`, measured off the baseline case. Not
guessed: the first draft of this file pinned `dsmcRhoN` three decades high and
rendered a black rectangle, which is how the next section was discovered.

### `dsmcRhoN` is a parcel COUNT, and the repository had it wrong

`dsmcRhoN` is drawn most prominently here because it **is** axis 1 of the sweep,
mapped over the domain. It is the one field whose absolute numbers matter rather
than its shape — the shape is `rhoN`'s, since the two differ by the cell volume.

Measured on `cases/cai2012/Cases/Kn100`, over 1 018 211 occupied cells:

```
rhoN * V / dsmcRhoN = 3.3100e+09 = nParticle    to 1.5e-9
dsmcRhoN peaks at 20.47                          configured: 20 per exit cell
```

`DSMCCloud::calculateFields` adds 1 per parcel and never divides by the cell
volume, where `rhoN` accumulates `nParticle/V`. **So `dsmcRhoN` is the parcel
count in the cell — dimensionless — and the colour bar is the parcels per cell,
read directly against `target_particles_per_cell` and
`min_particles_per_cell`.**

`plumetools/viz/catalog.py`, `plumetools/viz/slices.yaml` and
`cases/cai2012/viz.yaml` all described it as a number density in `m⁻³` whose
product with the cell volume was the occupancy. That reading is out by `1/V` — a
factor of 10⁶ on this mesh — and on a logarithmic scale it drew a perfectly
plausible picture with `[m^-3]` on the bar. All three are corrected.

### The dissection

[`plumetools/viz/dissect.py`](../../plumetools/viz/dissect.py) sits beside
`video.py` under the same contract — ffmpeg only, no ParaView, importable,
runnable under plain `python`. `./AllpostCases` runs it through
[`analyse.py`](analyse.py), which writes into `results/`:

| | |
|---|---|
| `contact-final.png` | the matrix at each case's own final frame — **the headline** |
| `contact-common.png` | the matrix at one time every case reached — the control |
| `convergence.yaml` | frame-to-frame PSNR per case, both against the previous frame and against the final one |
| `sweep-table.csv` | one row per case: occupancy, budget, error |
| `sweep-error.png` | the physical error across both axes, against `1/√N` |
| `overlap.yaml` | long runs against short ones, where they overlap |

The two sheets are easy to confuse and answer different questions. **At a
literal common time the three cases of a weight row are the same run truncated,
so their cells should be pixel-identical** — that sheet is a determinism check
you can see, not a comparison. The headline sheet uses each case's own final
frame, which is the answer that case actually delivers.

### An image metric is not a physical error

PSNR between two log-scaled, colour-mapped, palette-quantised PNGs is **a
perceptual proxy for "has this stopped changing"**. The scale is logarithmic and
clipped to the pinned range; the 8-bit palette hides any change below one colour
step and reports `inf` when two frames are identical (recorded as `saturated`,
never as a number); and the vacuum counts as many pixels as the plume.

It is reported **beside** the centreline error from
`plumetools.cai2012.post.centerline_metrics`, never on its own. The caveat is
attached to the written data as well as the report, and `convergence_report`
prints a warning in place of the physical column when `metrics.yaml` is missing.
Where the two disagree, that disagreement is itself a finding: an image that has
stopped changing while the centreline error is still falling means the change
has moved below a colour step and the picture can no longer show it.

## The overlap, and why it is checked

A 4.5-transit run **contains** the 1.5-transit answer as an earlier frame. Six
of the nine runs are partly redundant by construction, and within a weight row
the three cases differ *only* in `endTime`: same mesh, same particle weight,
same time step, same decomposition, same seeding.

So the short run is a **prefix** of the long one, and `fieldAverage` — which
starts at the same `timeStart` in both and stores its accumulators per time
directory — must have accumulated exactly the same numbers by any time both
reached. `results/overlap.yaml` asserts it, field by field, byte for byte,
falling back to a numeric comparison only when the hashes differ so a
disagreement is quantified rather than merely flagged.

An agreement is not evidence that the physics is right. It is evidence that the
pipeline is **reproducible**, which is the precondition for everything else
here: this study attributes its differences to statistics, and if the redundancy
does not hold, some of them are something else. A failure would be a more
important finding than anything else in the study, and `analyse.py` exits
non-zero on one.

## Layout

```
cases/cai2012-health/
  study.yaml            the two axes, and nothing else
  generate_cases.py     clone via CaseFoam, write two keys per case, cost it
  analyse.py            overlap, sweep table, contact sheets, curves
  viz.yaml              what gets drawn -- see above
  baseCase/
    case.yaml           THE authoritative source of physical inputs
    Allaudit.py         occupancy + per-frame centreline error   (new here)
    Allmesh Allrun Allpost Allclean postProcess.py runInflow.py
                        byte-identical copies of cases/cai2012/baseCase
  Cases/                generated; gitignored
    ppc005/ ppc020/ ppc040/
      s0p5/ s1p5/ s4p5/
        results/        centerline.csv, metrics.yaml, occupancy.yaml,
                        convergence.yaml, viz/
  manifest.yaml         generated; the auditable case -> inputs map
  results/              generated, except cost-model.yaml
```

The runner scripts are **copies**, not a fork: a health case *is* a `cai2012`
case, running the same `Allmesh`, `Allrun`, checks and post-processing with two
values changed. `tests/unit/test_cai2012_health.py` asserts the byte-identity,
so a fix applied to one and not the other is caught rather than discovered.

`Allaudit.py` is the one script that is genuinely new, because the two things it
measures — occupancy per region, and the centreline error at *every* sampled
frame rather than only the last — are what `cases/cai2012` does not do.

## Results

**Not yet run.** The machinery is validated end to end on the two cheapest
cases; the full matrix is 11.1 h of solver and ~2 h of rendering.

This section will carry the sweep table, the contact sheets and the answers to
the five questions the study exists to settle:

1. How many parcels per cell does the sampled solution actually have, per
   region, as opposed to the estimate `checks.py` prints?
2. How far below `target_particles_per_cell = 20` and the floor of 5 does the
   real occupancy fall once you are away from the exit cell?
3. Which axis buys more — more parcels, or more averaging? They should trade
   off, since the mean of *N* samples over *T* transits is the same statistical
   budget however it is split. Does the data agree, and where does it stop?
4. Where does the centreline error stop improving — where does statistical noise
   stop being the limiting factor and something systematic take over?
5. What does this say about the three cases in `cases/cai2012`? Are their
   published numbers statistically sound, and is the `Kn100` imagery noisy
   because the case is under-resolved or because instantaneous fields simply
   look like that?

Two figures are already measured, from `cases/cai2012/Cases/Kn100` itself:

| | |
|---|---|
| exit-cell occupancy, **measured** | 18.91 parcels/cell (median over 316 exit cells) |
| exit-cell occupancy, **estimated** by `checks.py` | 20.0 |
| median over every occupied cell | **0.373** |
| cells with no parcel at all | 22.3% of the mesh |

The estimate is good where it is made. The point of the study is what happens
everywhere else.

## Related documents

* [`cases/cai2012/README.md`](../cai2012/README.md) — the family this measures
* [`docs/cai2012-case.md`](../../docs/cai2012-case.md) — paper values,
  assumptions, and every difference from Cai's setup
* [`docs/viz-slices.md`](../../docs/viz-slices.md) — the imagery pipeline, and
  the opposite `prefer_mean` argument for the opposite purpose

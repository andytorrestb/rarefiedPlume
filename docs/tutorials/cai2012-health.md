# cai2012-health: how much statistical noise is in the answer

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | A 3 × 3 sweep over parcels per cell and averaging duration at fixed Kn = 100, judged through imagery and centreline error |
    | **Model** | A `cai2012` case with two values changed per member; the same runner scripts, byte for byte |
    | **Evidence** | Unit-tested study, audit and cost model; a recorded determinism (overlap) result on the cheapest pair. The sweep's conclusions are **in progress** and not established here |
    | **Runs on** | Python + **CaseFoam** for generation (`--dry-run` needs no CaseFoam); OpenFOAM v2512 + `plumeFieldInflow` for meshing and solving; ffmpeg for the contact sheets; `pvpython` for imagery |
    | **Start in** | `cases/cai2012-health/` |

## What the study is for

[`cai2012`](cai2012.md) measures how well `dsmcFoam` reproduces Cai's
collisionless solution; it never measures how much statistical noise its own
answer carries. This family is the post-run audit `cai2012.checks` says is
needed: what occupancy the sampled solution actually had, region by region,
and whether more parcels or more averaging buys a better answer.

The design — why the axes, why CaseFoam here and not in `cai2012`, why the write
interval, why `prefer_mean: true`, and what has been measured — is in
[cai2012-health, explained](../cai2012-health.md).

## The sweep

| Axis | Key in each generated `case.yaml` | Values | Directory |
|---|---|---|---|
| parcels per exit cell | `resolution.target_particles_per_cell` | 5 / 20 / 40 | `Cases/ppc005/`, `ppc020/`, `ppc040/` |
| averaging duration | `dsmc.sampling_domain_transits` | 0.5 / 1.5 / 4.5 transits | `…/s0p5/`, `s1p5/`, `s4p5/` |

`ppc020/s1p5` is `cai2012`'s Kn100 case re-run at a finer write interval, and is
the control. Everything else is fixed in `baseCase/case.yaml` for all nine:
the mesh (1 310 720 cells, 10 mm core), Kn = 100, the time step, two transits of
transient, and `dsmc.write_interval_s` = 3.536 × 10⁻⁴ s. `study.yaml` has **no
override mechanism**; `generate_cases.py` writes exactly the two axis keys and
rejects an `overrides:` block. `manifest.yaml` records the shared quantities per
case and prints `NOT SHARED` if any ever differs.

## Prerequisites

| Step | Needs |
|---|---|
| `generate_cases.py --dry-run` | Python with `plumetools`; reads `results/cost-model.yaml` |
| `generate_cases.py` | **CaseFoam** (`pip install -e ".[cases]"`, or `./Allsetup`); it builds the `Cases/<weight>/<sampling>` tree with CaseFoam's `tree` layout |
| mesh / run / post | as for [cai2012](cai2012.md#prerequisites) |
| `AllpostCases` sheets | ffmpeg (`--no-viz` skips everything that needs it) |
| `AllpostCases` imagery | `pvpython` (`PVPYTHON_PATH` to name one) |

CaseFoam is not installed in this checkout's environments, so full generation
was not exercised for this documentation; `--dry-run` was.

## Inputs

| File | Role |
|---|---|
| `study.yaml` | `weights` (name, `target_particles_per_cell`, `enabled`) and `samplings` (name, `sampling_domain_transits`, `enabled`); nothing else per case |
| `baseCase/case.yaml` | every other physical and numerical input, shared by all nine — a `cai2012` `case.yaml` with `exit.knudsen: 100` |
| `viz.yaml` | what is drawn: `prefer_mean: true`, a sampling window starting at `dsmc.average_start_s`, a pinned colour range per field |
| `results/cost-model.yaml` | **tracked** measured wall clock of every case run so far; turns `--dry-run`'s estimate from a prior into a fit |

Schemas: [configuration reference](../reference/configuration.md#cai2012-health).

## Cheap check (no OpenFOAM, no CaseFoam)

```bash
cd cases/cai2012-health
./generate_cases.py --dry-run          # the matrix, the predicted cost, the destinations
./generate_cases.py --dry-run --weight ppc005 --sampling s0p5 s1p5
```

It derives the template once, costs each selected case from the cost model, and
prints `--dry-run: nothing written.` The cost is labelled `prior` until enough
measurements exist to fit both of its terms (per parcel and per cell); a fit on
one weight row cannot separate them and falls back to the prior's ratio.

## Production run

```bash
cd cases/cai2012-health
./generate_cases.py                    # nine cases via CaseFoam, manifest.yaml
./AllmeshCases                         # mesh all nine independently; cell counts must agree
./AllrunCases                          # solve, cheapest first, timing each into results/cost-model.yaml
./AllpostCases                         # Allpost + Allaudit per case, imagery, then analyse.py
```

Exercise the whole chain on the two cheapest cases before committing the rest:

```bash
./AllrunCases --only Cases/ppc005/s0p5 --only Cases/ppc005/s1p5
```

`AllrunCases` options: `--only PATH` (repeatable), `--serial`, `-np N`,
`--no-solve`, `--keep-going`, and `--keep-clouds`. By default it **prunes the
reconstructed parcel clouds** (`<time>/lagrangian/`) from every time but the
newest as each case finishes — nothing downstream reads parcel positions, and
`processor*/`, which a resume reads, is untouched. Each case **continues** from
its newest complete checkpoint, exactly as in `cai2012`; finished cases report
and are skipped in seconds and are not re-timed.

`AllpostCases` options: `--no-viz`, `--viz-only`, `--analyse-only` (re-cut the
sheets and tables from what exists), `--no-audit` (skip the per-frame audit),
`--no-analyse`, `--stride N` (audit every Nth frame), `--viz-spec F`,
`--keep-going` (cases that have not run report and are skipped).

!!! danger "Destructive options"
    `./AllcleanCases --purge` deletes the whole `Cases/` tree and
    `manifest.yaml`; `./generate_cases.py --overwrite` replaces cases that
    contain results. Preview cleaning with `./AllcleanCases --dry-run`.

!!! warning "Template fixes do not reach generated cases"
    Cloning is a copy. A fix to a `baseCase` script reaches existing cases only
    if you copy it across (runner scripts are not physical inputs), for example
    `for c in Cases/*/*/; do cp -p baseCase/Allaudit.py "$c"; done`.
    `--overwrite` would destroy the results.

## Outputs

Per case, in `Cases/<weight>/<sampling>/results/`: `cai2012`'s
`centerline.csv`, `density_plane.csv`, `metrics.yaml` and figures, plus from
`Allaudit.py`:

| File | Contents |
|---|---|
| `occupancy.yaml` | measured parcels per cell by region (exit cells, inside Cai's lowest contour, the centreline tube, the whole domain), against the floor and target |
| `convergence.yaml` | the centreline error at every sampled frame, not only the last |
| `viz/` | the imagery `viz.yaml` asks for |

Study level, in `results/` (from `analyse.py`):

| File | Contents |
|---|---|
| `contact-final.png` | the matrix at each case's own final frame — the headline |
| `contact-common.png` | the matrix at one time every case reached — a determinism check: within a weight row the cells should be pixel-identical |
| `convergence.yaml` | frame-to-frame PSNR per case (against the previous and the final frame) |
| `sweep-table.csv` | one row per case: occupancy, statistical budget, error |
| `sweep-error.png` | the physical error across both axes against 1/√N |
| `overlap.yaml` | long runs against short ones where they overlap, field by field |

`dsmcRhoN` is a **parcel count per cell** (dimensionless), not a number
density: `DSMCCloud::calculateFields` adds one per parcel and never divides by
the volume. Its colour bar reads directly against the 5/20/40 targets.

**Judging the output.** `analyse.py` exits non-zero if a short run is not a
prefix of its long one in `overlap.yaml`: that would mean the pipeline is not
reproducible and the sweep's differences are not purely statistical. PSNR is a
perceptual proxy for "has this stopped changing", not a physical error; it is
reported beside the centreline error and never alone. An image that has stopped
changing while the centreline error still falls means the change has dropped
below one colour step.

## Resources

The cost model is in `results/cost-model.yaml` (tracked, measured on the
owner's runs) and `--dry-run` prints the current estimate. Recorded
observations, with their provenance, are in
[What it costs](../cai2012-health.md#what-it-costs): `reconstructPar` is
single-threaded and scales with the number of frames, so total time exceeds
solver time; disk, not CPU, is the tighter limit (pruning the clouds roughly
halves it); rendering is a separate, multi-hour cost.

## Variants and safe customisation

| Want | Do |
|---|---|
| Fewer cases | `enabled: false` on a weight or sampling row, or `--weight` / `--sampling` filters |
| A different fixed setting for all nine | edit `baseCase/case.yaml`; regenerate |
| Different colour ranges | `viz.yaml` `range:` per field — measure them off the control case; never per case |

## Figure placeholder

!!! figure "Figure placeholder — the 3 × 3 contact sheet at each case's final frame"
    **Subject:** `results/contact-final.png` for `dsmcRhoN` (parcels per cell).
    **Source run:** all nine cases solved; `./AllpostCases` (needs ffmpeg and `pvpython`).
    **Plotting command:** `./AllpostCases --analyse-only` re-cuts it from existing imagery (verified in source).
    **Destination:** `docs/assets/figures/cai2012-health-contact-final.png`
    **Suggested caption:** "Parcels per cell across the parcel × averaging matrix at Kn = 100, each case at its own final frame (pinned colour range)."

## Related

- [cai2012-health, explained](../cai2012-health.md) — design, costs, imagery, recorded results
- [cai2012 tutorial](cai2012.md) — the family this audits
- [Slice imagery](../viz-slices.md) — the imagery pipeline and the opposite `prefer_mean` argument
- API: [`cai2012.health`](../reference/api/cai2012/health.md), [`cai2012.audit`](../reference/api/cai2012/audit.md), [`viz.dissect`](../reference/api/viz/dissect.md)

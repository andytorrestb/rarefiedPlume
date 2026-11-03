# Kn extension: Cai's plume across seven decades of Knudsen number

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | Cai's circular argon plume (D = 0.2 m, S₀ = 2, T₀ = 300 K) at Kn = 100, 10, 1, 0.1, 0.01, 0.001, 0.0001, with `cai2012-health`'s highest-resolution settings and Cai's 10⁴ t₀ startup; resolution diagnostics on every case |
    | **Model** | `cai2012` direct-exit cases, generated and assessed by `plumetools.cai2012ext` |
    | **Evidence** | Generation, policy, feasibility, comparison, HPC wrappers and diagnostics are unit-tested. Preparation (dictionaries, idempotence, guards) was checked on 2026-10-26. **No production DSMC run of the seven-decade suite is recorded; nothing is validated** |
    | **Feasibility** | Kn ≥ 0.1: production *candidates* before running. Kn = 0.01, 0.001, 0.0001: **unsuitable for production** as configured |
    | **Runs on** | Python for generation/validation/smoke; OpenFOAM v2512 + `plumeFieldInflow` for mesh and run; PBS for the committed job scripts |
    | **Start in** | `cases/cai2012-kn-extension/` |

This page covers the **direct-exit** cases. The two other workflows of the
family have their own tutorials: the [interface track](cai2012-interfaces.md)
(removing the dense core at Kn ≤ 0.001) and [field scaling](cai2012-scaling.md)
(empirical whole-field maps). The method, inherited settings, metrics and open
problems are in [The Cai 2012 Knudsen extension](../cai2012-kn-extension.md).

## Which suite is active

`study.yaml` declares `meta.knudsen_suite: full_decades`, which selects the
**seven-decade** suite and writes it to `Cases_full_suite/` with
`manifest-full-suite.yaml`. A study without that declaration keeps the
historical **four-Kn** policy (Kn = 1e-4, 1e-3, 1, 10, in `Cases/` and
`manifest.yaml`), which refuses Cai's own 100, 0.1 and 0.01 at every step. The
profile travels from the study to the manifest and into each generated case's
metadata, and generation, the runners, the PBS guard and post-processing all
read it (`plumetools.cai2012ext.policy`). Older pages that describe "the four
Knudsen numbers" describe the historical policy.

The historical trees (`Cases/`, `manifest.yaml`) and any local graded archive
such as `Cases_run1_graded/` are left as they are; nothing in the active suite
imports their results.

## Inputs

| File | Role |
|---|---|
| `study.yaml` | the suite: `cases_dir`, `manifest_name`, the seven `kn_cases` (with per-Kn mesh and parcel-target `overrides`), one particle level `np1x`, `meta.knudsen_suite` |
| `baseCase/case.yaml` | the template, a `cai2012` `case.yaml` that must match `cases/cai2012-health/baseCase/case.yaml` **key by key** except the overrides declared, with reasons, in `policy.DELIBERATE_OVERRIDES` — checked on every generation |
| `baseCase/diagnostics.yaml` | rarefaction-diagnostics settings for `Alldiagnose.py` (thresholds, the cell budget for "cells needed") |
| `viz.yaml`, `viz-diagnostics.yaml` | imagery specs; per-Kn specs are written from the manifest |

Per-Kn differences in the committed suite:

| Case | Mesh | Exit-cell parcel target | Verdict before running |
|---|---|---|---|
| `Kn100_np1x`, `Kn10_np1x`, `Kn1_np1x` | uniform 20 mm (`core_cell_size_m: 0.02`, `outer_expansion: 1.0`, `max_cells: 4100000`) | 2000 | production candidate |
| `Kn0p1_np1x` | inherited graded planner | 40 | production candidate |
| `Kn0p01_np1x` | inherited graded planner | 40 | unsuitable: mesh budget / resolution |
| `Kn0p001_np1x`, `Kn0p0001_np1x` | inherited graded planner | 40 | unsuitable; known inflow-surface limitation |

The verdict is computed by `feasibility.assess` from the derived case, not
typed: cell/λ₀ and Δt/τ₀ for the three low-Kn cases are recorded in the
manifest's `feasibility` block (the case README records about 2.69 / 26.2 / 270
and 0.146 / 1.42 / 14.6). The criteria are not relaxed to make a case look
resolved. A *candidate* still needs measured occupancy, stationarity and
mesh/time-step evidence after it runs. The direct-exit inflow surface at
Kn ≤ 0.001 is a known limitation deferred to future work; the interface track
is the attempt to work around it.

Schema details: [configuration reference](../reference/configuration.md#cai2012-kn-extension).

## Generation, exactly

`generate_cases.py` stages every case in a scratch directory
(`.staging-*` beside `study.yaml`, removed on exit) using `cai2012`'s own
`apply_case_parameters`, compares it with the disk, and classifies it:

| Status | Meaning | Default action | With a flag |
|---|---|---|---|
| `new` | nothing there yet | create | — |
| `identical` | same inputs, same scripts | keep, byte for byte — mesh, dictionaries and results included | — |
| `stale` | physics identical; a runner script or the `meta` block differs from the template | keep, and report | `--refresh` copies the new scripts/meta in, **even into cases with results** |
| `conflict` | an input in `case.yaml` differs | **refuse**: nothing is written, the manifest is unchanged, exit 1 | `--regenerate` replaces it only if it has no solver output; `--overwrite` replaces it **even if it has results, destroying them** |

| Flag | Writes | Destroys |
|---|---|---|
| `--dry-run` | only the temporary `.staging-*` directory; prints the plan and the feasibility table | nothing |
| (none) | new cases, `manifest-full-suite.yaml` | nothing |
| `--refresh` | changed runner scripts and `meta` into stale cases | nothing physical |
| `--regenerate` | replaces conflicting cases **without** solver output | those unrun cases |
| `--overwrite` | replaces conflicting cases **with** solver output | their results — the script prints `DESTROYING` per case |

"Solver output" means a time directory in the case or in any `processor*`, or a
`log.dsmcFoam`. Exit status is 0 when the tree matches the study **even if some
cases are unsuitable** — that is a recorded verdict, not an error — and 1 on a
conflict, a hard DSMC check failure, or a value that should be shared and is
not.

## Cheap checks

All from `cases/cai2012-kn-extension/` (or give the script path from the
repository root). These need Python only.

```bash
python generate_cases.py --dry-run         # plan + feasibility table; nothing kept
python validate_cases.py --inputs-only     # read-only: case.yaml inputs against the study
python -m plumetools.cai2012ext.policy --study study.yaml --all   # every case path, policy enforced
python Allscale.py inventory               # writes results/scaling/full-suite/inventory.{json,csv}
```

`validate_cases.py` without `--inputs-only` also reads the generated
dictionaries, so it needs `./AllmeshCases --dict-only` first; `-v` prints every
check, `--manifest PATH` selects another manifest.

**`smoke_cases.py` is a preparation step, not a read-only check.** It runs
`generate_cases.py`, writes all seven cases' OpenFOAM dictionaries with the
`Allmesh --dict-only` writers, runs the independent verification, checks every
runtime-policy guard, re-runs generation to confirm nothing changes, checks
scaling registration, and writes `results/full-suite-preparation/smoke.json`
(or `--report PATH`). It never meshes, runs DSMC or submits a job, and it
**refuses** a tree that already holds a mesh or solver output — use
`validate_cases.py` there.

```bash
python smoke_cases.py                      # in a tree with no meshes or results
```

!!! warning "A passing smoke check is not a production check"
    It shows that the inputs are complete and mutually consistent and that the
    guards behave. It says nothing about whether a case is resolved: the three
    low-Kn cases pass it and remain unsuitable.

### Recorded preparation check (2026-10-26)

Moved from the case README, as recorded there: the local preparation produced
seven case configurations and 70 dictionary files (about 0.8 MB), passed 84
independent verification checks, and confirmed idempotence, all seven
runtime-policy gates and scaling registration, with zero mesh trees and no
solver output. The focused policy, generation, HPC-wrapper, diagnostics and
scaling tests then reported 237 passed and 3 skipped (the skips need Unix file
locking). That validates input preparation and bounded synthetic numerics; it
is not an OpenFOAM/DSMC run, and no scheduler command was run.

## Production run

### On a workstation

```bash
python generate_cases.py
./AllmeshCases                 # or --dict-only; --list / --index N / --case NAME; --keep-going
python validate_cases.py
cd Cases_full_suite/Kn1_np1x && ./Allrun     # or --serial / -np N; resumes like cai2012
cd ../.. && ./AllpostCases     # see below
```

The study runners read the case list through
`python -m plumetools.cai2012ext.policy … --all`, which enforces the suite's Kn
policy and lists every case including non-production ones.

### On PBS

Every job script sources `hpc_env.sh`, found in `./`, `../` or `../../` of
`$PBS_O_WORKDIR`. **That file is the owner's cluster environment** — Spack and
Lmod set-up, specific module versions, an OpenFOAM build under a home
directory, a Python virtual environment, and a ParaView path — and the job
headers name the owner's queues and notification address. Replace them before
submitting anywhere else; see [PBS jobs](../guides/pbs.md).

```bash
cd cases/cai2012-kn-extension
python smoke_cases.py                   # in a fresh tree only
qsub runMeshAll.sh                      # generate, validate --inputs-only, AllmeshCases, validate
python validate_cases.py
for case in $(python -m plumetools.cai2012ext.policy --study study.yaml --all); do
    (cd "$case" && qsub runCase.sh)     # one job per case, one node, dsmc.n_subdomains = 192 ranks
done
```

`runCase.sh` per case: `hpc guard` (the suite's Kn policy, from the case's own
`case.yaml`), `hpc label` (re-assess; write a `NON-PRODUCTION` marker into an
unsuitable case and at the top of `output.log` — it still runs), `hpc finished`
(whether the case was already at `endTime`), `./Allrun >> output.log 2> error.log`,
then `hpc record` (the wall clock into `results/cost-model.yaml`, under a
`fcntl` lock, unless the case was already finished). A job that runs out of
walltime is **resubmitted by hand** with the same `qsub`; `./Allrun` continues
from the newest complete checkpoint. Nothing resubmits itself.

Other job scripts: `baseCase/runMesh.sh` (guard, then `./Allmesh`),
`runPost.sh` (`AllpostCases --case <case> --no-viz --no-analyse`),
`runViz.sh` (`AllpostCases --case <case> --viz-only`), and at the study level
`runAnalyse.sh` (`AllpostCases --analyse-only`) and `runScaling.sh`.

!!! note "A stale comment in `runCase.sh`"
    Its header says the guard refuses Kn = 0.01, 0.1 or 100. That is the
    historical four-Kn policy; under `full_decades` the guard permits all seven
    (`plumetools.cai2012ext.hpc.guard`). The code is authoritative.

### Post-processing

`./AllpostCases` runs, per case with output: `Allpost` (the `cai2012` centreline,
contour and metrics, byte-identical script), `Alldiagnose.py` (rarefaction
diagnostics), `Allcompare.py` (r/D = 10 angular density and mass flux,
off-axis velocity, temperature, pressure, with exclusion counts),
`Allaudit.py` (occupancy per region, per-frame centreline error — the health
script, byte-identical), `Allconverge.py` (far-field stationarity across
frames); then imagery (one `pvpython` call per Kn) and `analyse.py`. Options:
`--case NAME` (repeatable), `--no-viz`, `--viz-only`, `--no-analyse`,
`--analyse-only`, `--no-diagnostics`, `--stride N`, `--keep-going`. Cases with
no solver output are skipped and reported.

!!! danger "Destructive options"
    `--overwrite` (above). `./AllcleanCases --purge` deletes the **active**
    cases tree and its manifest (it reads the active manifest, so historical
    trees are not touched). Preview with `./AllcleanCases --dry-run`.

## Outputs

| Path | Written by | Contents |
|---|---|---|
| `manifest-full-suite.yaml` | `generate_cases.py` | every case's derived inputs, `feasibility` (status, resolution, time step, cost, `unresolved` reasons), `production_status`, the inheritance report, git provenance, stray directories |
| `<case>/NON-PRODUCTION` | `hpc label` | the unresolved criteria, for an unsuitable case |
| `<case>/output.log`, `error.log` | `runCase.sh` | the label, then `Allrun`'s output |
| `<case>/results/` | `Allpost` and friends | `centerline.csv`, `density_plane.csv`, `metrics.yaml`; `centerline_extended.csv`, `meridional_plane.csv`, `comparison.yaml` (Allcompare); `diagnostics_summary.yaml`, `centreline_diagnostics.csv`, `centreline-diagnostics.png` (Alldiagnose); `occupancy.yaml`, `convergence.yaml` (Allaudit); `stationarity.yaml` (Allconverge) |
| `<case>/<time>/` | `Alldiagnose.py` | diagnostic volume fields (`mfp`, `KnD`, `KnCellMax`, `mfpResolved`, …) for ParaView |
| `results/full-suite/` | `analyse.py` | `study-table.csv`, `study-table.md`, `knudsen-comparison.png` |
| `results/cost-model.yaml` | `hpc record` | measured wall clock per case (tracked) |
| `results/full-suite-preparation/smoke.json` | `smoke_cases.py` | the preparation check's result |

Every artefact from an unsuitable case carries its label. Schemas and the
normalisation of each compared quantity: [Outputs](../reference/outputs.md#cai2012-kn-extension)
and [the extension's §7–8](../cai2012-kn-extension.md#7-what-is-compared-and-how).

**Judging the output.** Read `production_status` and `unresolved` first. For a
candidate, check measured occupancy (`occupancy.yaml`) against the target, the
stationarity verdict, and the diagnostics' fraction of cells coarser than λ/3
before reading any error against Cai. Cai's published errors are not acceptance
tolerances. A zero exit status from any step, or a clean smoke check, does not
establish production suitability.

## Resources

Configured: `dsmc.n_subdomains: 192` (one node per case); the PBS headers
request `nodes=1:ppn=192` and 120 h walltime for run/post/viz jobs and 1–4 h on
`shortq` for mesh, analysis and scaling jobs. Estimated parcels, cells and steps
per case are in the manifest's `feasibility.cost` (a cost model fitted from
`cai2012-health` measurements, or a prior; `--dry-run` prints which). No
measured run of the seven-decade suite is recorded in this repository.

## Figure placeholder

!!! figure "Figure placeholder — the Knudsen comparison"
    **Subject:** `results/full-suite/knudsen-comparison.png` — centreline n/n₀, U, T and r/D = 10 density and mass flux for every run Kn, against Cai's collisionless curve.
    **Source run:** at least the four candidate cases solved and post-processed; label any non-production case in the caption.
    **Plotting command:** `./AllpostCases --analyse-only` (or `python analyse.py`), verified in source.
    **Destination:** `docs/assets/figures/kn-extension-knudsen-comparison.png`
    **Suggested caption:** "Departure from Cai's collisionless solution across Knudsen number (non-production cases marked)."

## Related

- [The Cai 2012 Knudsen extension](../cai2012-kn-extension.md) — scope, inheritance, feasibility, comparisons, metrics, open problems
- [Interface track](cai2012-interfaces.md) · [Field scaling](cai2012-scaling.md)
- [PBS jobs](../guides/pbs.md) · [Restart and recovery](../guides/restart-and-recovery.md)
- API: [`cai2012ext.generation`](../reference/api/cai2012ext/generation.md), [`cai2012ext.policy`](../reference/api/cai2012ext/policy.md), [`cai2012ext.feasibility`](../reference/api/cai2012ext/feasibility.md), [`cai2012ext.hpc`](../reference/api/cai2012ext/hpc.md)

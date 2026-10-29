# Known gaps and defects

What is missing, broken or inconsistent, found while writing this
documentation (2026-10-28, commit `8585148`) or recorded earlier by the owner.
**None of these was fixed here**: the documentation work changed no code,
test, input or runner. Each item says where it is and what it affects.

## Blockers to running a workflow

| Where | Gap | Effect |
|---|---|---|
| `cases/3d-inflow` | `constant/polyMesh` is gitignored and absent; generating one needs OpenFOAM | `Allrun --no-solve` cannot run in a fresh checkout without OpenFOAM or a history checkout |
| `cases/3d-inflow/Allrun`, `Allmesh`, `.gitignore`, `tests/regression/test_3d_inflow.py` | name the tag `pointwise-mesh-3d-inflow`, which does not exist locally or on the remote | their recovery command fails; the working equivalent is `git checkout f83233f^ -- cases/3d-inflow/constant/polyMesh` (verified with `git ls-tree`) |
| `cases/ARCHIVE.md`, `tests/regression/README.md`, `tests/regression/capture_golden.py` | name the tag `pre-refactor-baseline`, which does not exist | the commit is `7997ccc` |
| `cases/iss_solar_panels*`, `solar_panel_particle_resolution_study/Fnum_*` | still use the MNF `dsmcFoam+`, per-case legacy `processInflowData.py` with hard-coded inputs, no run script | not runnable with the supported stack; their `clean.sh` is the unguarded `rm -r 0* …` (RP-01) |
| `solar_panel_particle_resolution_study` | every plotting script reads a CSV that was never committed (RP-06) | the study's figures cannot be regenerated |
| `cases/cai2012-health`, `cases/markelov1999` | generation needs CaseFoam (Markelov's `--dry-run` too) | blocked where CaseFoam will not install |
| `cases/cai2012-kn-extension/hpc_env.sh`, every `#PBS` header | the owner's site: module versions, home-directory paths, queues, notification address | must be replaced before use elsewhere ([PBS jobs](../guides/pbs.md)) |
| `cases/cai2012-kn-extension/scaling-registration.yaml` | registers `Cases_run1_graded/…`, a local untracked archive | a fresh clone's inventory reports those cases missing (correctly) |
| `Allsetup` | builds `libplumeDsmcBoundaryModels.so` but not `libhybridPlumeThermo.so`, and its closing hints name only markelov1999/cai2012 | the hybrid family needs `applications/thermophysicalModels/Allwmake` by hand |

## Scientific limits recorded by the families

| Where | Limit |
|---|---|
| Kn extension | Kn = 0.01, 0.001, 0.0001 unsuitable for production as configured; the direct-exit inflow surface at Kn ≤ 0.001 is a known limitation deferred to future work; no production run of the seven-decade suite is recorded |
| Kn-extension interface | the overlap validation cannot run as declared (no admissible sphere at Kn = 0.01 for B = 0.02–0.05); two decisions (core angle, source-validity radius) are left open; Markelov's angular law under-predicts axis density for Cai's supersonic exit by a factor 0.36–0.45 |
| hybrid | every recorded interface verdict is unsuitable; verification and production profiles prepared, not run; the `dsmcFoam` restart insertion burst is unexplained; the CFD background floor leaks 1–4 % of the source flow inward in the smoke runs |
| r-4d | under-resolved at the exit by necessity; the shipped cone mesh not run to its end time; the HPC variants never run |
| cai2012 | Kn = 0.01 does not meet Cai's Δx/λ₀ = 1 within the cell budget |
| markelov1999 | steady state not demonstrated at the shipped end time; occupancy target met at the cylinder only |
| virtual collisions | DSMC comparisons exploratory: spatial convergence of the comparison estimators not established |
| scaling | no whole-field DSMC fields exist to fit |
| 3d-inflow | HA-06 (`dsmcAxisymmetric` on a 3D mesh) and HA-07 (`type Unspecified` patch) unresolved; post-processing not wired (RP-17); outer boundary reflects under standard `dsmcFoam` |
| every family | no result physically validated |

## Defects and inconsistencies found in code, tests and comments

| Where | Finding |
|---|---|
| `tests/unit/test_viz_dissect.py::test_a_missing_cell_image_is_refused_before_ffmpeg_runs` | fails on a machine without ffmpeg: `dissect` checks for ffmpeg before checking the cell images, so the test gets "ffmpeg is not on PATH" instead of the expected "do not exist" (observed on native Windows, 2026-10-28). Either the order of checks or the test's assumption of ffmpeg is wrong |
| `cases/cai2012-kn-extension/baseCase/runCase.sh` (header comment) | says the guard refuses Kn = 0.01, 0.1, 100; under `meta.knudsen_suite: full_decades` the guard permits all seven (`plumetools.cai2012ext.hpc.guard`) |
| `cases/cai2012-kn-extension/AllmeshCases` (header comment) | says the policy listing refuses a manifest naming Kn = 0.01, 0.1 or 100 — true of the historical four-Kn suite only |
| `plumetools/markelov1999/resolution.py` `mean_free_path` docstring | writes the omitted VHS correction as (T_ref/T)^(ω−½); on λ it is (T/T_ref)^(ω−½), as `cai2012.gas` applies. The returned value is unaffected (the function deliberately applies no correction) |
| `plumetools/cai2012/dictionaries.py` | writes `NoBinaryCollision` when collisions are disabled; OpenFOAM v2512 registers that model as `none` (reported in [the hybrid page](../hybridplume.md#2-backend-decision); `nozzleplume.dictionaries` maps it). Not reachable with the committed configuration, which requires `VariableHardSphere` |
| `docs/cai2012-kn-extension.md` § "Observed in the frozen packages" | further defects in `cai2012` and other frozen packages, reported there earlier |
| `cases/3d-inflow/Allrun` header | describes the initialiser as `dsmcInitialise+`; the script derives it from `controlDict` (`dsmcFoam` → `dsmcInitialise`) |
| `cases/cai2012-virtual-collision/README.md` (before this change) | its `smoke` example wrote to `results/smoke/`, rewriting tracked files; the rewritten README and tutorial now say so |
| `.gitignore` | ignores `hybridPlumeThermo`'s build products but not `plumeFieldInflow`'s (`lnInclude/`, `Make/linux64*/`), which therefore show as untracked after a build |
| `cases/iss_solar_panels_3/iss_solar_panels_2/` | a nested copy of another case, with tracked `__pycache__` files |

## Documentation reconciled in this change

Stale statements corrected or superseded: the root README's four-Kn summary and
test counts, the case READMEs' solver name (3d-inflow) and resume description
(now per-rank checkpoints), `docs/environment.md`'s "`plumetools.run` is the
intended home" for provenance (per-family provenance now exists in
`cai2012ext.provenance` and `hybridplume.provenance`; there is no general
`plumetools.run`), and the missing-tag recovery commands in `docs/mesh.md`. The
detailed content of the cai2012-health and r-4d READMEs moved, preserved, to
[cai2012-health, explained](../cai2012-health.md) and
[r-4d, explained](../r-4d-case.md).

## Documentation tooling

The site depends on Material for MkDocs, which is in maintenance with an
announced end of life on 5 May 2027, and on MkDocs 1.6 (2.0 is not compatible).
Equations need a network connection to a CDN for MathJax. Module pages are
generated, so links to them work only in the built site. See
[Writing the documentation](documentation.md#tooling-status).

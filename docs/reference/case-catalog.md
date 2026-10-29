# Case catalog and coverage matrix

Every case directory under `cases/`, what it is, how far its evidence goes, and
what stands between it and a run. "Active" is the repository's label for work
that is maintained; it does **not** mean runnable — see the last two columns.

## Evidence levels

| Level | Meaning here |
|---|---|
| **Implemented** | the code and scripts exist and are read by this documentation |
| **Mechanically tested** | pure-Python tests in `tests/unit` (and the regression golden) exercise it; no OpenFOAM involved |
| **Solver-tested** | `tests/openfoam` (marker `needs_openfoam`) runs it through real OpenFOAM tools, with deliberately short runs |
| **Executed (recorded)** | the owner recorded a real run in a document; not re-run here |
| **Physically validated** | compared with experiment or an independent code with an accepted tolerance. **No case in this repository has reached this level** |

What each test establishes, and does not: [What the tests establish](../explanation/evidence.md).

## Active families

| Family | Directory | Model | Evidence | Runnable today with the supported stack? | External prerequisites | Known blockers / limits | Tutorial |
|---|---|---|---|---|---|---|---|
| Legacy reference | `cases/3d-inflow` | legacy source flow (`plumetools.sourceflow`) | mechanically tested; bit-exact regression golden vs the pre-refactor script | **only after obtaining a mesh** | OpenFOAM v2512 (`blockMesh`, `snappyHexMesh`, `dsmcFoam`) | `constant/polyMesh` is gitignored and absent; the restore tag in older text does not exist (use `git checkout f83233f^ -- …`); outer boundary reflects; post-processing not wired (RP-17); HA-06, HA-07 unresolved | [3d-inflow](../tutorials/3d-inflow.md) |
| Cai collisionless plume | `cases/cai2012` | 3-D `dsmcFoam`, nozzle exit as inlet | mechanically tested; solver-tested (short runs); recorded Kn 100 and 0.1 runs | yes | OpenFOAM v2512 + `libplumeDsmcBoundaryModels.so` | Kn = 0.01 under-resolved against Cai's Δx/λ₀ = 1 | [cai2012](../tutorials/cai2012.md) |
| Statistical health | `cases/cai2012-health` | `cai2012` case, two axes swept | mechanically tested; recorded cheapest-pair runs | yes, **with CaseFoam** for generation | as cai2012, plus CaseFoam, ffmpeg, `pvpython` | generation needs CaseFoam (not installed here); disk is the binding resource | [cai2012-health](../tutorials/cai2012-health.md) |
| Kn extension (direct) | `cases/cai2012-kn-extension` | `cai2012` cases via `cai2012ext` | mechanically tested; preparation recorded 2026-10-26; no suite run recorded | yes (generation, mesh, run); PBS scripts need site edits | as cai2012; a PBS cluster for the committed job scripts | Kn ≤ 0.01 unsuitable for production; Kn ≤ 0.001 inflow-surface limitation deferred; `hpc_env.sh` is the owner's site | [Kn extension](../tutorials/cai2012-kn-extension.md) |
| Kn extension (interface) | `cases/cai2012-kn-extension` (`Interface/`) | truncated domain + downstream angular source | mechanically tested; one development smoke run recorded | smoke case only | as cai2012; `snappyHexMesh` | overlap validation infeasible as declared; production cases refused at the run step | [Interface track](../tutorials/cai2012-interfaces.md) |
| Kn extension (scaling) | `cases/cai2012-kn-extension` (`Allscale.py`) | empirical coordinate maps of Cai's field | mechanically tested; synthetic smoke | yes (Python only) | none beyond Python | no whole-field DSMC fields exist to fit yet | [Field scaling](../tutorials/cai2012-scaling.md) |
| Virtual collisions | `cases/cai2012-virtual-collision` | reduced pair-event model vs Cai and DSMC | mechanically tested; committed bounded smoke result | yes (Python only) | DSMC references only for the comparison step | comparisons exploratory (spatial convergence not shown); the four saved export sets are local, not in the repository | [Virtual collisions](../tutorials/cai2012-virtual-collision.md) |
| Hybrid CFD → DSMC | `cases/cai2012-hybrid` | `rhoCentralFoam` precursor → interface → `dsmcFoam` | mechanically tested (analytical pieces); smoke runs recorded | yes | OpenFOAM v2512 + **both** custom libraries | every recorded interface verdict unsuitable; verification/production profiles prepared, not run; restart insertion burst unexplained | [Hybrid](../tutorials/cai2012-hybrid.md) |
| Plate in a cylinder wake | `cases/markelov1999` | paper-faithful source flow + `plumeFieldInflow` | mechanically tested; solver-tested | yes, **with CaseFoam** for generation | OpenFOAM v2512 + library + CaseFoam | steady state not demonstrated; occupancy target met at the cylinder only | [markelov1999](../tutorials/markelov1999.md) |
| Nozzle plume (R-4D) | `cases/r-4d` (`baseCase`, `baseCaseHPC`, `baseCaseFreePlumeHPC`) | `nozzleplume`: direct exit state, optional plate | mechanically tested; solver-tested (collisionless conservation); one full-length run recorded on an earlier mesh | yes | OpenFOAM v2512 + library; `pvpython` for imagery | under-resolved at the exit by necessity; the shipped cone mesh not run to its end; HPC variants never run | [r-4d](../tutorials/r-4d.md) |
| ISS solar panels | `cases/iss_solar_panels`, `_2`, `_3` | per-case legacy `processInflowData.py` | none | **no** | the MNF fork's `dsmcFoam+` | fork not available; no run script; inputs hard-coded; unguarded `clean.sh` (RP-01) | [ISS](../tutorials/iss-solar-panels.md) |
| Particle-weight study | `cases/solar_panel_particle_resolution_study/Fnum_*` | as ISS | none | **no** | `dsmcFoam+`; the CSV exports | as ISS, plus every plotting input CSV missing (RP-06) | [Particle resolution](../tutorials/particle-resolution.md) |

No other case directory exists under `cases/` in the current tree.

## Templates and generated trees

| Path | Tracked? | Role |
|---|---|---|
| `cases/{cai2012,cai2012-health,cai2012-kn-extension,markelov1999}/baseCase/` | yes | study templates; generators clone them |
| `cases/cai2012-kn-extension/baseCaseInterface/` | yes | overlay for interface cases |
| `cases/r-4d/baseCase*/` | yes | three complete case templates |
| `Cases/`, `Cases_full_suite/`, `Interface/`, `manifest*.yaml`, `results/` | no (gitignored) | generated; reproducible from the templates and study files |
| `cases/cai2012-health/results/cost-model.yaml`, `cases/cai2012-kn-extension/results/cost-model.yaml` | yes | measured wall clock — data, kept |
| `cases/cai2012-virtual-collision/results/smoke/` | yes | the committed bounded smoke result |
| `cases/cai2012-hybrid/runs/` | no | hybrid case directories |

## Archived families

Frozen historical record, never refactored or fixed:
`1d`, `2d-planar`, `2d-wedge`, `caseFoamEx`, `wake-cylinder`. They are not
supported quickstarts. What each was for and what is known to be wrong in it:
[Archive](../archive/index.md).

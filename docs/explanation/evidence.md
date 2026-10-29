# What the tests establish

A green test suite is a statement about code, not about physics. This page
says which kind of statement each tier makes, what was measured when this
documentation was written, and what no test in the repository establishes.

## The tiers

| Tier | Location | Runs | Establishes | Does not establish |
|---|---|---|---|---|
| 0 — unit | `tests/unit` | pure Python, no OpenFOAM (`pytest -m "not needs_openfoam"`) | that derivations, parsers, writers, generators, policies and analytical solutions do what their docstrings say, against independent oracles (closed forms, brute-force counts, round trips, synthetic meshes) | that OpenFOAM accepts the files, or that any solution is right |
| regression | `tests/regression` | pure Python; **needs the gitignored historical 3d-inflow mesh** (skips without it) | that `plumetools` reproduces the pre-refactor legacy inflow script **bit for bit** on that mesh | correctness: the frozen defects SM-01 … SM-10 are inside the golden ("unchanged", never "right") |
| 2 — OpenFOAM | `tests/openfoam` | `pytest -m needs_openfoam`, with OpenFOAM v2512 and the inflow library | that the generated dictionaries parse; meshes build and are the declared geometry; `topoSet` selects the predicted nozzle faces; `plumeFieldInflow` loads, injects on the named patch only and at the Maxwellian rate; open boundaries delete molecules; resumes behave; and, for the nozzle plume, that a collisionless jet onto a plate conserves momentum and energy | convergence, resolution, or agreement with experiment: the runs are deliberately a few hundred steps |

The recorded executions in the case documents ("executed (recorded)" in the
[case catalog](../reference/case-catalog.md)) sit above Tier 2 and below
validation: real runs, measured and reported with their noise floors, but not
compared with an independent reference under an accepted tolerance.

**No result in this repository is physically validated.** Cai's published
errors are carried as reference values, not tolerances; the Markelov family has
not been compared with DAC or experiment; the R-4D case has not been compared
with measurement or another code.

## Measured for this documentation (2026-10-28)

Run on native Windows (Python 3.12.10, numpy 2.4.4, scipy 1.18.0), at commit
`8585148`, without OpenFOAM, ffmpeg or CaseFoam:

```bash
python -m pytest tests/unit tests/regression -q -m "not needs_openfoam" -p no:cacheprovider
```

| Result | Count |
|---|---|
| passed | 2228 |
| skipped | 37 — among them all 16 regression tests (no historical 3d-inflow mesh), the CaseFoam-dependent `markelov1999` study tests, tests that drive the bash runners (on native Windows `bash` is WSL's and cannot open a Windows path; they run under WSL), and the `fcntl` file-locking tests (Unix only) |
| failed | 1 — `test_viz_dissect.py::test_a_missing_cell_image_is_refused_before_ffmpeg_runs`, which expects the missing-image error but gets "ffmpeg is not on PATH" on a machine without ffmpeg ([Known gaps](../development/known-gaps.md)) |

Wall time was about six minutes. The OpenFOAM tier was **not** run for this
documentation. Earlier counts in older pages (for example "475 tests, ~4 s") are
out of date; the refactor record's counts belong to its own commit.

## Per family

| Family | Tier 0 modules | Tier 2 | What the oracles are |
|---|---|---|---|
| legacy 3d-inflow | `test_sourceflow`, `test_config`, `test_geometry`, `test_blockmesh`, `test_snappy`, `test_mesh_parsing` + `tests/regression/test_3d_inflow.py` | — | equations E1–E8 term by term; the pre-refactor script on the historical mesh |
| cai2012 | `test_cai2012_*` (analytical, config, foamio, gas, geometry, health, inflow, mesh, particles, post, study) | `test_cai2012_openfoam.py` | closed-form limits of the collisionless solution; quadrature vs closed form on the axis to 1e-10; brute-force nozzle face counts; Kn ↔ n₀ round trips; byte identity of the health copies of the runner scripts |
| Kn extension | `test_cai2012ext_*` | (via the cai2012 tier) | policy and inheritance tables; idempotent generation against staged copies; feasibility verdicts; interface physics (∫ j dΩ = ṁ, core fraction, Bird's parameter by finite difference, Stewart's equations); overlap interval logic; synthetic-only scaling checks, labelled as such |
| hybrid | `test_hybridplume_*` | — | half-range moments vs quadrature, NSF flux signs, the collisionless reference vs Cai's closed form (1e-13), exact crossing accounting on synthetic clouds |
| virtual collisions | `test_virtualcollision_*` | — | physical limits and invariants of the VHS kernel and the reduced model; geometry recovery against independent counts |
| markelov1999 | `test_markelov_*` | `test_markelov_openfoam.py` | the paper's model equation by equation; flux verification by two independent routes |
| nozzle plume (r-4d) | `test_nozzleplume_*` | `test_nozzleplume_openfoam.py` | geometry and mesh-plan rules against brute-force counts; plate loads against surfaces with known answers; conservation of a collisionless jet onto a diffuse plate |
| shared | `test_rarefaction`, `test_foamio_*`, `test_casemanifest`, `test_iteration`, `test_allrun_restart`, `test_streaming_io`, `test_primitives`, `test_viz`, `test_viz_dissect` | — | diagnostics against closed forms; checkpoint decisions on synthetic rank layouts and with stub OpenFOAM tools; streamed reads byte-identical to whole-file reads |

## Reading a test claim

- **Independent oracle or not.** The project keeps checks independent of the
  code under check (`cai2012.audit`, the `verify` modules, the tests' oracles
  compute what they check without calling the code they check). A test that
  re-derives a value through the same function proves only determinism.
- **Synthetic data are labelled.** The scaling and virtual-collision tests
  state in their docstrings that their fields are not DSMC evidence; the same
  holds for the committed `results/smoke/`.
- **Mechanics vs physics.** A smoke run that completes shows the pipeline
  works. Whether the answer is right needs resolution, steady state, enough
  parcels and an independent reference — which the case pages discuss case by
  case.
- **Historical output is data, not validation.** Committed images and the
  archived cases' plots have provenance problems recorded in
  [the archive](../archive/index.md).

Where each model's equations live in code, and which tests cover them:
[Explanation index](index.md).

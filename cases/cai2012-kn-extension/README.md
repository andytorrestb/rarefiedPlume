# cases/cai2012-kn-extension — Cai's plume across seven decades of Knudsen number

The active study (`meta.knudsen_suite: full_decades` in [`study.yaml`](study.yaml))
generates **Kn = 100, 10, 1, 0.1, 0.01, 0.001, 0.0001** for Cai's circular argon
nozzle (D = 0.2 m, S₀ = 2, T = 300 K) into `Cases_full_suite/` and
`manifest-full-suite.yaml`, with `cai2012-health`'s highest-resolution settings
and Cai's 10⁴ t₀ startup. A study without that declaration keeps the historical
four-Kn policy (1e-4, 1e-3, 1, 10 in `Cases/` and `manifest.yaml`).

Before running, Kn ≥ 0.1 are production *candidates*; **Kn = 0.01, 0.001 and
0.0001 are unsuitable for production** as configured, and the direct-exit inflow
surface at Kn ≤ 0.001 is a known limitation deferred to future work. Unsuitable
cases still run, labelled NON-PRODUCTION. No production run of this suite is
recorded, and nothing is validated.

```bash
python generate_cases.py --dry-run          # plan + feasibility; keeps nothing
python validate_cases.py --inputs-only      # read-only input checks
python smoke_cases.py                       # PREPARATION: writes the tree, dictionaries, smoke.json
qsub runMeshAll.sh                          # cluster: generate, validate, mesh
for c in $(python -m plumetools.cai2012ext.policy --study study.yaml --all); do (cd "$c" && qsub runCase.sh); done
./AllpostCases                              # after the runs
```

`hpc_env.sh` and the `#PBS` headers are the owner's cluster settings; replace them
before submitting elsewhere. `--overwrite` and `AllcleanCases --purge` destroy
results.

- Direct cases: [`docs/tutorials/cai2012-kn-extension.md`](../../docs/tutorials/cai2012-kn-extension.md)
- Interface track (`generate_interfaces.py`, `Interface/`):
  [`docs/tutorials/cai2012-interfaces.md`](../../docs/tutorials/cai2012-interfaces.md)
- Field scaling (`Allscale.py`): [`docs/tutorials/cai2012-scaling.md`](../../docs/tutorials/cai2012-scaling.md)
- Method, inheritance, metrics, open problems:
  [`docs/cai2012-kn-extension.md`](../../docs/cai2012-kn-extension.md)

The preparation recorded on 2026-10-26 (seven configurations, 70 dictionary
files, 84 verification checks, idempotence, all seven policy gates, scaling
registration; no meshes, no solver output, no scheduler command) is described in
the tutorial's checks section; it validates input preparation, not a DSMC run.

The `--all` listing deliberately includes non-production cases. Without it,
policy listing returns production candidates and reports the refused cases.
Each case has `runMesh.sh`, `runCase.sh`, `runPost.sh`, and `runViz.sh`.
Completed runs are not automatically resubmitted. Post-process with each case's
`runPost.sh`, or use the study-level `AllpostCases` after reconstruction.
`analyse.py` writes the active study summary beneath `results/full-suite/`.

A `runCase.sh` that runs out of walltime during `reconstructPar` leaves a partial
time directory, and post-processing then stops on `no rhoNMean`. Resubmitting
does not repair it, because `reconstructPar -newTimes` skips existing directories.
`./AllreconstructCases --check` reports such cases; `qsub runReconstruct.sh`
finishes them.

## Scaling integration

The scaling inventory follows `study.yaml`'s `manifest_name` and `cases_dir`.
All seven new direct-exit cases are registered automatically. `scaling-registration.yaml`
adds the explicitly named local graded archives; it is not needed to register
new full-suite cases. Historical results remain separate by name and source.

```bash
python Allscale.py inventory
# Once reconstructed fields include actual C, V and density arrays:
python Allscale.py fit --full
python Allscale.py report
# Analysis-only PBS template, submitted manually if desired:
msub runScaling.sh
```

New scaling outputs use `results/scaling/full-suite/`. Missing results are
skipped honestly. Available low-Kn fields can enter exploratory fits but remain
excluded from suitable-reference fits while source/numerical limitations remain.
There is currently no whole-field scientific fit to the newly prepared suite.
Full fitting details and portable-field formats are in
[the scaling documentation](../../docs/cai2012-scaling.md). Existing comparison
normalizations and estimator definitions are in
[the study documentation](../../docs/cai2012-kn-extension.md).

## Preparation validation (2026-10-26)

The actual local preparation produced seven case configurations and 70
dictionary files (about 0.8 MB total), passed 84 independent verification checks,
and confirmed idempotence, all seven runtime-policy gates, and scaling
registration. There are zero generated mesh trees and no solver outputs.

The focused policy, generation, HPC-wrapper, diagnostics and scaling regression
run completed with **237 passed, 3 skipped in 62.70 seconds**. The three skips
require Unix file locking. This validates input preparation and bounded
synthetic numerical behavior; it does not claim a successful OpenFOAM/DSMC run
or a scientifically suitable low-Kn reference. No scheduler command was run.

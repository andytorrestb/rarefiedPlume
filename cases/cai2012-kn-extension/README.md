# Cai 2012 seven-decade Knudsen suite

The active study generates **Kn_e = 100, 10, 1, 0.1, 0.01, 0.001, 0.0001**
with the circular argon nozzle (`D_e=0.2 m`, `S0=2`, `T_e=300 K`), VHS/NTC
collisions enabled, and the same startup and averaging convention. Density is
derived from each Kn using the existing VHS mean-free-path convention.

The explicit `meta.knudsen_suite: full_decades` profile in `study.yaml` enables
all seven cases through generation, dictionary validation, runtime guards,
post-processing and scaling. Unmarked historical cases/manifests retain their
original four-Kn policy. No results from another study are imported.

Generated inputs live in `Cases_full_suite/` and `manifest-full-suite.yaml`.
The old `Cases/`, `manifest.yaml`, and graded archives are preserved.

## Prepare and smoke-test locally

From the repository root:

```powershell
python cases/cai2012-kn-extension/generate_cases.py --dry-run
python cases/cai2012-kn-extension/smoke_cases.py
python cases/cai2012-kn-extension/validate_cases.py
python cases/cai2012-kn-extension/Allscale.py inventory
```

`smoke_cases.py` generates the case files, writes all seven sets of OpenFOAM
input dictionaries using the same writers as `Allmesh --dict-only`, verifies
physical and numerical consistency, checks all seven runtime-policy guards,
confirms repeated generation preserves the files, and verifies that every case
is registered in scaling. It never builds a mesh, runs DSMC, or submits a job.
The machine-readable result is `results/full-suite-preparation/smoke.json`.
It refuses trees already containing a mesh or solver output; use the read-only
`validate_cases.py` for those.

Generation is idempotent. Conflicting physical inputs are reported without
replacement. Existing result-bearing cases are never overwritten unless a user
explicitly uses the destructive `--overwrite` option. No such option was used.
`AllcleanCases` reads the active manifest; its purge action removes only the
active cases tree and manifest, preserving shared/legacy results.

## Current numerical limitations

| Case | Exit-cell parcel target | Pre-run status |
|---|---:|---|
| `Kn100_np1x` | 2000 | production candidate |
| `Kn10_np1x` | 2000 | production candidate |
| `Kn1_np1x` | 2000 | production candidate |
| `Kn0p1_np1x` | 40 | production candidate |
| `Kn0p01_np1x` | 40 | unsuitable for production: mesh budget/resolution |
| `Kn0p001_np1x` | 40 | unsuitable for production; known inflow-surface limitation |
| `Kn0p0001_np1x` | 40 | unsuitable for production; known inflow-surface limitation |

**The direct-exit inflow surface at Kn_e <= 0.001 is a known limitation and is
left for future work.** The generated metadata, manifest, `KNOWN-LIMITATIONS`
file, scaling inventory and HPC job log retain that warning. No downstream
interface source was substituted. The nominal Kn still refers to the original
nozzle-exit state.

The numerical criteria were not relaxed: cell/lambda is about 2.69, 26.17 and
269.51 at Kn=.01, .001 and .0001; dt/tau is about .146, 1.42 and 14.60.
These dictionary files are mechanically valid but do not establish a resolved
physical simulation. The three cases retain NON-PRODUCTION labeling when run.
Pre-run candidates still require measured occupancy, stationarity and mesh/time
step evidence after HPC execution. Requested parcel targets differ by Kn and
are verified against each case's own weight; they are not physical invariants.

## Run later on HPC

These are manual commands for the cluster; no scheduler job was submitted by
this preparation. The existing PBS jobs inherit `hpc_env.sh` and use the case's
`dsmc.n_subdomains=192` for decomposition and MPI ranks.

```bash
cd cases/cai2012-kn-extension
python smoke_cases.py
# This next command actually meshes on HPC:
qsub runMeshAll.sh
```

After mesh preparation completes, validate the written dictionaries and review
the feasibility warnings. To submit all seven as separate jobs, including the
explicitly exploratory low-Kn cases:

```bash
python validate_cases.py
for case in $(python -m plumetools.cai2012ext.policy --study study.yaml --all); do
    (cd "$case" && qsub runCase.sh)
done
```

The `--all` listing deliberately includes non-production cases. Without it,
policy listing returns production candidates and reports the refused cases.
Each case has `runMesh.sh`, `runCase.sh`, `runPost.sh`, and `runViz.sh`.
Completed runs are not automatically resubmitted. Post-process with each case's
`runPost.sh`, or use the study-level `AllpostCases` after reconstruction.
`analyse.py` writes the active study summary beneath `results/full-suite/`.

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
qsub runScaling.sh
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

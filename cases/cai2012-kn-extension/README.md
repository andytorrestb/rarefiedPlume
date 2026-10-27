# Cai 2012 seven-decade Knudsen suite

The active study generates Kn_e = 100, 10, 1, 0.1, 0.01, 0.001, 0.0001.
The source is the circular argon nozzle (D_e=0.2 m, S0=2, T_e=300 K), with
VHS/NTC collisions enabled and the shared Cai startup/averaging convention.

`meta.knudsen_suite: full_decades` explicitly enables the seven values through
generation, dictionaries, runtime guards and direct-exit post-processing.
Unmarked historical studies/manifests/cases retain their four-Kn policy.
Generated inputs use `Cases_full_suite/` and `manifest-full-suite.yaml`.
Existing `Cases/`, `manifest.yaml`, archives and results are preserved.

## Prepare inputs

From this directory, using the project Python environment:

```bash
python generate_cases.py --dry-run
python generate_cases.py
bash AllmeshCases --dict-only
python validate_cases.py
```

These commands generate small inputs and dictionaries without a mesh, solver
run or scheduler submission. Generation is idempotent. Conflicting physical
inputs are reported rather than overwritten; result-bearing trees are retained.

## Source and numerical limitations

The direct-exit inflow surface at Kn_e <= .001 is a known limitation deferred
for future work. Generated metadata, the manifest, KNOWN-LIMITATIONS and HPC
logs keep that warning. No downstream interface source was substituted.

Kn=.01, .001 and .0001 remain unsuitable for production under the inherited
mesh/time-step requirements. Cell/lambda is approximately 2.69, 26.17 and 269.51;
dt/tau is .146, 1.42 and 14.60. Their dictionary files can verify mechanically
while the physical simulation remains unresolved. No numerical criterion was
relaxed. Kn>=.1 cases are pre-run candidates requiring measured quality evidence.
The exit-cell parcel target is 2000 for Kn>=1 and 40 for Kn<=.1; each case's
weight is verified against its own target.

## Later HPC execution

The existing PBS environment and 192-subdomain case jobs are retained. Submit
manually on HPC; no submission is part of input generation:

```bash
qsub runMeshAll.sh
# After meshing completes and numerical limitations have been reviewed:
python validate_cases.py
for case in $(python -m plumetools.cai2012ext.policy --study study.yaml --all); do
    (cd "$case" && qsub runCase.sh)
done
```

`--all` deliberately lists the non-production low-Kn cases as exploratory runs.
Without it, listing returns candidates and reports unsuitable cases. Post-process
with each case's runPost.sh or the study-level AllpostCases. The active study
summary is written under `results/full-suite/`. AllcleanCases scopes cleaning to
the active tree/manifest; shared and historical results are preserved by purge.

Comparison definitions are in [the study documentation](../../docs/cai2012-kn-extension.md).

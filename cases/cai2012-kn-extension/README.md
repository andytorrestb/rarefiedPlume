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

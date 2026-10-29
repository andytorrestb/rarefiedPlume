# cases/cai2012-hybrid — one-way CFD → DSMC free plume, full 3-D

Cai & Wang (2012)'s circular argon source (D = 0.2 m, S₀ = 2), computed as a
`rhoCentralFoam` precursor, a fixed continuum-breakdown surface, and standard
`dsmcFoam` + `plumeFieldInflow` beyond it. Everything is 3-D Cartesian. Smoke
runs are intentionally under-resolved; verification and production profiles are
prepared, not run; nothing is validated.

| file | what |
|---|---|
| `hybrid.yaml` | every input, labelled [PAPER] / [REPO] / [BIRD] / [ASSUMPTION] / [NUMERICAL]; profiles overlay `base` |
| `Allrun` | the OpenFOAM stages, resumable: `setup`, `cfd`, `post-cfd`, `dsmc`, `post-dsmc`, `crossings`, `all` |

```bash
# OpenFOAM v2512 sourced, BOTH libraries built in this environment:
#   applications/thermophysicalModels/Allwmake
#   applications/dsmcBoundaryModels/Allwmake
C=runs/kn0p01c                     # runs/ is gitignored
./Allrun $C setup --profile smoke-kn0p01-coarse -np 4   # Python only
./Allrun $C cfd -np 4
./Allrun $C post-cfd               # verdict in results/interface.json
./Allrun $C dsmc -np 4             # refuses an unsuitable interface unless --allow-unsuitable
./Allrun $C post-dsmc              # results/summary.json, *.png
```

Every stage resumes after an interruption; nothing is deleted (to start over,
remove the case directory). Batch scripts:
`python -m plumetools.hybridplume hpc <case> --scheduler pbs|slurm|none --np N`.

- Tutorial (profiles, stages, suitability, resume, outputs, HPC):
  [`docs/tutorials/cai2012-hybrid.md`](../../docs/tutorials/cai2012-hybrid.md)
- Method, decisions and measured smoke results: [`docs/hybridplume.md`](../../docs/hybridplume.md)

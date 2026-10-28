# cases/cai2012-hybrid — one-way CFD → DSMC free plume, full 3-D

Cai & Wang (2012)'s circular argon source (D = 0.2 m, S0 = 2), computed as a
`rhoCentralFoam` precursor, a fixed continuum-breakdown surface, and standard
`dsmcFoam` + `plumeFieldInflow` beyond it. Everything is 3-D Cartesian.
Method, decisions and measured results: [`docs/hybridplume.md`](../../docs/hybridplume.md).

| file | what |
|---|---|
| `hybrid.yaml` | every input, labelled [PAPER] / [REPO] / [BIRD] / [ASSUMPTION] / [NUMERICAL]; profiles overlay `base` |
| `Allrun` | the OpenFOAM stages, resumable: `setup`, `cfd`, `post-cfd`, `dsmc`, `post-dsmc`, `crossings`, `all` |

Profiles:

| profile | purpose | where |
|---|---|---|
| `smoke-kn0p01`, `smoke-kn0p01-coarse`, `smoke-kn1em3` | hybrid machinery, **intentionally under-resolved** | laptop |
| `smoke-kn0p01-coarse-full` | full DSMC on the same mesh, for the overlap comparison | laptop |
| `smoke-collisionless` | collisionless 3-D DSMC vs the analytical reference | laptop |
| `verification-*` | Cai Kne 100 / 0.1 / 0.01 at verification resolution | HPC |
| `production-kn1em3`, `production-kn1em4` | the extensions | HPC |

```bash
# OpenFOAM v2512 sourced, both libraries built in THIS environment:
#   applications/thermophysicalModels/Allwmake
#   applications/dsmcBoundaryModels/Allwmake
C=runs/kn0p01c                     # runs/ is gitignored
./Allrun $C setup --profile smoke-kn0p01-coarse -np 4
./Allrun $C cfd -np 4              # mesh + precursor (mesh only for full-DSMC profiles)
./Allrun $C post-cfd               # cfd-report, breakdown, interface (verdict in results/interface.json)
./Allrun $C dsmc -np 4             # refuses an unsuitable interface unless post-cfd had --allow-unsuitable
./Allrun $C post-dsmc              # dsmc-post, reference, report -> results/summary.json, *.png
./Allrun $C crossings -np 4        # optional: exact particle crossings over a per-step window
```

Every stage resumes after an interruption: rerun the same command. Nothing is
deleted; to start over, remove the case directory. HPC batch scripts:
`python -m plumetools.hybridplume hpc <case> --scheduler pbs|slurm --np N`.

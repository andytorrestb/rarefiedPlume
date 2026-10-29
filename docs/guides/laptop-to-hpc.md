# From a laptop to a cluster

How work is split between a Windows laptop (native Python and WSL) and a Linux
cluster running PBS, using only the controls the repository already has. The
cluster specifics in this repository are the owner's; every site value below
is a placeholder for yours.

## What runs where

| Task | Native Windows Python | WSL / Linux shell | Needs OpenFOAM |
|---|---|---|---|
| `pytest -m "not needs_openfoam"` (unit + regression) | yes | yes | no |
| Load and validate any `case.yaml`, `study.yaml`, `hybrid.yaml` from Python | yes | yes | no |
| `generate_cases.py --dry-run` (cai2012, health, Kn extension), `generate_interfaces.py --dry-run` | yes, as `python script.py …` | yes | no |
| `validate_cases.py --inputs-only` (Kn extension), `Allscale.py inventory`/`smoke`, `assess.py` | yes | yes | no |
| `python -m plumetools.hybridplume setup …`, `python -m plumetools.nozzleplume describe <case>` | yes | yes | no |
| `./Allmesh --dict-only`, `./Allmesh --describe`, `./AllmeshCases --dict-only` | no (bash scripts) | yes | no |
| Meshing, solving, `Allpost` | no | yes | **yes** |
| Slice imagery (`render_slices.py`) | with a Windows ParaView's `pvpython` (not exercised here) | with Linux `pvpython` | no |
| `hpc record` (cost-model bookkeeping) | no (`fcntl` locks are Unix only) | yes | no |

The bash runners are `#!/bin/bash`; on Windows run them under WSL (or Git Bash
for the pure-Python ones). The repository's `.venv/` is a **Linux** virtual
environment for WSL; native Windows needs its own Python with the package
installed. Line endings are not a concern: `plumetools` writes `\n`
explicitly.

## Set up each environment

Each environment — the laptop's native Python, WSL, the cluster — needs its own
install, and each OpenFOAM installation its own library builds:

```bash
./Allsetup                    # pip install -e ".[test,plots,cases]", build libplumeDsmcBoundaryModels.so, verify
./Allsetup --python-only      # skip everything that needs OpenFOAM
./Allsetup --check            # verify only
./Allsetup --extras test,plots   # without CaseFoam, if it will not install
applications/thermophysicalModels/Allwmake   # libhybridPlumeThermo.so -- Allsetup does NOT build it
```

A library compiled against one OpenFOAM build is not portable to another: build
both on the cluster after sourcing its OpenFOAM, and again in WSL. Details and
the ParaView side: [Installation](../getting-started/installation.md).

## Prepare on the laptop, run on the cluster

1. **Generate and check inputs locally.** Each family's cheap checks
   (dry-run generation, `--dict-only`, `--describe`, validation) need no
   OpenFOAM. They prove the inputs are consistent and print the estimated
   cells, parcels, steps and memory — read those before asking for a node.
2. **Optionally run a smoke case locally**, where a family provides one: the
   Kn-extension interface `smoke_Kn0p001_B0p05`, the hybrid `smoke-*`
   profiles. They check mechanics, never physics. Never use all cores of the
   laptop: pass `-np` (or rely on the profile's `resources.max_ranks`), and
   leave memory headroom — the owner's WSL VM is 15 GB on a 32 GB host.
3. **Move the inputs, not the outputs.** On the cluster, update the checkout,
   rebuild the libraries, regenerate the cases there from the same study files
   (generation is deterministic and idempotent), and mesh there: a mesh made on
   a laptop is not guaranteed portable, and the HPC `case.yaml` variants
   (`r-4d/baseCaseHPC`) exist precisely because production meshes differ.
4. **Submit** with the family's PBS scripts ([PBS jobs](pbs.md)), or run the
   same `Allrun` commands interactively on a compute node.
5. **Bring back small files.** Copy `results/` directories, manifests and the
   evidence files (`interface-evidence.yaml`, `overlap.yaml`), not time
   directories, unless a figure needs fields; the families' post-processing
   runs where the fields are.

## Resource controls that exist

| Control | Where | Effect |
|---|---|---|
| `dsmc.n_subdomains` | `case.yaml` (all DSMC families) | written into `decomposeParDict`; the rank count `Allrun` uses |
| `-np N` | `Allrun`, study runners | overrides the rank count for one run via `system/decomposeParDict.Allrun` |
| `--serial` | `Allrun` | one process |
| `resources.max_ranks`, `resources.max_memory_GB` | `hybrid.yaml` profiles | the hybrid's rank default and memory budget |
| `mesh.max_cells` | Cai and nozzle-plume `case.yaml` | a hard cell budget; the mesh is coarsened to fit and every report says so |
| `resources.memory_ceiling_bytes` | `interfaces.yaml` | the interface assessment's memory limit |
| `resources.cores` (≤ 8), byte caps | `scaling.yaml` | bounds on the scaling fits |
| `OMP_NUM_THREADS` etc. | environment | clamped to 1–8 by `assess.py`, pinned to 1 by `Allscale.py` |

Open MPI counts **physical cores**, not hardware threads, when it decides how
many slots exist; that is why the Cai cases ask for 12 on a 12-core machine.

## Monitoring a run

| Look at | For |
|---|---|
| `log.dsmcFoam`, `Number of dsmc particles` | the parcel count levelling off — the only steady-state evidence the scripts give |
| `python -m plumetools.foamio.checkpoints report <case>` | the newest complete checkpoint on every rank, and any partial one |
| `python -m plumetools.foamio.checkpoints finished <case>` | exit 0 at `endTime` |
| `output.log`, `error.log` (Kn extension PBS) | the `NON-PRODUCTION` label, then `Allrun`'s output |
| `qstat` | the scheduler's view (site-specific) |

Restart after a walltime kill by resubmitting the same job; see
[Restart and recovery](restart-and-recovery.md).

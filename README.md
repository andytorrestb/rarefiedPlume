# rarefiedPlume

DSMC test cases and the `plumetools` package for rarefied plume expansion and
impingement, on OpenFOAM's standard `dsmcFoam` (verified against **v2512**; the
MNF `dsmcFoam+` fork is not required). Every active case is generated from one
configuration file; post-processing compares DSMC with analytical solutions
where they exist.

> **No result in this repository is physically validated.** The tests show the
> code does what it says; recorded runs show the pipelines work. Each family's
> tutorial states what is and is not established.

## Quick start

```bash
./Allsetup                              # pip install -e ".[test,plots,cases]", build the inflow library, verify
./Allsetup --python-only                # without OpenFOAM
python -m pytest -m "not needs_openfoam"   # the fast suite; no OpenFOAM needed
```

Then pick a family below. Before any run, its cheap checks (dry-run
generation, `--dict-only`, `--describe`) show what the inputs imply and cost
without OpenFOAM.

## Case families

| Family | What | Tutorial |
|---|---|---|
| [`cases/cai2012`](cases/cai2012/README.md) | Cai & Wang (2012) collisionless circular argon plume, Kn 100 / 0.1 / 0.01 × four particle levels | [docs/tutorials/cai2012.md](docs/tutorials/cai2012.md) |
| [`cases/cai2012-health`](cases/cai2012-health/README.md) | statistical noise of the Kn = 100 answer: parcels × averaging sweep | [docs/tutorials/cai2012-health.md](docs/tutorials/cai2012-health.md) |
| [`cases/cai2012-kn-extension`](cases/cai2012-kn-extension/README.md) | the same plume at seven decades of Kn (100 … 1e-4); interface track; field scaling | [direct](docs/tutorials/cai2012-kn-extension.md) · [interface](docs/tutorials/cai2012-interfaces.md) · [scaling](docs/tutorials/cai2012-scaling.md) |
| [`cases/cai2012-virtual-collision`](cases/cai2012-virtual-collision/README.md) | a reduced virtual-collision model against Cai and DSMC | [docs/tutorials/cai2012-virtual-collision.md](docs/tutorials/cai2012-virtual-collision.md) |
| [`cases/cai2012-hybrid`](cases/cai2012-hybrid/README.md) | one-way 3-D CFD → DSMC hybrid of Cai's source | [docs/tutorials/cai2012-hybrid.md](docs/tutorials/cai2012-hybrid.md) |
| [`cases/markelov1999`](cases/markelov1999/README.md) | AIAA 99-3455: N₂ plume past a cylinder onto a plate, four reservoir pressures | [docs/tutorials/markelov1999.md](docs/tutorials/markelov1999.md) |
| [`cases/r-4d`](cases/r-4d/README.md) | an R-4D exit-plane plume onto a plate, or free | [docs/tutorials/r-4d.md](docs/tutorials/r-4d.md) |
| [`cases/3d-inflow`](cases/3d-inflow/README.md) | the legacy source-flow reference case (frozen model; no mesh in a fresh clone) | [docs/tutorials/3d-inflow.md](docs/tutorials/3d-inflow.md) |
| `cases/iss_solar_panels{,_2,_3}`, `cases/solar_panel_particle_resolution_study` | unconverted legacy copies; need `dsmcFoam+`; **not runnable today** | [ISS](docs/tutorials/iss-solar-panels.md) · [particle resolution](docs/tutorials/particle-resolution.md) |

Archived, frozen historical record: `1d`, `2d-planar`, `2d-wedge`, `caseFoamEx`,
`wake-cylinder` — see [`cases/ARCHIVE.md`](cases/ARCHIVE.md). Which cases run
today, with what, and what blocks the rest:
[case catalog](docs/reference/case-catalog.md).

## What the active cases need

| | |
|---|---|
| OpenFOAM | **v2512**, standard distribution (`dsmcFoam`, `dsmcInitialise`, `blockMesh`, `snappyHexMesh`, …) |
| Custom code | `libplumeDsmcBoundaryModels.so` (`applications/dsmcBoundaryModels/Allwmake`, built by `Allsetup`); `libhybridPlumeThermo.so` for the hybrid only (`applications/thermophysicalModels/Allwmake`) |
| Python | ≥ 3.9 with numpy, scipy ≥ 1.6, PyYAML; extras `test`, `plots`, `cases` (CaseFoam, for two generators) |
| Imagery | ParaView's `pvpython` (not pip-installable); ffmpeg optional |

Historical results were produced with OpenFOAM v1706 and `dsmcFoam+`, whose
build is not recoverable: [docs/environment.md](docs/environment.md).

## Documentation

The documentation is a local MkDocs site; every page is also readable on
GitHub under [`docs/`](docs/index.md). To build it (Python ≥ 3.11, separate
environment, no simulation stack needed):

```bash
python3 -m venv .venv-docs
.venv-docs/bin/python -m pip install -r docs/requirements.txt
.venv-docs/bin/python -m mkdocs serve          # http://127.0.0.1:8000
.venv-docs/bin/python -m mkdocs build --strict # into site/ (gitignored)
```

Starting points: [installation](docs/getting-started/installation.md) ·
[a first session](docs/getting-started/first-run.md) ·
[repository map](docs/getting-started/repository-map.md) ·
[command line](docs/reference/cli.md) ·
[configuration files](docs/reference/configuration.md) ·
[outputs](docs/reference/outputs.md) ·
[models and evidence](docs/explanation/index.md) ·
[architecture](docs/plumetools-architecture.md) ·
[known gaps](docs/development/known-gaps.md).
The Python API reference is generated in the built site.

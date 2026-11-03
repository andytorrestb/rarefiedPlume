# Python API reference

Every module of `plumetools`, grouped by family. Each family page describes its
modules — purpose, principal interfaces, dependencies, the workflows that use
them — and gives tested usage examples; each module page renders the module's
own docstrings.

!!! info "How these pages are built"
    Module pages are generated at build time by `docs/_tools/hooks.py` from the
    files under `plumetools/`, and rendered by mkdocstrings from a **static**
    reading of the source: nothing is imported, so neither numpy nor ParaView is
    needed to build them. The reStructuredText in the docstrings (`.. math::`,
    `:func:` roles, simple tables) is translated to Markdown at render time by
    `docs/_tools/griffe_rst.py`; the source is unchanged.
    `docs/_tools/check_api.py` fails the build if a module is added or removed
    without its family page being updated. Private helpers (leading `_`) are
    not rendered unless a family page explains one.

!!! warning "No stability guarantee"
    `plumetools` is a research package at version 0.1.0. "Public" below means
    *meant to be called from outside the module*, not a promise that a
    signature will not change. The 2026-10-27 refactor kept every moved name
    importable from where it was ([record](../../refactor-2026-10-27.md)), which
    is the strongest compatibility statement the project makes.

## Families

| Family | Package | What it is | Status |
|---|---|---|---|
| [Legacy source flow](legacy/index.md) | top-level modules | the 3d-inflow source-flow model, its `case.yaml` loader and orchestration | frozen by the regression golden |
| [cai2012](cai2012/index.md) | `plumetools.cai2012` | Cai & Wang 2012 collisionless circular plume: config, gas, analytical solution, mesh, study, checks, post, health sweep, audit | active, output pinned byte for byte |
| [cai2012ext](cai2012ext/index.md) | `plumetools.cai2012ext` | the Kn extension: policy, generation, feasibility, comparisons, HPC wrappers, interface track, scaling | active |
| [markelov1999](markelov1999/index.md) | `plumetools.markelov1999` | AIAA 99-3455: paper-faithful source flow, geometry, mesh, flux verification, resolution, post | active, frozen package |
| [nozzleplume](nozzleplume/index.md) | `plumetools.nozzleplume` | general nozzle-exit plume with an optional plate (r-4d) | active |
| [hybridplume](hybridplume/index.md) | `plumetools.hybridplume` | one-way 3-D CFD → DSMC hybrid | active |
| [virtualcollision](virtualcollision/index.md) | `plumetools.virtualcollision` | reduced virtual-collision model and its assessment | active |
| [rarefaction](rarefaction/index.md) | `plumetools.rarefaction` | mean free path, Knudsen numbers and mesh resolution from a sampled run | shared tool |
| [viz](viz/index.md) | `plumetools.viz` | ParaView slice imagery, video, dissection | shared tool |
| [foamio](foamio/index.md) | `plumetools.foamio` | OpenFOAM field and dictionary I/O, checkpoints | shared layer |
| [mesh](mesh/index.md) | `plumetools.mesh` | `polyMesh` parsing and patch geometry | shared layer |
| [Shared workflow modules](shared/index.md) | `casemanifest`, `iteration` | manifests and fresh mesh iterations | shared tools |

How the layers depend on each other: [Architecture](../../plumetools-architecture.md).

## The package root

::: plumetools
    options:
      show_root_heading: false
      members: false

The root package defines `__version__` and nothing else a caller needs; import
from the family packages.

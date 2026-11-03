# rarefiedPlume

DSMC test cases and the `plumetools` package for rarefied plume expansion and
impingement, on OpenFOAM's standard `dsmcFoam` (verified against **v2512**).
A plume leaves a nozzle exit or an analytical source surface and expands into
vacuum, sometimes past a body or onto a plate; the repository generates every
case from one configuration file, runs it, and compares it with analytical
solutions where they exist.

!!! warning "Read the evidence before the results"
    No result in this repository is physically validated. The tests establish
    that the code does what it says; a handful of recorded runs show the
    pipelines work and how their answers behave. Each tutorial opens with what
    is and is not established for its family, and
    [What the tests establish](explanation/evidence.md) explains the levels.

## Start here

<div class="grid cards" markdown>

- **Install** — the package, OpenFOAM v2512, the two custom libraries, ParaView.
  [Installation](getting-started/installation.md)
- **First session** — the fast tests, what a case implies without running it,
  an optional small run. [A first session](getting-started/first-run.md)
- **Find your way** — layout, conventions, what is generated.
  [Repository map](getting-started/repository-map.md)
- **What runs today** — every case, its evidence and its blockers.
  [Case catalog](reference/case-catalog.md)

</div>

## Case families

| Family | Question it answers | Tutorial |
|---|---|---|
| `cai2012` | does 3-D `dsmcFoam` approach Cai's collisionless circular plume, and depart from it the right way as Kn falls? | [cai2012](tutorials/cai2012.md) |
| `cai2012-health` | how much statistical noise does that answer carry? | [cai2012-health](tutorials/cai2012-health.md) |
| `cai2012-kn-extension` | the same plume across seven decades of Kn — and which cases are even feasible | [direct cases](tutorials/cai2012-kn-extension.md) · [interface track](tutorials/cai2012-interfaces.md) · [field scaling](tutorials/cai2012-scaling.md) |
| `cai2012-virtual-collision` | can a reduced collision model correct the collisionless field? | [virtual collisions](tutorials/cai2012-virtual-collision.md) |
| `cai2012-hybrid` | a one-way CFD → DSMC hybrid in full 3-D | [hybrid](tutorials/cai2012-hybrid.md) |
| `markelov1999` | a plume past a cylinder onto a plate (AIAA 99-3455) | [markelov1999](tutorials/markelov1999.md) |
| `r-4d` | an R-4D exit plume, with or without a plate | [r-4d](tutorials/r-4d.md) |
| `3d-inflow` | the legacy source-flow reference case | [3d-inflow](tutorials/3d-inflow.md) |
| legacy copies | ISS solar panels, a particle-weight study — not runnable today | [ISS](tutorials/iss-solar-panels.md) · [particle resolution](tutorials/particle-resolution.md) |

Archived cases (`1d`, `2d-planar`, `2d-wedge`, `caseFoamEx`, `wake-cylinder`)
are historical record: [Archive](archive/index.md).

## The rest of the site

- **Guides** — [laptop to cluster](guides/laptop-to-hpc.md), [PBS jobs](guides/pbs.md),
  [restart and recovery](guides/restart-and-recovery.md),
  [post-processing](guides/postprocessing.md), [troubleshooting](guides/troubleshooting.md).
- **Reference** — [command line](reference/cli.md),
  [configuration files](reference/configuration.md), [outputs](reference/outputs.md),
  [Python API](reference/api/index.md).
- **Explanation** — the models, their assumptions, their code and their
  evidence: [Explanation](explanation/index.md).
- **Development** — [architecture](plumetools-architecture.md),
  [extending the workflows](development/extending-workflows.md),
  [writing these docs](development/documentation.md),
  [known gaps and defects](development/known-gaps.md).

This site is built locally (`python -m mkdocs build --strict`) and is not
published.

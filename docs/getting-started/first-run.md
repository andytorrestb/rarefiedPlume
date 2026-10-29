# A first session

A guided hour for someone new to the repository: check the install, look at
what a case's inputs imply, and — only if OpenFOAM is available — run a small
case end to end. Nothing here meshes an expensive case, submits a job, or
touches existing results. Commands start in the repository root.

## 1. The fast test suite

```bash
python -m pytest -m "not needs_openfoam" -q -p no:cacheprovider
```

No OpenFOAM is needed. Expect a few dozen skips (the regression golden without
its mesh, CaseFoam, Unix-only pieces); see
[What the tests establish](../explanation/evidence.md) for what a pass means
and does not mean.

## 2. What a case implies, without running it

Each family can print or check its derived inputs without OpenFOAM. None of the
following writes into the repository except where noted.

```bash
# the Cai matrix: 3 Knudsen numbers x 4 particle levels; lists the cases, writes nothing
python cases/cai2012/generate_cases.py --dry-run

# the R-4D plume: derived exit state, mesh plan, cost, resolution warnings; writes nothing
python -m plumetools.nozzleplume describe cases/r-4d/baseCase

# the Kn extension: the plan and the production verdict per Kn
# (stages in a temporary .staging-* directory beside study.yaml, removed on exit)
python cases/cai2012-kn-extension/generate_cases.py --dry-run
```

Read the feasibility table the last command prints: cell size over mean free
path, time step over collision time, estimated parcels — and which cases are
unsuitable for production, and why. That table is the first thing to look at
before any run in this repository.

## 3. The physics in a few lines of Python

<!-- docs-example: first-run-analytical -->
```python
from plumetools.cai2012 import analytical

# Cai's collisionless centreline for D = 0.2 m, S0 = 2: density, velocity, temperature
profile = analytical.centerline_profile([0.05, 0.2, 1.0], nozzle_radius_m=0.1, speed_ratio=2.0)
for x_over_D, n in zip(profile["x_over_D"], profile["n_over_n0"]):
    print(f"X/D = {x_over_D:4.2f}   n/n0 = {n:.4f}")
```

The derivation behind it is on the [cai2012 API page](../reference/api/cai2012/index.md#the-analytical-solution).

## 4. Set up a hybrid case (Python only, writes one directory)

```bash
python -m plumetools.hybridplume setup cases/cai2012-hybrid/runs/first \
    --config cases/cai2012-hybrid/hybrid.yaml --profile smoke-kn0p01-coarse --np 2
```

This writes `cases/cai2012-hybrid/runs/first/` (gitignored): the resolved
configuration, `results/setup.json` with resource estimates, and the CFD case
dictionaries. Delete the directory when done. It shows how a profile overlays
the base configuration and what a case will cost before anything runs.

## 5. Optional: a real run (OpenFOAM v2512 and the library)

Only with OpenFOAM sourced and `libplumeDsmcBoundaryModels.so` built
([Installation](installation.md)). The cheapest end-to-end DSMC case with a
documented size is the Kn-extension interface **development smoke** case,
recorded at about 37 000 cells and 2 × 10⁴ parcels, under two minutes per run
on the owner's laptop:

```bash
cd cases/cai2012-kn-extension
python generate_interfaces.py                   # writes Interface/ and interface-manifest.yaml (gitignored)
cd Interface/smoke_Kn0p001_B0p05
bash ./Allmesh
bash ./Allrun --serial
bash ./Allpost
```

It is labelled NON-PRODUCTION and checks mechanics only. The steps after it,
and what each output means, are in the
[interface-track tutorial](../tutorials/cai2012-interfaces.md#the-development-smoke-run-openfoam-minutes).

## Where to go next

- Pick the family you will work on in the [tutorials](../tutorials/index.md);
  each states its evidence, prerequisites and blockers first.
- The [case catalog](../reference/case-catalog.md) says which cases run today
  and what each needs.
- Before running on a cluster: [From a laptop to a cluster](../guides/laptop-to-hpc.md).

# rarefiedPlume

Test cases and utilities for a plume-impingement model for rarefied flows, built
on OpenFOAM's `dsmcFoam+`.

A rarefied plume expands from a hemispherical inflow surface into a vacuum. The
inflow conditions come from an analytical source-flow model evaluated per face;
DSMC solves the field from there.

## Quick start

```bash
pip install -e ".[test]"          # numpy, scipy, pyyaml (+ pytest)
pytest -m "not needs_openfoam"    # 124 tests, ~2 s, no OpenFOAM required

cd cases/3d-inflow
topoSet && dsmcInitialise+        # needs OpenFOAM v1706 + dsmcFoam+
python runInflow.py               # writes 0/boundaryU, boundaryT, boundaryNumberDensity_Ar
dsmcFoam+
```

All physical inputs live in [`cases/3d-inflow/case.yaml`](cases/3d-inflow/case.yaml).
No Python needs editing to change a parameter.

## Layout

```
plumetools/          the analytical model and OpenFOAM I/O, as a package
  sourceflow.py        equations E1-E8, pure functions
  geometry.py          centroids, normals, spherical coordinates
  mesh/                polyMesh parsing
  foamio/              0/ field writers and the blockMeshDict generator
  config.py            case.yaml loading and validation
  inflow.py            orchestration
cases/               standard OpenFOAM cases -- see the scope split below
docs/                workflow, model reconstruction, mesh, environment
tests/unit/          Tier 0: pure Python, no OpenFOAM
tests/regression/    bit-exact freeze of the pre-refactor behaviour
util/                one-off analysis scripts
```

Cases keep the ordinary OpenFOAM layout (`constant/`, `system/`, `0/`).
`plumetools` is a library they call, not a framework that replaces them.

## Active vs archived cases

| Status | Cases |
|---|---|
| **Active** | `3d-inflow` (reference) · `iss_solar_panels{,_2,_3}` · `solar_panel_particle_resolution_study/Fnum_*` |
| **Archived** | `1d` · `2d-planar` · `2d-wedge` · `caseFoamEx` · `wake-cylinder` |

Archived cases are one-off studies kept for historical record. They are frozen:
not refactored, not fixed, not deleted. Several produced results under labels
that do not describe what was actually run — [`cases/ARCHIVE.md`](cases/ARCHIVE.md)
records what is known to be wrong in each, so archived output is never mistaken
for validated output.

Only `cases/3d-inflow` has been migrated to `plumetools` so far. The other seven
active cases still carry their own copy of the legacy script.

## Read this before trusting the model

**The source-flow model has no recoverable citation.** An exhaustive search of
every comment, commit message and document found exactly one literature
reference, and it covers the *archived* 1d/2d-wedge validation case (Stewart &
Lumpkin 2012, via Jonathan S. Pitt, Aegis Aerospace, for LENS). The 3D model —
the limiting angle, the `(γ+0.41)/(γ−1)` angular exponent, the normalisation
integral — appears nowhere. It was reconstructed from the code.

Several known defects are **deliberately frozen** so the extraction could be
proven behaviour-preserving. Every one is a visible `LEGACY` value in `case.yaml`
rather than a buried literal. The most consequential:

- **SM-01** — the per-case stagnation pressure was never read. `calculateRhoN`
  hard-coded 475 psi while each case set its own value and then ignored it.
- **SM-02** — density used 300 K while velocity and temperature used the
  configured temperature; a factor of 10 apart in the resolution-study cases.
- **SM-03** — the model computes for N₂ while the solver simulates argon.
  v_ℓ = 788.16 m/s instead of 558.89 m/s.
- **SM-04** — θ is measured from +z, not from the +x plume axis.

Full reconstruction, confidence classification and open questions:
[`docs/source-flow-model.md`](docs/source-flow-model.md).

## Environment

| | |
|---|---|
| OpenFOAM | build **v1706** (recorded in `cases/wake-cylinder/5psi/sample_out.txt`) |
| Solver | `dsmcFoam+` / `dsmcInitialise+` — the MNF fork; version not recorded |
| Meshing | Pointwise **V18.5R2** (binary `.pw`; `cases/3d-inflow` no longer needs it) |
| Post-processing | ParaView **5.10.0** |
| Python | ≥ 3.9; numpy, scipy ≥ 1.6, pyyaml |

The Python library versions used to produce the historical results are **unknown
and unrecoverable**, so bit-level reproduction of archived output is not
achievable. See [`docs/environment.md`](docs/environment.md).

## Documentation

- [`docs/mesh.md`](docs/mesh.md) — geometry, the 5-block O-grid, and why
  generation had to wait for the centroid fix
- [`docs/source-flow-model.md`](docs/source-flow-model.md) — the equations, units,
  conventions, and what could not be verified
- [`cases/3d-inflow/README.md`](cases/3d-inflow/README.md) — the reference case,
  including two unresolved settings its own descendants later changed
- [`cases/ARCHIVE.md`](cases/ARCHIVE.md) — frozen cases and their known defects
- [`tests/regression/README.md`](tests/regression/README.md) — what the golden
  does and does not prove

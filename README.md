# rarefiedPlume

Test cases and utilities for a plume-impingement model for rarefied flows, built
on OpenFOAM's `dsmcFoam` (verified against **v2512**).

A rarefied plume expands from a hemispherical inflow surface into a vacuum. The
inflow conditions come from an analytical source-flow model evaluated per face;
DSMC solves the field from there.

## Quick start

```bash
./Allsetup                        # pip install, build the inflow library, verify
pytest -m "not needs_openfoam"    # 475 tests, ~4 s, no OpenFOAM required

cd cases/3d-inflow
./Allrun --no-solve               # generates 0/ inflow fields; needs no OpenFOAM
./Allrun                          # full run; needs OpenFOAM v2512
./Allclean                        # reset (--dry-run to preview)
```

`./Allsetup --python-only` skips everything needing OpenFOAM; `./Allsetup --check`
verifies an existing setup without installing or building anything. By hand it is:

```bash
pip install -e ".[test]"                           # numpy, scipy, pyyaml (+ pytest)
cd applications/dsmcBoundaryModels && ./Allwmake   # the custom inflow model
```

That second step is easy to miss and pip never does it — see
[Environment](#environment).

All physical inputs live in [`cases/3d-inflow/case.yaml`](cases/3d-inflow/case.yaml).
No Python needs editing to change a parameter.

## Layout

```
plumetools/          the analytical model and OpenFOAM I/O, as a package
  sourceflow.py        equations E1-E8, pure functions -- the LEGACY model
  geometry.py          centroids, normals, spherical coordinates
  mesh/                polyMesh parsing
  foamio/              0/ field writers, mesh dictionary generators
    primitives.py        searchable sphere/cylinder/box, shared by both families
  config.py            case.yaml loading and validation
  inflow.py            orchestration
  markelov1999/        the AIAA 99-3455 family: paper-faithful model, geometry,
                       mesh, flux verification, resolution, checks, post
applications/        custom OpenFOAM code
  dsmcBoundaryModels/  plumeFieldInflow -- per-face, patch-selected DSMC inflow
cases/               standard OpenFOAM cases -- see the scope split below
docs/                workflow, model reconstruction, mesh, environment
tests/unit/          Tier 0: pure Python, no OpenFOAM
tests/regression/    bit-exact freeze of the pre-refactor behaviour
tests/openfoam/      Tier 2: needs_openfoam, deselected by default
util/                one-off analysis scripts
```

Cases keep the ordinary OpenFOAM layout (`constant/`, `system/`, `0/`).
`plumetools` is a library they call, not a framework that replaces them.

## Active vs archived cases

| Status | Cases |
|---|---|
| **Active** | `markelov1999` (AIAA 99-3455) · `3d-inflow` (reference) · `iss_solar_panels{,_2,_3}` · `solar_panel_particle_resolution_study/Fnum_*` |
| **Archived** | `1d` · `2d-planar` · `2d-wedge` · `caseFoamEx` · `wake-cylinder` |

Archived cases are one-off studies kept for historical record. They are frozen:
not refactored, not fixed, not deleted. Several produced results under labels
that do not describe what was actually run — [`cases/ARCHIVE.md`](cases/ARCHIVE.md)
records what is known to be wrong in each, so archived output is never mistaken
for validated output.

`cases/markelov1999` and `cases/3d-inflow` use `plumetools`. The other seven
active cases still carry their own copy of the legacy script.

## The AIAA 99-3455 case family

`cases/markelov1999` reproduces the configuration of Lumpkin, Stewart & Markelov,
*Study of 3D Rarefied Flow on a Flat Plate in the Wake of a Cylinder* (1999): a
rarefied N₂ plume past a finite cylinder onto a flat plate 6 inches behind it, at
four reservoir pressures.

```bash
./Allsetup --extras test,cases                     # + CaseFoam, for case generation
cd cases/markelov1999
./generate_cases.py && ./AllmeshCases && ./AllrunCases && ./AllpostCases
```

Case generation uses [CaseFoam](https://github.com/DLR-RY/caseFOAM) to clone the
base case into the hierarchy, as `util/caseFoam/` has since before the refactor.
The physical values are applied afterwards by structured YAML editing rather than
through CaseFoam's `caseData`, because that mechanism reaches a non-OpenFOAM file
only through `'#!stringManipulation'` — whitespace-sensitive substitution.

It is **new work**, not a fix of the archived `wake-cylinder` lineage, and it
does not inherit that lineage's defects. Its source-flow model is a separate,
explicitly named variant — `markelov1999_axisymmetric` — because
`plumetools/sourceflow.py` is frozen to reproduce those defects and disagrees
with the paper in six substantive places (SF-1 … SF-6 in
[`docs/markelov1999-case.md`](docs/markelov1999-case.md)).

**No result in it has been validated** against DAC or experiment.

It also carries the repository's first custom OpenFOAM code:
[`applications/dsmcBoundaryModels/plumeFieldInflow`](applications/dsmcBoundaryModels/plumeFieldInflow),
an `InflowBoundaryModel` that injects across a named patch list with a per-face
number density — the two things stock `FreeStream` cannot do, and the reason
`3d-inflow` has to reflect at its outer boundary.

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

### What the active cases need today

| | |
|---|---|
| OpenFOAM | **v2512**, the standard distribution — verified by reading its source and running the cases |
| Solver | `dsmcFoam` / `dsmcInitialise`. **Not** the MNF fork's `dsmcFoam+` |
| Custom code | `libplumeDsmcBoundaryModels.so`, built by `applications/dsmcBoundaryModels/Allwmake` |
| Python | ≥ 3.9; numpy, scipy ≥ 1.6, pyyaml |
| Post-processing | ParaView **5.10.0** (`paraview.simple`; not pip-installable) |

The fork is not required because the one thing it provided that this work needs —
a per-face, patch-selected inflow — is supplied instead by
[`plumeFieldInflow`](applications/dsmcBoundaryModels/plumeFieldInflow), a peer of
OpenFOAM's own `FreeStream` that links against the standard distribution. Run
`./Allsetup` to install and build all of it, or `./Allsetup --check` to find out
which piece is missing. Details in
[`docs/solver-compatibility.md`](docs/solver-compatibility.md).

### What produced the historical results

| | |
|---|---|
| OpenFOAM | build **v1706** (recorded in `cases/wake-cylinder/5psi/sample_out.txt`) |
| Solver | `dsmcFoam+` / `dsmcInitialise+` — the MNF fork; version not recorded |
| Meshing | Pointwise **V18.5R2** (binary `.pw`; `cases/3d-inflow` no longer needs it) |

This row is provenance, not an instruction: the archived cases still name
`dsmcFoam+` in their `controlDict`s and are frozen that way
([`cases/ARCHIVE.md`](cases/ARCHIVE.md)). The Python library versions behind
those results are **unknown and unrecoverable**, so bit-level reproduction of
archived output is not achievable. See
[`docs/environment.md`](docs/environment.md).

## Documentation

- [`docs/markelov1999-case.md`](docs/markelov1999-case.md) — the AIAA 99-3455
  family: the model, the geometry, and every assumption with its reasoning
- [`docs/plume-field-inflow.md`](docs/plume-field-inflow.md) — the custom DSMC
  inflow boundary model, and why standard `dsmcFoam` needs one
- [`cases/markelov1999/README.md`](cases/markelov1999/README.md) — running the
  case family
- [`docs/mesh.md`](docs/mesh.md) — geometry, the 5-block O-grid, and why
  generation had to wait for the centroid fix
- [`docs/source-flow-model.md`](docs/source-flow-model.md) — the equations, units,
  conventions, and what could not be verified
- [`cases/3d-inflow/README.md`](cases/3d-inflow/README.md) — the reference case,
  including two unresolved settings its own descendants later changed
- [`cases/ARCHIVE.md`](cases/ARCHIVE.md) — frozen cases and their known defects
- [`tests/regression/README.md`](tests/regression/README.md) — what the golden
  does and does not prove

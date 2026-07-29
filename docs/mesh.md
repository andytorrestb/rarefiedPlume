# Mesh generation for `cases/3d-inflow`

## Why this exists

Every mesh in this repository was exported from **Pointwise V18.5R2**. The `.pw`
project files under `cases/*/mesh/` are **binary** — `file` reports `data`, and the
header is `PWI0`. They cannot be diffed, scripted, or regenerated without an
interactive, licensed Pointwise session. The committed `constant/polyMesh` is
therefore the only usable artifact, and step 1 of the workflow was not
reproducible at all (finding RP-16).

Worse, `cases/3d-inflow/system/blockMeshDict` was **the archived 2d-wedge
dictionary**, copied in by a directory copy and never replaced. It describes a 2D
axisymmetric wedge — `front`/`back` patches of type `wedge`, an `axis` patch of
type `empty`, one block `hex (0 1 2 0 5 3 4 5) (132 1 164)`, inflow radius 1.05 m
expanding to 20 m. Running `blockMesh` — the first thing anyone attempting
reproduction would try — would have **silently overwritten `constant/polyMesh`**
and replaced it with an unrelated 2D wedge (finding FD-08).

`plumetools/foamio/blockmesh.py` generates a correct dictionary from `case.yaml`,
making the geometry text, diffable, parametric and version-controlled.

## The geometry

Measured from the committed mesh, not assumed:

| | |
|---|---|
| Domain | box, x ∈ [0, 5] m, y, z ∈ [−2.5, 2.5] m |
| `inflow` | hemisphere, radius **exactly 0.5000000000 m**, centred at the origin, x ≥ 0 |
| `sym` | the x = 0 plane: a 5 × 5 square minus the r = 0.5 circle |
| `vacuum` | the five outer box faces (x = 5, y = ±2.5, z = ±2.5) |
| Committed cells | 99 240 tetrahedra (204 717 faces, 192 243 internal) |
| Committed inflow | 2044 triangles |

```
        z
        ^          vacuum (5 outer box faces)
        |   +---------------------------+
        |   |                           |
     ---+---|-- ) inflow                |----> x
        |   |  (hemisphere, R = 0.5 m)  |
        |   +---------------------------+
            ^ sym (the x = 0 plane)
```

The plume axis is **+x**. That is inferred from the one-sided hemisphere, from
`revolutionAxis "x"` in `dsmcProperties`, and from the first sampling line — it is
stated nowhere in the original code.

## Topology — a 5-block O-grid

Take the +x half of a cube inscribed in the sphere and connect each of its five
outward faces to the corresponding box face. With `R` the sphere radius, `H` the
box half-width and `L` the box length:

* `a = R/√3 = 0.2886751346` — a cube corner projected onto the sphere
* `b = R/√2 = 0.3535533906` — an equator point at 45°

**16 vertices**

| Index | Position | Role |
|---|---|---|
| 0–3 | `(a, ±a, ±a)` | inner cap corners, on the sphere |
| 4–7 | `(0, ±b, ±b)` | inner equator ring, sphere ∩ x = 0 |
| 8–11 | `(L, ±H, ±H)` | outer, x = L |
| 12–15 | `(0, ±H, ±H)` | outer, x = 0 |

**5 blocks**

| Block | Inner face (sphere) | Outer face (box) | Also contributes |
|---|---|---|---|
| cap | the four `(a, ±a, ±a)` | the four `(L, ±H, ±H)` | — |
| +y | `(a,a,±a)`, `(0,b,±b)` | y = +H | `sym` |
| −y | `(a,−a,±a)`, `(0,−b,±b)` | y = −H | `sym` |
| +z | `(a,±a,a)`, `(0,±b,b)` | z = +H | `sym` |
| −z | `(a,±a,−a)`, `(0,±b,−b)` | z = −H | `sym` |

The cap block does not touch x = 0; each side block contributes one x = 0 face to
`sym`. The five inner faces tile the hemisphere exactly — the +x half of a cube's
six faces is the +x face whole, half of each of ±y and ±z, and no −x face.

**12 curved edges** — 4 cap, 4 meridian (cap corner → equator point in the same
octant), 4 equator. The `arc` interpolation point is computed analytically:

```
P_mid = R * (P1 + P2) / |P1 + P2|
```

so the whole dictionary is a pure function of `(R, H, L, n_tangential, n_radial,
radial_grading)`.

**Orientation is derived, not hand-written.** Blocks are flipped if their
inner-to-outer normal points the wrong way, and each boundary face is wound
outward from the domain. All five blocks turn out to need flipping relative to
the naive vertex order — getting that wrong yields negative-volume cells, and it
is not obvious by inspection. Boundary faces are classified *geometrically* (on
the sphere → `inflow`; at x = 0 → `sym`; on a box face → `vacuum`) rather than by
index lists, so a change to the block table cannot silently mis-assign a patch.

## Resolution

Inflow faces = `5 · n_tangential²`; cells = `5 · n_tangential² · n_radial`.

The default `n_tangential: 20` gives **2000 inflow quads**, close to the Pointwise
mesh's 2044 triangles, which keeps the two comparable. `n_radial: 24` with
`radial_grading: 10.0` clusters cells toward the sphere, giving 48 000 cells.

**This is a starting resolution, not a converged one.** Mesh-convergence
requirements are explicitly out of scope here.

## ⚠ blockMesh gives hexahedra — the centroid consequence

The committed Pointwise mesh is **tetrahedral** with **triangular** inflow faces.
`blockMesh` produces **hexahedra** with **quadrilateral** inflow faces.

The pre-refactor `calculateCentroid` divided the vertex sum by the literal `3.0`
regardless of vertex count. On a quad that gives exactly **4/3** of the true
centroid. This is not hypothetical: it is what silently corrupted the archived
wake-cylinder cases, whose inflow faces are quads — measured code centroid mean
|r| **0.3946 m** against a true **0.2960 m**, ratio 1.3331 (finding AR-01).

`plumetools.geometry.centroid` divides by `len(vertices)`, which is
arithmetically identical on triangles (so the golden is untouched) and correct on
quads. **Enabling mesh generation before that fix would have reproduced AR-01 in
the reference case.** Do not reintroduce the `3.0`.

Two other consequences:

* A generated mesh emits `type patch;`, permanently sidestepping **HA-07** (the
  committed mesh declares `type Unspecified;`, a Pointwise export artifact that
  is not a registered OpenFOAM patch type).
* `mesh.sphere_radius_m` and `geometry.sphere_radius_m` must agree —
  `plumetools.config` enforces it. The model uses the radius as a *fixed* value,
  so a mismatch would evaluate the source flow on a different sphere than the one
  meshed (finding SM-09).

## Curvature: `arc` versus `project`

The generator uses **`arc` edges only**. `arc v0 v1 (px py pz)` is stable across
every OpenFOAM version, and the interpolation points are exact.

Its limitation: with arc edges alone, block *faces* on the sphere are ruled
surfaces, so interior surface points bulge very slightly inside the true sphere
(second order in the block angle). Projecting the faces onto a
`searchableSphere` would place every surface point exactly on the sphere.

`mesh.projection: searchable_sphere` is accepted by the config schema but raises
`NotImplementedError`. **The exact `project` directive syntax could not be
verified against OpenFOAM v1706 in this environment**, and shipping an unverified
directive that `blockMesh` might reject — or worse, silently ignore — would be
less useful than not shipping it. Enabling it needs a working DSMC environment
and a `needs_openfoam` test.

## Running it

```bash
cd cases/3d-inflow
./Allmesh --yes
```

`Allmesh` writes `system/blockMeshDict`, runs `blockMesh` and `checkMesh`, then
runs `python -m plumetools.verify_mesh .`. It requires an explicit `--yes`, and
`Allrun` never calls it.

**It overwrites `constant/polyMesh`.** Restore the Pointwise mesh with:

```bash
git checkout -- constant/polyMesh
```

## Verification

`checkMesh` covers mesh *quality*. `plumetools.verify_mesh` covers whether the
mesh is the geometry the model is about to be evaluated on — which `checkMesh`
cannot know:

* every inflow vertex lies on the sphere at `geometry.sphere_radius_m`
* inflow face count is `5 · n_tangential²`, and all faces are quads
* all inflow centroids have x ≥ 0 (a hemisphere, not a sphere)
* centroids sit inside the sphere, as faceting requires
* inflow face normals are consistently oriented
* no patch has `type Unspecified`

## Why `constant/polyMesh` stays tracked

Even though it is now generable, the committed `polyMesh` remains under version
control, and the `.gitignore` rules that were supposed to exclude it have been
removed rather than repaired. Two reasons:

1. It is the **historical Pointwise mesh** — the one every committed result was
   produced on. A generated mesh is an equivalent geometry, not the same mesh.
2. It is the substrate of the regression golden.

(Those ignore rules never worked anyway: `constant/polyMesh/*` contains a slash,
so it anchored to the repository root and matched nothing under `cases/*/`.
Confirm with `git check-ignore -v cases/3d-inflow/constant/polyMesh/points`,
which exits 1 — finding RP-15.)

The generated `system/blockMeshDict` **is** gitignored: it is derived from
`case.yaml`, which is the file under review.

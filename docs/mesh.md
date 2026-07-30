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

## Curvature: `arc` versus `project` — and why `arc` is not enough

`mesh.projection` selects how the inflow surface is made curved.

### `arc` — exact edges, **wrong faces**

`arc v0 v1 (px py pz)` curves the twelve block edges, with interpolation points
computed exactly as `R·(P₁+P₂)/|P₁+P₂|`.

But blockMesh fills a face *interior* by transfinite interpolation of its four
edges. Curving the edges therefore does not curve the face. At the face centre
that interpolation reduces to

```
centre = Σ(edge midpoints)/2 − Σ(corners)/4
```

which lands **inside** the sphere:

| Face | \|r\| at centre | Deficit |
|---|---|---|
| cap (spans 109.5°) | 0.418432 | **16.3 %** |
| each side face | 0.479287 | 4.1 % |

On a 0.5 m sphere that is **82 mm** of missing radius at the cap centre. The
inflow patch is not a hemisphere; it is five inward-dished panels joined at exact
edges.

**Refining does not help.** The face interior is determined by the four edges
alone, so increasing `n_tangential` puts *more* points on the same wrong surface.
This is a topological error, not a resolution error — which is why it survives
any amount of mesh refinement and why `checkMesh` is perfectly happy with it.

### `searchable_sphere` — exact surface *(default)*

Declares the sphere as an analytic primitive and projects both the edges **and**
the five inflow faces onto it:

```
geometry
{
    inflowSphere
    {
        type    searchableSphere;
        centre  (0 0 0);
        radius  0.5;
    }
}

edges ( project 0 1 (inflowSphere) ... );        // 12
faces ( project (0 1 2 3) inflowSphere ... );    //  5
```

Every surface point then lies on the sphere, at any resolution. Same topology,
same cell count, same face count — it costs nothing.

Projecting the edges *without* the faces would still leave the dished interiors,
so the `faces` section is the part that actually matters. A unit test asserts
both sections are present and that the projected quads are exactly the five the
`inflow` patch declares.

> **Version note.** v1706 spells the sphere centre `centre`; v2006 and later also
> accept `origin`. If your build rejects `centre`, that is the key to change.
> The projection syntax has **not** been executed against a real blockMesh here —
> only its structure is tested. `arc` remains available as a fallback.

`verify_mesh` measures the face-centroid deficit and reports it, so the
difference is visible without inspecting the mesh by eye.

## Alternative: `mesh.type: snappy_hex_sphere`

A second generator, using snappyHexMesh with the same `searchableSphere`
primitive instead of an O-grid. Set `mesh.type: snappy_hex_sphere` and `Allmesh`
writes three dictionaries and runs a two-stage pipeline:

```
system/blockMeshDict       a plain graded box -- the background mesh
system/snappyHexMeshDict   the searchableSphere, castellate and snap controls
system/meshQualityDict     quality limits, #included by the above

blockMesh                  # background box
snappyHexMesh -overwrite   # carve the cavity, snap to the sphere
```

The hemisphere falls out of the geometry rather than being constructed: the
sphere is centred on the `x = 0` face of the box, so only its `+x` half
intersects the mesh. `locationInMesh` is placed in the fluid, off every symmetry
plane, so snappyHexMesh keeps the region connected to it and discards the cavity
interior. The geometry entry is *named* `inflow`, because snappyHexMesh names the
patch it creates after the surface. `patchInfo` sets `type patch` — never `wall`,
which would make the plume source a no-slip surface.

Controls: `background_cell_size_m` (≤ R/4, else the cavity falls between cells),
`refinement_level`, `n_cells_between_levels`.

**Prefer the projected O-grid unless something rules it out.** snappyHexMesh
gives up several things the O-grid has:

| | projected O-grid | snappyHexMesh |
|---|---|---|
| Cells | pure hexahedra | polyhedra near the surface |
| Inflow faces | quads, exactly `5·n²` | polygons, count set by refinement |
| Surface accuracy | exact | within the snapping tolerance |
| Radial grading | structured, controllable | octree levels only |
| Comparability with the 2044-face Pointwise mesh | direct | not like-for-like |

Two consequences worth knowing. First, snapped faces have more than four
vertices, so the legacy `/3.0` centroid divisor would be wrong on *every one of
them* — `plumetools.geometry.centroid` divides by `len(vertices)`, which is what
makes this mesh type usable at all (finding AR-01). Second, the sphere meets
`x = 0` tangentially, and snappyHexMesh is least reliable where a surface grazes
a domain boundary; inspect the `sym`/`inflow` intersection before trusting it.

## Running it

```bash
cd cases/3d-inflow
./Allmesh --dict-only     # write the dictionaries, no OpenFOAM needed
./Allmesh --yes           # generate, mesh, checkMesh, verify
```

`Allmesh` dispatches on `mesh.type`, runs the matching pipeline, then runs
`python -m plumetools.verify_mesh .`. It requires an explicit `--yes`, and
`Allrun` never calls it.

**It overwrites `constant/polyMesh`.** Restore the Pointwise mesh with:

```bash
git checkout -- constant/polyMesh
```

## Verification

`checkMesh` covers mesh *quality*. `plumetools.verify_mesh` covers whether the
mesh is the geometry the model is about to be evaluated on — which `checkMesh`
cannot know:

* every inflow vertex lies on the sphere at `geometry.sphere_radius_m` — exactly
  for an O-grid, within 2 % for a snapped mesh
* inflow face count is `5 · n_tangential²` and all faces are quads (O-grid only;
  a snapped mesh reports its polygon mix instead)
* all inflow centroids have x ≥ 0 (a hemisphere, not a sphere)
* the **face-centroid deficit**, which is what distinguishes a projected surface
  from a dished `arc` one — the vertex check alone cannot see it
* inflow face normals are consistently oriented

`type Unspecified` is reported as a **warning**, not a failure. It is a Pointwise
export artifact (HA-07) and a property of a mesh that already exists; failing on
it would make the tool useless on the committed mesh, which is exactly where
someone would first run it. Generated meshes never emit it, and a unit test
enforces that.

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

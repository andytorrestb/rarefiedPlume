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

`plumetools` generates the dictionaries from `case.yaml`, making the geometry
text, diffable, parametric and version-controlled.

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

## The pipeline — `mesh.type: snappy_hex_sphere` *(what the case uses)*

**`blockMesh` builds the background box. `snappyHexMesh` carves the inflow
cavity.** The hemisphere appears in no `blockMeshDict`.

```
system/blockMeshDict       a plain uniform box — the background mesh
system/snappyHexMeshDict   the searchableSphere, castellate and snap controls
system/meshQualityDict     quality limits, #included by the above

blockMesh                  # background box: 20 × 20 × 20 = 8000 cells of 0.25 m
snappyHexMesh -overwrite   # carve the cavity, refine 4 levels, snap to the sphere
```

The background `blockMeshDict` is genuinely plain: eight vertices, one
`hex (0 1 2 3 4 5 6 7) (20 20 20) simpleGrading (1 1 1)`, an empty `edges` block,
and two patches — `sym` on `x = 0` and `vacuum` on the other five faces. There is
no `geometry` section and no `inflow` patch; snappyHexMesh creates that patch when
it carves.

### Why the hemisphere moved out of `blockMeshDict`

Building it there requires the 5-block O-grid below, and that topology has a
quality floor that refinement cannot lift. Five blocks meet at the cube corners
projected onto the sphere, and `radial_grading` shears every cell in the graded
direction — both are properties of the block topology, so raising `n_tangential`
puts more cells on the same badly conditioned arrangement.

snappyHexMesh keeps the background hexes axis-aligned and orthogonal everywhere
except the layer it snaps, so the bulk of the mesh is a perfect Cartesian grid and
non-orthogonality is confined to the cells touching the sphere.

**Both meshed and measured** with `checkMesh` — O-grid at `n_tangential: 20` /
`n_radial: 24` / `radial_grading: 10`, snappy at `0.25 m` / level 4:

| | O-grid | snappyHexMesh |
|---|---|---|
| Max non-orthogonality | 66.1 (limit 70) | **36.9** |
| Mean non-orthogonality | 31.2 | **11.9** |
| Max skewness | 1.79 | **0.76** |
| Max aspect ratio | 4.05 | 3.74 |
| Cells | 48 000 | **32 272** |
| Inflow faces | 2000 quads | 5452 quads/pentagons/hexagons |
| Total volume | 124.738729 m³ | 124.738350 m³ |

The mean is the number that matters most: 31.2 → 11.9 is the graded O-grid's
shear disappearing. Cell mix is 25 368 hexahedra, 792 prisms and 6112 polyhedra —
the polyhedra are the snapped layer.

### The cavity falls out of the geometry

The sphere is centred at the origin, which lies **on** the `x = 0` face of the box
(`x ∈ [0, L]`). Only its `+x` half intersects the mesh, so that half is what gets
carved — the hemisphere is not constructed, it is what remains.

The region inside the cavity is bounded by the sphere on the `+x` side and by the
`sym` patch at `x = 0`, which makes it disconnected from the fluid.
`locationInMesh` sits in the fluid, so snappyHexMesh keeps that region and
discards the cavity interior.

The geometry entry is *named* `inflow`, because snappyHexMesh names the patch it
creates after the surface. `patchInfo` sets `type patch` — never `wall`, which
would make the plume source a no-slip surface.

### Resolution

Only `background_cell_size_m` sizes the box; the octree sets the surface:

```
surface cell size = background_cell_size_m / 2**refinement_level
                  = 0.25 / 2**4 = 0.015625 m
```

so a **coarse background with more levels** is what replaces the O-grid's
`radial_grading` — cells go where the plume is dense, and the near-vacuum far
field stays cheap. At the case's settings that is 8000 background cells refining
to **32 272 total, below the O-grid's 48 000**, with a surface 2.5× finer.
snappyHexMesh takes about 3 s to do it.

The background only has to *find* the sphere, not resolve it, so
`background_cell_size_m` may be as large as `R` (two cells across the diameter).
Above that the generator refuses, because the cavity would fall between cells.

> **Do not raise `refinement_level` past 4 without lowering `controlDict`'s
> `deltaT`.** At `1e-5` s and ~800 m/s a particle already crosses half a
> 0.015625 m cell per step; one more level puts it over a full cell, which DSMC
> does not tolerate.

`n_cells_between_levels: 2` sets the buffer between octree levels. Fewer means a
cheaper but more abrupt transition.

### `locationInMesh` — a cell centre, not a fraction

A `locationInMesh` on a face, an edge, or a plane of symmetry is a classic way to
make snappyHexMesh keep the wrong region or fail outright, and *off-axis is not
sufficient on its own*: any fixed fraction of the domain lands exactly on a cell
face whenever it happens to be a multiple of the background cell size. The
previous `0.5 * L` did, for every even `nx` — including the default 40.

The seed is therefore snapped to the **centre of the background cell containing
it**. A cell centre is interior by construction at any cell size, so the point
tracks the grid instead of drifting onto it. A unit test sweeps six cell sizes and
asserts the seed sits at exactly 0.5 of the way through its cell on every axis.

### The rim, and what to inspect

The sphere is *bisected* by `x = 0`, so it crosses that plane at right angles —
the well-conditioned case, not a tangency. But the rim where `inflow` meets `sym`
is still a surface / domain-boundary intersection, which is where snapping is
least predictable. Feature snapping is deliberately **off** (a sphere has no
feature edges, and enabling it tends to pull rim points off the surface rather
than onto it), so that circle is the first thing to look at in ParaView.

In the meshed result it behaves: the lowest inflow face centroid sits at
x = 0.0067 m, so the patch closes cleanly onto `sym` rather than wrapping past it,
and `sym` comes out as 1420 faces — a 5 × 5 square minus the r = 0.5 circle, as it
should be. snappyHexMesh does emit a handful of `Displacement ... points through
the surrounding patch faces` warnings while snapping the rim, then reports
`Finished meshing without any errors`.

Snapping to a **primitive** rather than an STL turns out to cost nothing in
accuracy: the worst inflow-vertex deviation from r = 0.5 m is **8.3e-16 m**, which
is round-off, not tolerance. The remaining 0.04 % face-centroid deficit
(0.499803 … 0.499919) is ordinary faceting.

## Alternative: `mesh.type: block_mesh_ogrid`

The original generator, and no longer what this case uses. It builds the
hemisphere directly in `blockMeshDict` as a 5-block O-grid, giving an exact
sphere, pure hexahedra, and a face count of exactly `5·n²` — at the quality cost
tabulated above.

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

### Resolution (O-grid)

Inflow faces = `5 · n_tangential²`; cells = `5 · n_tangential² · n_radial`.

The default `n_tangential: 20` gives **2000 inflow quads**, close to the Pointwise
mesh's 2044 triangles, which keeps the two comparable. `n_radial: 24` with
`radial_grading: 10.0` clusters cells toward the sphere, giving 48 000 cells.

**This is a starting resolution, not a converged one.** Mesh-convergence
requirements are explicitly out of scope here.

## ⚠ A generated mesh is not tetrahedral — the centroid consequence

The committed Pointwise mesh is **tetrahedral** with **triangular** inflow faces.
`blockMesh` produces **hexahedra** with **quadrilateral** inflow faces, and
snappyHexMesh produces **polyhedra** with **polygonal** ones — faces with more
than four vertices are the norm, not the exception, on a snapped surface.

The pre-refactor `calculateCentroid` divided the vertex sum by the literal `3.0`
regardless of vertex count. On a quad that gives exactly **4/3** of the true
centroid. This is not hypothetical: it is what silently corrupted the archived
wake-cylinder cases, whose inflow faces are quads — measured code centroid mean
|r| **0.3946 m** against a true **0.2960 m**, ratio 1.3331 (finding AR-01).

`plumetools.geometry.centroid` divides by `len(vertices)`, which is
arithmetically identical on triangles (so the golden is untouched) and correct on
quads and polygons alike. **Enabling mesh generation before that fix would have
reproduced AR-01 in the reference case**, and the snapped mesh would have hit it
on nearly every face. Do not reintroduce the `3.0`.

Two other consequences:

* A generated mesh emits `type patch;`, permanently sidestepping **HA-07** (the
  committed mesh declares `type Unspecified;`, a Pointwise export artifact that
  is not a registered OpenFOAM patch type).
* `mesh.sphere_radius_m` and `geometry.sphere_radius_m` must agree —
  `plumetools.config` enforces it. The model uses the radius as a *fixed* value,
  so a mismatch would evaluate the source flow on a different sphere than the one
  meshed (finding SM-09).

## O-grid curvature: `arc` versus `project` — and why `arc` is not enough

`mesh.projection` selects how the **O-grid's** inflow surface is made curved. It
does nothing under `snappy_hex_sphere`, which snaps to the primitive directly.

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

### `searchable_sphere` — exact surface *(the O-grid default)*

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

**Verified against a real blockMesh** (2026-08-01), at the default
`n_tangential: 20` / `n_radial: 24`:

| | `arc` | `searchable_sphere` |
|---|---|---|
| Worst inflow-vertex deviation from r = 0.5 m | 0.0816 m | **9.51e-12 m** |
| Worst face-centroid deficit | 16.31 % of R | **0.16 % of R** |
| Inflow faces | 2000 quads | 2000 quads |
| Mesh | 51 025 points, 146 960 faces | unchanged |

The residual 0.16 % is ordinary faceting — a flat quad's centroid necessarily
sits inside the sphere its corners lie on — and it shrinks with `n_tangential`,
unlike the 16.31 % which did not.

> **Version note.** The `centre` key works on the build used here. v2006 and
> later also accept `origin`; if a build rejects `centre`, that is the key to
> change. `arc` remains available as a fallback.

`verify_mesh` measures the face-centroid deficit and reports it, so the
difference is visible without inspecting the mesh by eye.

## Choosing between the two

| | snappyHexMesh *(current)* | projected O-grid |
|---|---|---|
| Max / mean non-orthogonality | **36.9 / 11.9** | 66.1 / 31.2 |
| Cells | 32 272; polyhedra near the surface, hexes elsewhere | 48 000 pure hexahedra |
| Inflow faces | 5452 polygons, count set by refinement | 2000 quads, exactly `5·n²` |
| Surface accuracy | 8.3e-16 m from the primitive | exact |
| Near-source clustering | octree levels | structured `radial_grading` |
| Comparability with the 2044-face Pointwise mesh | not like-for-like | not like-for-like |

The O-grid wins on face-count predictability; snappyHexMesh wins on cell quality
and cell count, which is why this case uses it — and because the surface is a
primitive rather than an STL, it gives up essentially nothing on accuracy. Switch
back by setting
`mesh.type: block_mesh_ogrid` in `case.yaml` and re-running `./Allmesh --yes` —
nothing else in the case needs to change, and `Allmesh` picks up the different
pipeline on its own.

## Running it

```bash
cd cases/3d-inflow
./Allmesh --dict-only     # write the dictionaries, no OpenFOAM needed
./Allmesh --yes           # generate, mesh, checkMesh, verify
```

`Allmesh` dispatches on `mesh.type`, runs the matching pipeline, then runs
`checkMesh` and `python -m plumetools.verify_mesh .` and prints the quality
numbers. It requires an explicit `--yes`, and `Allrun` never calls it.

**It overwrites `constant/polyMesh`.** Restore the Pointwise mesh with:

```bash
git checkout -- constant/polyMesh
```

### What `Allmesh` clears first

`blockMesh` rewrites `points`, `faces`, `owner`, `neighbour` and `boundary` —
and nothing else. Everything else in `constant/polyMesh` addresses the *previous*
mesh by index, so `Allmesh` deletes it before meshing:

| Removed | Why it would break the next mesh |
|---|---|
| `sets/` | Face labels into the old mesh. `runInflow.py` **prefers** `sets/<patch>` over the patch's face range, and only a count mismatch is caught — equal counts would apply the inflow model to the wrong faces. |
| `cellLevel`, `pointLevel`, `level0Edge`, `refinementHistory`, `surfaceIndex` | snappyHexMesh's own octree bookkeeping from a previous `-overwrite` run, sized for that mesh. |
| `cellZones`, `faceZones`, `pointZones` | Pointwise export leftovers, addressed by index. |

All of them are regenerated by the pipeline or by `topoSet`, so removing them
costs nothing. Because `sets/inflow` is deleted, **`./Allrun` must be re-run after
`./Allmesh`** — it runs `topoSet` and then `runInflow.py`, in that order.

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

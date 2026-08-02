# cases/3d-inflow — reference plume source-flow case

The reference workflow for the 3D inflow lineage. A rarefied plume expands from a
hemispherical inflow surface into a vacuum box; `dsmcFoam+` solves the DSMC field.

Seven other cases (`iss_solar_panels{,_2,_3}`, `solar_panel_particle_resolution_study/Fnum_*`)
were copied from this one and share its model. This is the case the regression
golden is bound to.

## Geometry

Measured from the committed mesh, not assumed:

| | |
|---|---|
| Domain | box, x ∈ [0, 5] m, y, z ∈ [−2.5, 2.5] m |
| `inflow` | hemisphere, radius **exactly 0.5 m**, centred at the origin, opening toward **+x** |
| `sym` | the x = 0 plane (5 × 5 square minus the r = 0.5 circle), symmetry |
| `vacuum` | the five outer box faces |
| Cells | 99 240 tetrahedra (204 717 faces, 192 243 internal) |
| Inflow faces | 2044 triangles |

The plume axis is **+x**. That is inferred from the one-sided hemisphere, from
`revolutionAxis "x"` in `dsmcProperties`, and from the first sampling line — it is
stated nowhere in the original code.

## Running

```bash
./Allrun                 # decomposition read from system/decomposeParDict
./Allrun --serial        # single process
./Allrun -np 4           # override the subdomain count
./Allrun --no-solve      # stop after generating 0/ — needs no OpenFOAM
./Allclean               # remove results; --dry-run to preview
```

`Allrun` performs, in this order:

```
topoSet                  # -> constant/polyMesh/sets/inflow   (skipped if unavailable)
dsmcInitialise+          # -> 0/
python runInflow.py      # overwrites 0/boundaryU, boundaryT, boundaryNumberDensity_Ar
dsmcFoam+                # or decomposePar -> mpirun -parallel -> reconstructPar
```

**The order is load-bearing.** `dsmcInitialise+` *creates* `0/`; `runInflow.py`
then overwrites three files inside it. Reversed, the inflow conditions are
silently discarded — the commented-out `touch 0/boundary*` lines in the archived
`caseFoamEx/runCases.py` were an attempt to work around exactly that.

Every command's exit status is checked, so a failed step aborts with the tail of
its log rather than letting the run continue on stale data. Output goes to
`log.<application>`, including `log.dsmcFoam+`, which is what the `monitor`
gnuplot scripts grep — no pre-refactor run script ever created it.

`Allclean` never touches `constant/polyMesh`, and matches time directories by
parsing their names as numbers rather than globbing `0*`, so a `0.org` cannot be
caught by accident. It replaces the old `clean.sh`, which ran
`rm -r 0* boundaries fieldMeasurements postProcessing processor*` with no `cd`
guard — from the repository root that deleted matching paths *from the repository
root* (finding RP-01). Recover it with
`git show pre-refactor-baseline:cases/3d-inflow/clean.sh` if you need it.

All physical inputs live in [`case.yaml`](case.yaml). No Python needs editing to
change a parameter — that is the point, and it is what broke down in the cases
copied from here.

## Meshing

Regenerating the mesh is opt-in and destructive, and `Allrun` never does it:

```bash
./Allmesh --dict-only     # write mesh dictionaries only, no OpenFOAM needed
./Allmesh --yes           # generate, mesh, checkMesh, verify -- OVERWRITES constant/polyMesh
git checkout -- constant/polyMesh    # restore the Pointwise mesh afterwards
```

`mesh.type: snappy_hex_sphere` — **`blockMesh` builds the background box only;
`snappyHexMesh` carves the inflow cavity.** The hemisphere is in no
`blockMeshDict`.

```
blockMesh                  # 20 x 20 x 20 = 8000 uniform cells of 0.25 m
snappyHexMesh -overwrite   # carve the cavity, 4 octree levels, snap to the sphere
```

Building the hemisphere in `blockMeshDict` instead needs a 5-block O-grid
(`mesh.type: block_mesh_ogrid`), and that topology has a quality floor refinement
cannot lift: five blocks meet at the projected cube corners and
`radial_grading: 10` shears every cell in the graded direction. Measured with
`checkMesh`:

| | O-grid | snappyHexMesh |
|---|---|---|
| Max / mean non-orthogonality | 66.1 (limit 70) / 31.2 | **36.9 / 11.9** |
| Max skewness | 1.79 | **0.76** |
| Cells | 48 000 | **32 272** |
| Inflow faces | 2000 quads | 5452 polygons |

Surface resolution comes from the octree, not the box:
`background_cell_size_m / 2**refinement_level = 0.25 / 16 = 0.015625 m`, 32 cells
across the radius. A coarse background with more levels is what replaces the
O-grid's `radial_grading` — fewer cells overall, on a 2.5× finer surface.
Snapping to a `searchableSphere` **primitive** rather than an STL keeps the
inflow vertices 8.3e-16 m from the true sphere.

**`./Allmesh` deletes `constant/polyMesh/sets/`** along with the other index-based
bookkeeping, because those labels address the old mesh. Re-run `./Allrun`
afterwards: it redoes `topoSet` and then `runInflow.py`, in that order.

See [docs/mesh.md](../../docs/mesh.md).

Post-processing is **not** wired up yet (finding RP-17): the shipped
`system/sampleDict` was copied from the archived 2d-wedge case and samples from
r = 1.005 m, outside this case's 0.5 m inflow surface. The `sampling:` block in
`case.yaml` is validated but not yet consumed.

## ⚠ Unresolved — needs a DSMC environment

Two settings in this case are inconsistent with what its own descendants later
did. Neither has been changed, because confirming the consequence requires
running OpenFOAM.

### HA-06 — `coordinateSystem dsmcAxisymmetric` on a 3D mesh

`constant/dsmcProperties` sets:

```
coordinateSystem   dsmcAxisymmetric;
axisymmetricProperties { maxRadialWeightingFactor 1000; revolutionAxis "x"; polarAxis "z"; }
```

This mesh is fully three-dimensional, the inflow model is φ-dependent, velocities
have three components, and there is a `sym` patch. The setting was inherited from
the archived `2d-wedge` case by directory copy (commit `51a75fe`).
`iss_solar_panels_2` — same lineage, same model — switched to `dsmcCartesian`
(commit `115d4f5`).

If dsmcFoam+ applies radial weighting on a 3D mesh, results produced here and in
`iss_solar_panels` are affected. **Verify before trusting historical output.**

### HA-07 — `type Unspecified;` on the inflow patch

`constant/polyMesh/boundary` declares:

```
inflow { type Unspecified; nFaces 2044; startFace 192243; }
```

`Unspecified` is a Pointwise V18.5R2 export artifact, not a registered OpenFOAM
patch type. `iss_solar_panels_2` and all three `Fnum_*` cases carry `type patch;`
on the *same* mesh.

Either v1706 tolerated it, or these cases were run from a locally corrected
`boundary` that was never committed — in which case **this case cannot be run
exactly as committed**. This is the most likely blocker to a first successful
reproduction. A generated mesh emits `type patch;` and avoids the question
entirely.

## Known frozen defects

Every value marked `LEGACY` in `case.yaml` reproduces pre-refactor behaviour
bit-for-bit. Several are known to be wrong and are frozen deliberately so the
extraction could be proven behaviour-preserving. The most consequential:

| ID | Frozen behaviour |
|---|---|
| **SM-01** | `p0_pa` is 475 psi (3.275 MPa). The original set 34 500 Pa (5 psi) in `physical_props` and then never read it — the hard-coded value won. |
| **SM-02** | Density uses `legacy.rhoN_T0_K` = 300 K while velocity and temperature use `stagnation.T0_K`. Identical here; a factor of 10 apart in the resolution-study cases. |
| **SM-03** | The model runs on N₂ (γ = 1.4, M_w = 28.0134) while the solver simulates argon. v_ℓ = 788.16 m/s instead of 558.89 m/s. |
| **SM-04** | θ is measured from **+z**, not from the +x plume axis, and the angular dependence is a separable product of two angles in different planes. |
| **SM-05** | The normalisation constant normalises a *different* angular function from the one the density uses. |
| **AD-03** | `cylinder` and `plate` boundary entries are written into every field file; this mesh has neither. Loading `case.yaml` warns about it. |

Full reconstruction, confidence classification and open questions:
[docs/source-flow-model.md](../../docs/source-flow-model.md).

**The source-flow model has no recoverable citation.** The only literature
reference in the repository covers the archived 1d/2d-wedge validation case.

## Files

| Path | |
|---|---|
| `case.yaml` | all physical and mesh parameters — the only file to edit |
| `Allrun` / `Allclean` | run the case / reset it |
| `Allmesh` | regenerate `constant/polyMesh` from `case.yaml`; destructive, opt-in |
| `runInflow.py` | generates `0/` boundary fields; identical in every case |
| `constant/polyMesh/` | the Pointwise mesh, deliberately tracked (see docs/mesh.md) |
| `mesh/*.pw` | binary Pointwise project files; not reproducible |
| `processInflowData.py.legacy` | the pre-refactor script, kept for provenance. **Do not run.** |
| `getMeshStats.py`, `printInflow.py` | superseded by `plumetools`; left untouched |
| `system/blockMeshDict.2d-wedge-INVALID` | wrong-geometry dict that would destroy the mesh (FD-08) |

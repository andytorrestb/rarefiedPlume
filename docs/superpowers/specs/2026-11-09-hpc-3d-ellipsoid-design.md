# HPC 3D CFD–DSMC plume with a conservative fitted ellipsoid

Status: proposed for owner approval, 2026-11-09.

## Intent and agreed choices

Extend the working fitted-interface workflow to a fully 3D circular-source plume.
The owner selected native snappyHexMesh meshing around a searchableSphere
ellipsoid, exit Kn = 0.001, a domain extending 8 D downstream and ±5 D in both
transverse directions, and an HPC node with 192 cores and 1.5 TB RAM. The normal
submitted job must execute CFD, fitting, mesh preparation, DSMC initialization,
startup, sampling and post-processing. It must not end at prepare-dsmc.

Keep development checks minimal because the owner has 14% weekly usage remaining.
Do not run either solver locally, a production mesh, a convergence study, a
benchmark, a parameter sweep or the full repository test suite. Prepare the job
for the owner to submit through the existing PBS environment; remote submission
is outside this change. Preserve the working planar case and other 3D profiles.

## Approaches

1. **Native snappyHexMesh with a fitted searchableSphere ellipsoid (selected).**
   Reuse the repository's native cavity-meshing pattern and the existing physics,
   breakdown and DSMC boundary model. This keeps the extension focused.
2. **Projected multiblock 3D mesh.** More control of structured cells, but adds
   substantial geometry construction and validation. Rejected in favor of the
   smaller native-mesher extension.

The fitted surface has equal transverse radii for the circular source, but both
solvers and diagnostics operate on full 3D fields. Do not rotate the planar CFD
solution, average away azimuthal violations, or use a wedge or radial weights.

## Source, domain and initial profile

Retain D = 0.2 m, T = 300 K, speed ratio 2, the existing argon VHS parameters and
consistent CFD transport; set exit Kn = 0.001. The box is x = 0…1.6 m and
y,z = −1…1 m. Source, exitPlane and farField keep their existing physical roles.

Proposed initial profile name: hpc-kn1em3-ellipsoid. Builder proposals for its
Cartesian CFD mesh: D/40 core cells, core extending to x = 2 D and |y|,|z| = D,
and 80 graded cells in each outer segment. The current planner gives
160 × 240 × 240 = 9,216,000 CFD cells. This is a starting HPC case, not a
mesh-converged production result. Preserve the stable Minmod reconstruction,
background-floor policy and viscous timestep factor 0.25; do not inherit the
aggressive factor 10 from the old verification profile.

Retain the configured startup/sampling-transit and parcel-occupancy policies.
DSMC timestep and particle weight must be recomputed from the actual DSMC mesh,
mapped CFD states, collision time and tracking widths. Refinement controls must
support finer cells around the interface and uncovered source; report local
cell-size/mean-free-path and occupancy metrics without claiming convergence.

## Conservative 3D fit

Use the unchanged Kn_GLL and abs(B_exp) tests, configured thresholds and
source-connected 3D continuum component. Invalid cells and disconnected pockets
are excluded, including interior holes. The lip diagnostic exclusion is not a
fitting exemption.

Fit ((x−xc)/a)^2 + (y^2+z^2)/b^2 ≤ 1, clipped to x ≥ 0, with 0 ≤ xc < a.
Maximize its downstream retained volume, subject to exclusion of every rejected
CFD voxel and source intersection within the physical disk. Use the same
quarter-core-cell default inward clearance and a positive outer-domain gap.
Report the achieved local optimum and optimizer termination; do not claim a
global optimum. If no feasible source-connected ellipsoid exists, report it and
stop, without substituting an unsafe surface or an unrequested full-DSMC case.

Use whole voxel footprints, not just cell centres. For equal transverse axes,
constraints can be reduced to the nearest excluded radial footprint in each
axial slab using all azimuths; independently verify the final result against the
3D excluded voxels in bounded-memory chunks. This avoids millions of redundant
optimizer constraints without relaxing conservatism.

Bind breakdown, fitting and transfer to field/configuration/mesh hashes, including
checks for data changed during calculation. Reuse the existing provenance guard
rather than accepting cached classifications from a different solution.

## Curved mesh and field transfer

Create a separate 3D DSMC mesh outside the ellipsoid, using searchableSphere with
radius (a b b), native snapping and a background box covering the enlarged
domain. Reuse existing snappy dictionary/quality conventions where applicable.
Use configured local refinement and bounded global cell counts; do not run an
automatic mesh-search loop or introduce a meshing dependency.

Parallel meshing must not leave processor directories that are mistaken for
initialized DSMC checkpoints. Use an isolated meshing workspace, or a verified
mesh-only lifecycle, reconstruct before actual-face mapping, then initialize and
decompose the DSMC case with its final fields. Record mesh identity so retries
cannot silently combine a new fit with an old particle cloud.

Generalize the containing-cell sampler to the rectilinear 3D CFD mesh, including
coordinate-roundoff tolerance. Supply bounded cell states at actual interface
faces, preserve their CFD-cell provenance, derive n = rho/m consistently, and
use actual areas and correctly oriented normals. Do not weaken the existing
coincident-face mapping used by the original subsetMesh path.

Check actual snapped facets against the certified continuum region and report
geometry deviations. If a facet enters an excluded CFD voxel or exceeds the
allowed geometry tolerance, stop before DSMC initialization. Certifying an
analytic enclosure of all facet vertices against the excluded voxels is an
acceptable conservative implementation; checking face centres alone is not.
Preserve direct injection from uncovered source faces.
Measure both meshes' source areas and covered/direct-source coverage against
the intended disk, with an explicit geometry tolerance; do not renormalize
density or omit source flux to disguise area error. Use actual 3D geometry for
closure, fluxes, normal Mach, return fraction and outward overlap layers.

## Suitability and complete execution

The larger domain and lower exit Kn are intended to improve the interface; they
do not predetermine its suitability. Preserve the existing measured criteria.
The default job continues through dsmcInitialise, parallel dsmcFoam startup and
sampling, and post-dsmc when the surface passes. A failed criterion remains an
explicit stop with diagnostics. The existing explicitly requested smoke-only
override remains available and labelled; the new profile does not silently
enable it. Infeasible geometry cannot be overridden.

Post-processing must use the actual remeshed DSMC coordinates, not assume CFD
cell/face identity. Reuse existing moment readers and flux accounting. Provide
an ellipsoid/breakdown overlay or section and .foam entry points. Comparisons
with the circular-source reference can be labelled appropriately for a
collisional 3D case; planar-only assumptions must not leak into this path.

## HPC resources and delivery

Add a dedicated PBS run3DCase.sh and the new profile; leave the current local 2D
launcher and existing PBS profile choices intact. Reuse the known longq,
120-hour ceiling, cluster environment and library build steps. Locate paths from
the submitted checkout; do not embed Windows or temporary worktree paths.

Use 192 MPI ranks. The physical node limit is 1.5 TB; reserve headroom rather
than consuming every reported byte. Proposed initial explicit PBS memory request
is 1200gb, below the stated node capacity, with a configurable estimate ceiling
and overhead allowance. State units explicitly because the current resource
estimator labels GiB as GB. Account for both meshes, analysis temporaries, MPI
and parcels; report estimates before launch and stop if the configured budget
is exceeded. Do not claim measured cluster RAM or runtime.

The default launcher must call the complete workflow and resume complete
checkpoints after interruption. Document the exact qsub command, profile,
geometry/refinement settings, memory controls, outputs and explicit smoke flag.

## Essential verification only

- Focused numerical fit/3D sampling checks: excluded voxels, correct containing
  cells, valid face states, and provenance refusal. Reuse existing small helpers.
- One small synthetic 3D native mesh-and-transfer check, including checkMesh,
  actual facet location, normals, source coverage and geometry closure. No full
  HPC mesh or solver run.
- Shell syntax and a stubbed PBS/runner check proving that the default path
  reaches DSMC initialization/solve and post-processing, with 192-rank settings.
- Run only directly affected existing tests where shared code changed. One
  focused final review; no repeated passing checks without a new concern.

Remaining refinement levels, snapping tolerances, footprint implementation and
memory-accounting details are builder choices within the stated constraints.

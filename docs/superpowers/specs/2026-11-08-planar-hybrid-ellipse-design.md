# Planar CFD–DSMC prototype with a conservative elliptical interface

Status: proposed for owner approval, 2026-11-08.

## Intent and agreed scope

Replace the jagged continuum-breakdown interface with a conservative elliptical
approximation and a mesh that follows it. Supply a complete, inexpensive planar
CFD → breakdown → fitted interface → DSMC workflow for manual assessment on a
laptop. The owner accepts planar plume physics and selected the constrained fit.
Reuse the current case's physical source state, gas, breakdown definitions,
thresholds, source connectivity and one-way inflow model.

This is a prototype, not a validation study. Existing 3D profiles remain available.
Do not run simulations, convergence studies or parameter sweeps during development.
Do not change unrelated local work or existing generated runs.

## Approaches

1. **Constrained ellipse fit (selected).** Maximize retained continuum area subject
   to the original continuum-valid region and a small inward clearance. This
   makes the conservative requirement an explicit constraint.
2. **Least-squares fit followed by shrinking.** Less fitting machinery, but
   potentially sacrifices considerably more useful continuum to satisfy the
   same constraints. The owner chose option 1.

Merely selecting cells inside an ellipse would preserve the Cartesian staircase.
The selected approach therefore includes a separate mesh following the curve and
sampling CFD states onto its new interface faces.

## Planar physics and proposed laptop defaults

Use x along the plume, y across it, and exactly one cell through a uniform depth
in z. Front and back are `empty` in both solvers and all fields. The source is a
slit of height D, rather than a circular disk. Preserve the existing exit n, T,
U and Kn based on D; preserve VHS argon and the existing CFD transport model.
DSMC retains three molecular velocity components and translational degrees of
freedom, while tracking only in x and y. There is no radial weighting.

Implementer proposals for the initial opt-in profile: Kn = 0.01, x extent 2 D,
y extent ±1 D, core spacing D/20 and depth D/20. This gives about 1,600 CFD
cells before any mesh-specific adjustments. Keep the established stable CFD
limiter, pressure floor and timestep cap; use the existing DSMC occupancy and
collision/transport timestep policy. Start with serial execution, with up to
four ranks available. Runtime estimates are estimates, not promises.

Use actual extruded volumes to set parcel weight and occupancy. Report fluxes
for the chosen depth, optionally also per unit depth. Do not label them circular
nozzle totals or compare with the current 3D disk collisionless solution.

## Conservative fit

Compute the same `Kn_GLL` and `abs(B_exp)` tests on unaveraged CFD values. Use the
same mean-free-path convention, configured thresholds and source-connected
component. Treat invalid cells and disconnected pockets as outside the allowable
continuum region. Do not average the two sides of the plume to hide a violation,
loosen thresholds or use the lip diagnostic exclusion as a fitting exemption.

Fit an axis-aligned ellipse with parameters `(xc, a, b)`, centred at `(xc, 0)`:

`((x - xc)/a)^2 + (y/b)^2 <= 1`, restricted to x >= 0.

Maximize the area of that downstream, clipped ellipse. Require positive axes,
intersection with the source at x = 0, no outer-boundary contact, and its
intersection with x = 0 to lie within the physical source slit. A partially
covered slit is allowed: its uncovered portions inject the original source
state directly into DSMC, as in the current workflow.

Conservatism means the retained ellipse lies inside the union of accepted CFD
cell footprints, not merely that selected cell centres look valid. Check
ellipse intersection with excluded cell footprints, including interior holes,
and with domain limits. Apply a configurable inward clearance; the proposed
default is 0.25 core-cell widths from excluded cells. Source contact itself is
exempt from domain clearance. This is conservative relative to the resolved
CFD-cell classification, not a claim about unresolved physical gradients.

Use existing SciPy for a small deterministic constrained optimization, followed
by an independent geometric feasibility check. Report its termination and
achieved area; do not claim a proven global optimum. If no feasible,
source-connected ellipse exists, report that result and stop interface setup.
Do not silently substitute an unsafe ellipse or a full-DSMC run.

## Geometry and meshing

Represent the fitted ellipse using OpenFOAM's native `searchableSphere`
ellipsoid support (three radii). Use its planar section as the geometric
reference for the curved interface. The final computational boundary is an
extrusion of that section, with identical x/y points at both depth planes.
It must not retain the depth-dependent curvature of a 3D ellipsoid.

Generate a small, deterministic mesh outside the clipped ellipse, using curved
block edges and uniform extrusion. Prefer native `blockMesh` projection with
plane constraints or a midplane construction followed by extrusion. The builder
may select the least complex supported construction after checking local
OpenFOAM syntax; no 3D snappyHexMesh pipeline or new meshing dependency is needed.
The ellipse is faceted only at the intended interface resolution, with no
cell-selection staircase. Verify the actual interface facets against the same
conservative region, not just the ideal curve.

The CFD precursor can remain a simple planar Cartesian mesh. Keep patch roles
`source`, `exitPlane`, `farField` and `interface`; add explicit front/back empty
patches. Use the real mesh geometry for normals, face areas, volumes and closure.

## Field transfer and suitability

The new mesh does not share CFD face or cell identities. Give planar setup and
mapping an explicit path instead of weakening the existing 3D coincident-face
checks. Locate interface sample points in accepted CFD cells and use a bounded
local reconstruction that remains finite and positive. No extrapolation from
outside the accepted region and no replacement with invented exit states.
Store the contributing CFD cells and interpolation information for inspection.

Derive `n = rho/m` and `p = n k T` consistently. Feed the actual new face areas
and inward DSMC normals to the existing half-range inflow model. Compute and
report normal Mach, returning mass fraction, overlap availability, invalid
states, geometric closure and flux diagnostics on the new surface. Retain the
current suitability policy and explicit smoke-only override. A conservative
ellipse alone does not imply suitable one-way coupling or exact NSF flux transfer.

Initialize the DSMC domain with the existing near-vacuum procedure. Resource
estimates use its own cell geometry and sampled finite CFD states where needed;
they must not require a one-to-one cell map. Planar post-processing must use its
own mesh and coordinates rather than assuming a subset of the CFD cells.

## Integration and deliverables

Add an opt-in planar configuration and fitted-ellipse interface mode through the
existing hybrid CLI, with focused planar fitting, mesh and mapping helpers.
Reuse the existing gas, breakdown, half-range, inflow and DSMC field writers.
Update the staged runner to build the separate planar DSMC mesh instead of
calling `subsetMesh` on this path. Leave default 3D behavior intact.

Deliver one laptop profile and exact commands for setup, CFD, interface/mesh
preparation, DSMC and basic planar post-processing. Save fit parameters, fit
clearance, retained-area statistics, source coverage, geometry and suitability
diagnostics. Provide a basic overlay of the raw breakdown boundary and fitted
ellipse, plus a `.foam` entry point for manual inspection. Existing 3D-only
analytical comparisons must be marked not applicable to planar runs.

## Minimal verification

- One focused numerical check of a known ellipse and an excluded interior
  pocket: a fitted ellipse must stay conservative or return infeasible.
- One small planar setup/mapping check covering empty patches, actual face
  normals, positive inflow and source coverage, plus relevant existing tests
  only if shared behavior changes.
- If the local OpenFOAM environment is accessible, run mesh generation and
  `checkMesh` once on a tiny synthetic feasible ellipse. Check matching depth
  planes, positive volumes, interface orientation and no unintended openings.

No CFD or DSMC solver runs, benchmarks, sensitivity sweeps or full test suite.
The owner will run the supplied case and assess the result. Document any native
check that could not be run; do not claim the simulation itself was validated.

## References checked

- Existing implementation: `plumetools/hybridplume/{breakdown,interface,mesh,
  mapping,workflow,cfd,dsmc}.py` and `cases/cai2012-hybrid/Allrun`.
- [OpenFOAM ellipsoidal searchableSphere support](https://www.openfoam.com/news/main-news/openfoam-v20-12/post-processing).
- [OpenFOAM DSMC reduced-dimensional tracking with retained 3D molecular velocity](https://api.openfoam.com/2506/DSMCParcel_8C_source.html).

Remaining meshing syntax and bounded-reconstruction details are builder choices
within these constraints. Axisymmetric DSMC, production-resolution accuracy,
two-way coupling and running the simulations are outside this prototype.

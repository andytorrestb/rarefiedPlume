# Fitted capsule interface on a larger, coarser, stabilised 3D CFD

Status: design approved by the owner in conversation, 2026-11-11; written spec
awaiting review.

## Intent and agreed choices

Job 79505 (profile `hpc-kn1em3-ellipsoid`) stopped at the DSMC stage because the
fitted interface was unsuitable: `return_mass_fraction = 0.9873` (limit 0.05).
Diagnosis from that run:

- The fitted spheroid has `xc ≈ a`, so it meets the exit plane at a point and
  covers 0.13 % of the source. DSMC injects 99.9 % of the jet directly; that gas
  re-enters the spheroid through its backward-facing flanks (98.4 % of the
  reverse flux comes from the 10.8 % of interface area with flow into the CFD
  region, near x ≈ 0.15 D).
- The retained continuum region next to the exit is a column of radius ≈ 0.5 D
  (the full source radius) out to x ≈ 3 D. A spheroid anchored on the source
  cannot fill it without entering lip voxels; a cylinder can.
- A throwaway probe (analytic surfaces on the job's CFD state, same half-range
  flux code) gave return fractions of 0.049–0.055 for R = 0.49 D cylinder-plus-
  tip shapes, 0.22 for R = 0.40 D. About two thirds of the remainder comes from
  the cylinder wall within 0.25 D of the exit, where U_r ≈ 0. The tip contributes
  almost nothing. Passing the 5 % limit is therefore expected to be borderline.

The owner also judged the CFD domain much too small and the cell count far too
high, and asked for the coarsest mesh that gives stable results. A mesh
sensitivity study is explicitly out of scope and recorded as follow-up work.

Agreed choices:

1. New interface mode `fitted-capsule`: a cylinder from the exit plane capped by
   a spheroidal tip **tangent** to it (tip transverse semi-axis = cylinder
   radius). Fitted parameters (R, L, a).
2. Selection: among all certified-conservative capsules, the **largest volume
   whose pre-mesh flux screen passes** the configured return-fraction and
   subsonic-area limits.
3. Architecture: a shared shape interface; the existing spheroid becomes one
   implementation (bit-identical behaviour), the capsule another. All downstream
   utilities are reused unchanged.
4. CFD domain 30 D × 15 D (x_max × half-width). The DSMC background mesh keeps
   sharing this box.
5. Outer grading cap raised as needed for a continuous core/outer transition.
6. Stability via a kinematic-viscosity cap in the project's own transport
   library, only allowed to act outside the retained continuum region.
7. Existing profiles and the `fitted-ellipsoid` mode are unchanged.

Success: a new profile `hpc-kn1em3-capsule` runs CFD to completion on the coarse
mesh with the stability criteria below, fits and meshes a capsule, and reaches
the DSMC stage. Passing the suitability check without the smoke-only override is
the goal, but the design does not guarantee it (see probe numbers above).

## Part 1 — CFD domain, mesh and stability

### 1a. Domain and mesh (configuration only)

Profile `hpc-kn1em3-capsule` inherits `hpc-kn1em3-ellipsoid` except:

```yaml
domain: {x_max_over_D: 30.0, half_width_over_D: 15.0}
mesh: {geometry: cartesian, core_cell_over_D: 0.05, core_x_over_D: 3.0,
       core_half_over_D: 1.5, outer_cells: 40, max_outer_expansion: 60.0}
cfd: {end_time_transits: 1.0, viscosity_cap: {enabled: true, factor: 1.0}}
interface: {mode: fitted-capsule, ...}   # Part 2
```

Planner output for these values (computed with `mesh.plan` on 2026-11-11):

| quantity | value |
|---|---|
| cells | 1.96 M (job 79505: 9.22 M) |
| max neighbour-cell ratio | 1.107 (first outer cell = core cell) |
| outermost cell | ≈ 2.7 D axial, ≈ 1.1 D lateral |
| source staircase area error | 0.59 % (limit 1 %) |

0.05 D is the coarsest core cell that satisfies the 1 % source-area check
(0.04 D: 1.4 %, 0.0625 D: 3.5 %, 0.08 D: 2.2 %, 0.1 D: 1.9 %). With the default
cap of 8 the outer blocks cannot start at the core cell size and the first outer
cell jumps 4× at the core boundary; a cap of 60 removes the jump (total
last/first ratio needed ≈ 53 axially).

One transit of the 30 D domain (≈ 8.5 ms) equals ≈ 3.7 transits of the old 8 D
domain, longer than the 2.75 used for job 79505.

### 1b. Viscosity cap

Root cause of the job-79485 crash, confirmed on job 79505's final write: the
region beside and behind the plume (x 0–2.3 D, r 1.7–5.4 D) drains to
≈ ρ_bg/100. The tensor-product mesh carries the core's x spacing to large r, so
those cells have both the smallest dx and the lowest density; the explicit
viscous number `μ dt/(ρ dx²)` reaches ≈ 30. Coarsening does not help (the dt cap
scales with dx² too) and a larger domain enlarges the drained region.
`rhoCentralFoam` in v2512 has no fvOptions/fvConstraints hook.

Design:

- C++: in `applications/thermophysicalModels/hybridPlumeThermo`, add a transport
  wrapper `cappedTabulatedTransport` around `tabulatedTransport`:
  `mu_eff(p,T) = min(mu_tab(T), nu_max · p/(R T))`;
  `kappa_eff = kappa_tab(T) · mu_eff/mu_tab` (Prandtl number preserved).
  `nu_max` is read from the transport dictionary. Register one more `psiThermo`
  combination (`transport cappedTabulated`) exactly like the existing file.
- Python (`cfd.py`): when `cfd.viscosity_cap.enabled`, write
  `nu_max = factor · viscous_dt_factor · h² / dt_visc` (= factor · ν_bg, i.e. no
  cell can exceed the viscous number the dt cap was designed for) and select the
  capped transport. Report `nu_max` in `setup.json`. Default disabled.
- Accountability (CFD report / breakdown step): count and write a
  `viscosityCapped` cell field (`mu_tab > nu_max ρ`) at the breakdown time.
  **Error** if any capped cell is inside the retained continuum region: the cap
  may only act where CFD is already invalid and DSMC takes over.

### 1c. "Stable" acceptance for the coarsest run

1. Completes `end_time_transits` with no crash and no T < 0.
2. CFD report mass `residual_over_reference` ≤ 1e-3.
3. Zero capped cells in the retained region (1b).
4. Retained-region volume fraction changes ≤ 2 % between the last two writes
   (new small check in the breakdown step; reported, and a failure is an error).

## Part 2 — Fitted capsule interface

### 2a. Shape interface (`shapes.py`)

Two frozen dataclasses with a common method set:

| method | purpose |
|---|---|
| `q(points)` | scaled inside measure; 1 on the surface, < 1 inside |
| `deviation(vertices)` | max \|√q − 1\| for the snapped-facet tolerance |
| `outward_gradient(points)` | reversed-normal check in `_geometry` |
| `verify(lower, upper, retained, clearance)` | certification against excluded voxels |
| `grown(factor)` | enclosure for snapped-facet certification |
| `section(n)` | z = 0 section points for plots |
| `volume()`, `base_radius()` | reporting, source coverage |
| `to_json()`, `from_json()`, `params()` | persistence and fit identity |

- `Spheroid(xc, a, b)` wraps the current code in `ellipsoid.py`
  (`verify` = `verify_ellipsoid`). Fitting (`fit_ellipsoid`) is unchanged.
- `Capsule(R, L, a)` = cylinder {0 ≤ x ≤ L, r ≤ R} ∪ spheroid centred (L,0,0)
  with semi-axes (a, R, R), clipped at x ≥ 0.
  `q = (r/R)²` for x ≤ L, `((x−L)/a)² + (r/R)²` for x > L.
  `verify` = box-to-cylinder test for the base plus `verify_ellipsoid` on the
  tip spheroid. `grown(g)` scales R and a by g and keeps L.

`spatial.py` (`certify_facets`, `_geometry`, `plots`, mesh seed check) and
`fitted._identity` consume a shape instead of `fit["xc"], fit["a"], fit["b"]`.
The fit identity hashes `shape.params()`; for the spheroid this is the same
`[xc, a, b]` list as today, so existing identities do not change.

### 2b. Fit and flux screen (`capsule.py`)

Exact sweep, no optimizer:

- Feasible R: `0 < R ≤ min(source radius, domain radial limit) − clearance`,
  sampled at `interface.fit.screen.radius_samples` values from R_max downward.
- For each R, `L_max(R)` = the axial start of the first excluded slab whose
  radial footprint distance is < R (with clearance), capped by the domain.
- For each (R, L) on `length_samples` values in (0, L_max], the tip length
  `a_max(R, L)` is the minimum over excluded slabs beyond L of
  `dx / √(1 − ρ²/R²)` (ρ < R), capped by `x_max − L` minus domain clearance.
- Every candidate (R, L, a_max) is certified with `Capsule.verify` before it
  counts as feasible.

Flux screen per feasible candidate: sample the analytic surface (cylinder wall +
tip; the base disc is the CFD source, not interface), locate states with
`fitted.sample_containing_cells` (rectilinear grid built once and cached), and
integrate `halfrange.incoming_flux` forward and reverse to get the return
fraction and subsonic area fraction.

Selection: largest volume with screened return fraction ≤
`max_return_mass_fraction` and subsonic area ≤
`max_subsonic_normal_area_fraction`. If none passes, take the lowest-return
feasible candidate and set `screen.status = "no_passing_candidate"`; the actual-
mesh suitability check stays authoritative and refuses without the smoke-only
override. All candidates and the selection are written to
`results/capsule.json` and summarised in `interface.json`.

New config keys (defaults shown):

```yaml
interface:
  fit:
    screen: {radius_samples: 12, length_samples: 24,
             axial_points_per_D: 40, azimuthal_points: 90}
```

Existing `clearance_cells`, `geometry_relative_tolerance` and
`max_source_area_error` are reused unchanged.

### 2c. Meshing (`capsule_mesh.py`)

- Write a watertight, outward-oriented STL `constant/triSurface/interface.stl`:
  the cylinder extends to x = −δ (δ = one background cell) and is capped there,
  so the domain's exit plane cuts it cleanly. Tessellate so the chord error is
  ≤ 0.1 × `geometry_relative_tolerance`.
- snappyHexMeshDict references it as `triSurfaceMesh` named `interface` at
  `dsmc.mesh.surface_level`; source refinement annuli, background block, seed
  logic (seed check via `shape.q`) and the `maxLocalCells` handling are reused
  from `ellipsoid_mesh.py` (factored into shared helpers, not copied).
- Feature edges: a `surfaceFeatureExtract` step for the rim circle where the
  cylinder meets the exit plane and the cylinder–tip junction, referenced in
  snappy's `features`. Allrun runs it before snappyHexMesh for spatial fitted
  modes that provide an STL.
- If snapped facets exceed the tolerance, the existing certification refuses
  loudly. No silent degradation.

### 2d. Routing

- `config.INTERFACE_MODES` gains `fitted-capsule`; cartesian geometry required.
- A single helper `config.is_spatial_fit(cfg)` replaces every
  `== "fitted-ellipsoid"` test in `workflow.py`, and the Allrun condition uses
  the same predicate (via `cfg_value`), so the modes cannot drift apart.
- `spatial.interface` dispatches on the mode to the spheroid or capsule fitter;
  the fit file is `ellipsoid.json` or `capsule.json` respectively.

### 2e. Reused unchanged

`fitted._assess` (all suitability criteria), source-coverage and closure checks,
`mapping`, `interface_fields.npz`, `interface.vtk`, DSMC setup and fields,
resource estimate and memory guard, `dsmc_post`, crossings, report, and the
192-rank PBS path (`PROFILE=hpc-kn1em3-capsule ./run3DCase.sh`).

## Testing

Unit tests (no solver runs):

- Capsule geometry: `q`, gradient direction, volume, `section`.
- Fit on synthetic voxel sets with known answers: unbounded continuum column
  (R → source limit, L → domain limit); a step at x0 (`L_max = x0`, expected
  `a_max`); exclusions touching only the tip.
- Certification: the fitted capsule passes; the same capsule grown by 1 % fails.
- Screen selection on a synthetic flow (axial core plus radial outflow): picks
  the largest passing candidate; reports `no_passing_candidate` correctly.
- STL: closed, consistently outward-oriented, chord error within bound.
- Spheroid regression: `ellipsoid.json`, certification output and fit identity
  identical before and after the refactor.
- Config/workflow routing for `fitted-capsule`; `is_spatial_fit`.
- Viscosity cap: `nu_max` value written to the dictionary; capped-cell
  accountability raises when a capped cell is retained.
- C++: the new thermo library compiles (Allwmake).
- Existing hybridplume ellipsoid/interface/workflow unit tests keep passing.

## Order of work

1. Part 1 (profile, viscosity cap, stability checks) → owner runs CFD and
   breakdown on HPC.
2. Part 2 (shapes, capsule fit/screen, STL meshing, routing) → fit on the new
   breakdown.
3. DSMC on the owner's submission.

## Risks

- **DSMC in the 30 D × 15 D box**: ≈ 3.4 M background cells at 0.2 D before
  refinement. The existing resource estimator and memory guard report or refuse.
  If refused, a coarser or graded DSMC background is the next step; not built now.
- **Coarser core can move the breakdown surface** (gradient-based indicators are
  resolution-dependent). Covered by the follow-up sensitivity study.
- **Borderline return fraction** (≈ 0.05 expected): intrinsic to the near-exit
  wall where U_r ≈ 0; the screen picks the best certified shape but cannot
  guarantee a pass.
- **Snappy at the exit-plane rim** (90° edge at r = R): mitigated by feature
  edges; refused by certification if it fails.

## Follow-up (out of scope)

- Mesh sensitivity study: core cell (0.05 D vs finer), outer cell count/cap,
  domain size, and the viscosity cap's effect on the breakdown surface.
- Graded DSMC background mesh if the 30 D box is too expensive.

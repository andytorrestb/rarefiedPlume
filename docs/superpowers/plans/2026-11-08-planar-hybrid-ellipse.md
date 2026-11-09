# Planar Hybrid Ellipse Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Deliver a laptop-scale planar CFD–DSMC case with a conservative fitted ellipse and a mesh following its interface.

**Architecture:** Keep the existing hybrid CLI and 3D defaults. Add opt-in planar source/mesh handling, an isolated constrained fit, and a separate planar DSMC meshing/mapping path. Reuse the existing gas, breakdown, half-range fluxes, boundary model and field writers.

**Tech Stack:** Existing Python, NumPy, SciPy, PyYAML and optional Matplotlib; installed OpenFOAM v2512 (`blockMesh`, `checkMesh`, `rhoCentralFoam`, `dsmcFoam`). No dependencies to install.

**Spec:** `docs/superpowers/specs/2026-11-08-planar-hybrid-ellipse-design.md` (owner approved).

## Global Constraints

- Exactly one cell through a uniform depth in z; front and back are `empty`.
- Initial profile: Kn = 0.01, x extent 2 D, y extent ±1 D, core spacing D/20 and depth D/20.
- Preserve the original breakdown thresholds and source connectivity; invalid cells and disconnected pockets are excluded.
- Proposed default inward clearance: 0.25 core-cell widths from excluded cells.
- A partially covered slit injects the original source state over uncovered portions.
- No CFD or DSMC solver runs, benchmarks, sensitivity sweeps or full test suite.
- Do not change unrelated local work or existing generated runs.

## Review Focus

- An excluded pocket inside an otherwise valid region must constrain the entire ellipse (Task 1).
- A source with no feasible symmetric source-connected ellipse must stop with a useful report (Task 1).
- Geometry projection must preserve identical front/back x/y coordinates and positive volumes (Task 2).
- Remeshed faces must sample their containing valid CFD cells; unmatched or invalid samples must fail (Task 3).
- Partial source coverage must preserve the full prescribed slit, and planar reports must exclude 3D disk references (Tasks 2–4).

## File Structure

- New `plumetools/hybridplume/ellipse.py`: ellipse geometry, constrained fit and independent footprint feasibility.
- New `plumetools/hybridplume/planar_mesh.py`: curved planar DSMC block mesh.
- New `plumetools/hybridplume/planar.py`: planar interface sampling/assessment, DSMC preparation and basic post-processing.
- Modify `config.py`, `physics.py`, `mesh.py`, `cfd.py`, `workflow.py`, `report.py` in the same package only where an explicit planar branch is needed.
- Modify `cases/cai2012-hybrid/{hybrid.yaml,Allrun,README.md}`: opt-in laptop profile and staged commands.
- New `tests/unit/test_hybridplume_planar.py`; adapt the existing 3D-profile test to distinguish the new explicitly planar profile.

### Task 1: Conservative ellipse fitting

**Files:** create `plumetools/hybridplume/ellipse.py` and the numerical fit cases in `tests/unit/test_hybridplume_planar.py`.

**Interfaces:**
- `fit_ellipse(lower: ndarray, upper: ndarray, retained: ndarray, *, source_radius: float, clearance: float) -> dict`: 2D cell footprint bounds in metres; return status, `(xc, a, b)`, clipped area and feasibility/optimizer diagnostics.
- `footprint_clearance(parameters: ndarray, lower: ndarray, upper: ndarray) -> ndarray`: minimum squared ellipse-coordinate distance to each rectangular footprint.
- `ellipse_points(fit: dict, count: int) -> ndarray`: downstream arc from upper to lower x=0 intersection.

- [x] Write and run failing tests for a resolved synthetic ellipse, an excluded interior pocket, and missing central source connectivity. Assert accepted fits never intersect an excluded footprint; infeasible inputs return no ellipse.
- [x] Implement dimensionless SciPy constrained optimization over `(xc, a, b)`, maximizing the analytically clipped area. Require `0 <= xc < a`, source intersection within the slit, and separation from excluded footprints inflated by the clearance. Use a small feasible source seed, bounded deterministic optimization and independent final feasibility; record local termination without claiming global optimality.
- [x] Run only these fit tests and verify they pass. Commit the fitting deliverable and tests.

### Task 2: Opt-in planar setup and curved mesh

**Files:** modify `config.py`, `physics.py`, `mesh.py`, `cfd.py`, `hybrid.yaml`; create `planar_mesh.py`; extend the same focused test file and adjust `test_hybridplume_workflow.py::test_every_profile_validates_and_is_full_3d`.

**Interfaces:**
- Configuration: `mesh.geometry: planar`, `mesh.depth_over_D: 0.05`, `interface.mode: fitted-ellipse`, `interface.fit.clearance_cells: 0.25`; profile `planar-laptop-kn0p01`, serial by default and at most four ranks.
- `render_dsmc_mesh(cfg: CaseConfig, source: SourceState, fit: dict) -> str`: native blockMeshDict for the exterior of the clipped ellipse.
- Existing `mesh.plan`, `mesh.predicted_source`, `mesh.render_block_mesh_dict`, `mesh.render_topo_set_dict` and `cfd.render_fields` retain their signatures.

- [x] Write and run failing planar setup assertions: 40×40×1 cells for the initial profile, source area `D * depth`, explicit empty patches in mesh and CFD fields, and rejected incompatible geometry/mode settings. Keep assertions that existing profiles remain 3D.
- [x] Add optional planar depth to the hybrid source representation so all existing ideal source fluxes use the slit area. Add one-cell z planning, slit selection and empty boundaries while preserving default 3D behavior.
- [x] Implement a small multi-block exterior mesh from the upper x=0 ellipse intersection around its downstream tip to the lower intersection. Split the source-side radial edges at y=±D/2 so remaining slit faces and exit-plane faces have exact separate patches. Use actual boundary lengths to choose subdivisions.
- [x] Use `searchableSphere` with a vector radius and projected curved edges constrained to front/back planes. Compensate the ellipsoid radii so both depth-plane intersections have the fitted in-plane semi-axes; reject any depth-plane mismatch. Keep mesh construction separate from the Cartesian CFD mesh.
- [x] Run the focused setup assertions. Generate one tiny feasible CFD/curved-DSMC mesh fixture and run native `blockMesh`/`checkMesh`; require positive cells, two solution directions and correct patches. Do not invoke either solver. Commit the meshing deliverable.

### Task 3: Actual interface states and full staged execution

**Files:** create `planar.py`; modify the planar dispatch points in `workflow.py` and the mesh preparation branch in `Allrun`; extend the same focused test file.

**Interfaces:**
- `planar.interface(case_dir: Path, time: str | None = None, allow_unsuitable: bool = False) -> dict` consumes saved breakdown arrays, generates fit/mesh inputs and saves fit/geometry reports.
- `planar.dsmc_setup(case_dir: Path, np_ranks: int | None = None) -> dict` writes the curved mesh input and existing DSMC dictionaries using its own geometry estimates.
- `planar.dsmc_fields(case_dir: Path) -> dict` reads the built curved mesh and writes the actual face-ordered inflow fields and suitability report.
- `sample_containing_cells(mesh: MeshGeometry, positions: ndarray, allowed: ndarray) -> ndarray` returns CFD cell labels and refuses points outside the allowed cell footprints.

- [x] Write and run failing mapping assertions with a small synthetic planar CFD field: each face samples its containing valid cell, positive states and correct inward DSMC normals are preserved, invalid/outside samples fail, and covered plus uncovered slit lengths equal D.
- [x] Implement bounded cell reconstruction using the containing accepted CFD cell (zero slope for this prototype). Save those cell labels and unit weights; do not use nearest-centre guesses or loosen the 3D coincidence checks. Derive number density and pressure consistently.
- [x] Evaluate actual curved faces and facets against the fit and accepted footprints. Compute normal Mach, forward/reverse half-range and NSF flux diagnostics, source coverage, closure, invalid states and outward overlap availability. Apply the existing configured suitability limits and explicit smoke-only override. Final acceptance happens after the actual DSMC mesh exists.
- [x] Derive DSMC parcel weight, timestep and schedule using its own mesh volumes and in-plane tracking sizes with existing collision/transport policy. Reuse existing DSMC dictionaries, near-vacuum initialization and field writers. No CFD-to-DSMC cell identity map is created.
- [x] Branch the staged runner so planar `dsmc` builds the curved mesh, prepares/assesses fields and only then initializes/runs DSMC; preserve the original 3D `subsetMesh` path. An infeasible fit or unapproved unsuitable interface must stop before solver launch. Run the focused mapping assertions and relevant existing hybrid interface/workflow tests. Commit the integrated workflow.

### Task 4: Manual assessment output and handoff

**Files:** modify `planar.py`, planar dispatch in `workflow.py`/`report.py`, `Allrun` and `cases/cai2012-hybrid/README.md`.

**Interfaces:**
- `planar.dsmc_post(case_dir: Path, time: str | None = None) -> dict` reads moments on the actual DSMC mesh and exports a basic planar profile/field summary.
- `planar.plot_interface(case_dir: Path) -> list[Path]` draws the raw accepted-region boundary and fitted ellipse over the breakdown field.
- Planar `reference` records `not_applicable` with the circular-disk model reason; it never evaluates that model.

- [x] Add a focused assertion that planar reference/post/report dispatch never requires `dsmc_to_cfd_cells.npy` or calls the disk analytical reference.
- [x] Implement basic planar post-processing, source-depth labels, raw/fitted interface overlay and `.foam` entry points. Reuse existing moment readers; skip 3D azimuthal and disk-reference comparisons.
- [x] Document exact serial commands for setup, CFD, fit/interface, DSMC and post-processing, plus optional rank/allow-unsuitable flags and the fit clearance setting. Clearly label the prototype as planar and under-resolved.
- [x] Run only the new focused test file and the existing hybrid tests affected by shared changes; do not repeat passing native mesh checks unless the mesh implementation changes. Run `git diff --check` and verify the generated laptop inputs without launching solvers.
- [x] Request one whole-change review using the native execution workflow, address actionable findings, and provide the ready-to-run case path, commands, checks performed and the statement that solvers were left for the owner. Commit the final deliverable.

## Plan Self-Review

The four tasks cover the complete approved data path: original breakdown logic,
conservative fit, genuinely planar curved mesh, actual-face transfer, unchanged
suitability limits and manually runnable solvers. The five review concerns are
pinned to the focused checks above. The plan intentionally has no solver runs,
new dependencies, parameter studies or general refactor.

## Execution record (2026-11-08)

Completed by native execution in the attached managed worktree.

- Final focused check: 67 tests passed across the planar, existing workflow,
  interface and source-physics test files. The full repository suite was
  intentionally omitted under the owner's minimal-testing instruction.
- Native checkMesh: the prepared 1,600-cell CFD mesh and one 3,960-cell synthetic
  curved DSMC mesh both reported Mesh OK. Checked matching depth planes,
  positive volumes, interface orientation and exact total slit area.
- Synthetic actual-face transfer: 120 faces received positive inlet states;
  containing-cell provenance, source coverage, unsuitable refusal and an explicit
  smoke-only override were checked. Native Allrun prepare-dsmc also passed.
- One final independent review found stale same-time breakdown provenance and
  a zero-clearance outer-contact case. Both were reproduced and fixed. Breakdown
  now binds field/configuration/mesh hashes; a quarter-smallest-cell minimum
  outer gap keeps zero-clearance fits separate from far boundaries.
- No CFD/DSMC solver, benchmark, convergence study or parameter sweep was run.
  Solver stability, physical accuracy, runtime and native mesher behavior beyond
  the synthetic fixture remain for the owner's manual assessment.

Implementation choices: retain the managed worktree to preserve unrelated local
work; use an equivalent log-area objective to avoid numerical collapse (the
result is still a local optimum); retain a positive outer-domain gap even when
excluded-cell clearance is zero. These may leave a smaller continuum region
than an ideal global fit. Run the supplied commands from this worktree.

Prepared case: cases/cai2012-hybrid/runs/planar-laptop-ellipse.
Run instructions: cases/cai2012-hybrid/README.md, planar prototype section.

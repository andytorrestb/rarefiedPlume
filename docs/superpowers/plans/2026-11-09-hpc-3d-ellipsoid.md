# HPC 3D Ellipsoid Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a fully 3D, Kn = 0.001 fitted-ellipsoid plume whose 192-core PBS job runs through DSMC sampling and post-processing.

**Architecture:** Preserve the existing planar and subsetMesh paths. Share the existing fitted-interface sampling/provenance/assessment helpers, add a 3D spheroid fit and native snappy mesh path, and expose one new HPC profile and launcher. The DSMC mesh is independent of the Cartesian CFD mesh.

**Tech Stack:** Existing Python/NumPy/SciPy/PyYAML, optional Matplotlib, OpenFOAM v2512 and the existing two compiled libraries; PBS/OpenMPI. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-11-09-hpc-3d-ellipsoid-design.md` (approved).

## Global Constraints

- Exit Kn = 0.001; D = 0.2 m, T = 300 K and speed ratio 2; unchanged argon VHS physics.
- Domain 8 D downstream and ±5 D in both transverse directions; full 3D, equal transverse fit radii, no wedge or planar-field rotation.
- Initial CFD mesh: D/40 core cells, x core 2 D, transverse core half-width D, 80 outer cells; 9,216,000 predicted cells.
- Keep Minmod and viscous timestep factor 0.25; existing coupling thresholds remain active.
- 192 MPI ranks on a 1.5 TB node; initial explicit PBS memory request 1200gb, with units and overhead stated.
- Default submitted workflow includes dsmcInitialise, parallel dsmcFoam startup/sampling and post-dsmc.
- No local solvers, production mesh, full test suite, benchmark, convergence study or parameter sweep.
- One small native 3D mesh/transfer check; only directly affected numerical/runner checks and one final review.

## Review Focus

- An excluded off-axis voxel or interior pocket must invalidate any enclosing fit, without azimuthal averaging (Task 1).
- Missing or changed snapshot provenance must stop fitting/transfer, including same-time CFD replacement (Task 1).
- Snapped facets outside the accepted region must stop before DSMC initialization (Task 2).
- Parallel mesher processor directories must never masquerade as a DSMC restart, and excessive memory estimates must stop execution (Task 2).
- The default PBS path must reach DSMC solve/post-processing with 192 ranks while preserving unsuitable-surface refusal and explicit smoke-only override (Task 3).

## File Structure

- New `plumetools/hybridplume/ellipsoid.py`: 3D conservative fit and voxel verification.
- New `plumetools/hybridplume/fitted.py`: shared fitted-interface sampling, provenance and face assessment extracted from `planar.py`, with backward-compatible planar aliases.
- New `plumetools/hybridplume/ellipsoid_mesh.py`: native background/snappy dictionaries and mesh metadata.
- New `plumetools/hybridplume/spatial.py`: 3D fitted workflow, actual-face mapping, own-mesh post-processing and meridional plotting.
- Modify `config.py`, `workflow.py`, `resources.py`, `report.py` in that package for explicit new-mode dispatch and memory guarding.
- Modify `cases/cai2012-hybrid/{hybrid.yaml,Allrun,README.md}` and add `run3DCase.sh`.
- New focused tests `tests/unit/test_hybridplume_ellipsoid.py`; reuse existing hybrid/planar helpers. No tests unrelated to changed behavior.

### Task 1: Full-3D fit and shared fitted-interface numerics

**Files:** `ellipsoid.py`, `fitted.py`, existing `planar.py`, `config.py`, `workflow.py`, `hybrid.yaml`, and the focused test file.

**Interfaces:**
- `fit_ellipsoid(lower: ndarray, upper: ndarray, retained: ndarray, *, source_radius: float, clearance: float) -> dict`: 3D Cartesian voxel bounds in metres; return status, xc/a/b, clipped volume, source intersection, clearances and optimizer/verification diagnostics.
- `verify_ellipsoid(fit: dict, lower: ndarray, upper: ndarray, retained: ndarray) -> dict`: bounded-memory verification against every excluded voxel.
- `fitted.sample_containing_cells(mesh, positions, allowed, axes=(0, 1)) -> ndarray`: explicit active axes; the spatial path passes `(0, 1, 2)`.
- Shared `breakdown_identity`, snapshot verification and face-assessment helpers retain their current contracts; planar exports remain compatible.
- New mode `interface.mode: fitted-ellipsoid` with `mesh.geometry: cartesian`; profile `hpc-kn1em3-ellipsoid` carries the Global Constraints.

- [ ] Write focused failing fit/sampling checks. Include a known excluded axis-crossing pocket whose first x bound is 0.6 m and clearance 0.01 m:

  ```python
  assert fit["xc"] + fit["a"] <= 0.59 + 1e-6
  assert containing_cells.tolist() == [0, 7]  # hand-labelled 2x2x2 CFD fixture
  ```

  Also cover an excluded off-axis voxel and refusal after a signed CFD snapshot changes. Reuse the existing roundoff/provenance fixtures where practical.
- [ ] Run only those new checks and confirm failure comes from the missing 3D behavior.
- [ ] Implement log clipped-volume optimization, with volume `pi*a*b*b*(2/3+s-s**3/3)`, `s=xc/a`. Reduce rejected-voxel constraints to the minimum rejected radial footprint in each x slab using all azimuths; independently verify the result against all 3D voxels in chunks. Preserve source intersection, excluded-cell clearance and a positive outer-domain gap; report infeasibility without fallback.
- [ ] Extract only the numerics/provenance/assessment shared by the two remeshed paths. Extend the containing-cell locator to explicit 3D axes, retaining coordinate clustering and invalid/outside refusal. Apply snapshot hashes before/after calculation to the new mode; leave the original subsetMesh path unchanged.
- [ ] Add and validate the new profile/mode. Keep direct-injection source handling compatible with absent or zero-face source patches when the fitted core covers the whole source.
- [ ] Run the new numerical checks and directly affected planar/configuration checks once; require passing results, then commit.

### Task 2: Native ellipsoid mesh, transfer and resource gate

**Files:** `ellipsoid_mesh.py`, `spatial.py`, shared dispatch in `workflow.py`, memory helper in `resources.py`, `Allrun`, and the focused tests.

**Interfaces:**
- `ellipsoid_mesh.write_mesh_inputs(case_dir: Path, cfg, source, fit: dict) -> list[Path]`: background blockMesh, searchableSphere snappy dictionary, source-patch creation and bounded meshing controls.
- `spatial.interface(case_dir, time=None, allow_unsuitable=False) -> dict` and `spatial.dsmc_setup(case_dir, np_ranks=None) -> dict`: signed fit/mesh inputs and provisional setup.
- `spatial.dsmc_fields(case_dir) -> dict`: actual mesh certification, face-ordered canonical inflow, suitability and final DSMC schedule/fields.
- `resources.peak_memory_estimate(cfg, cfd_cells: int, dsmc_estimate: dict, *, ranks: int) -> dict`: estimates both meshes, analysis temporaries, MPI and parcels; report GiB and refuse estimates over the configured allowance before initialization.

- [ ] Write failing checks for a snapped-facet enclosure entering an excluded voxel, a deliberately over-budget parcel estimate, and the meshing/restart distinction. They must demonstrate refusal before any DSMC initialization.
- [ ] Render a genuine 3D `(a b b)` searchableSphere cavity with configured interface/source refinement. Reuse the established snappy conventions without changing `cases/3d-inflow`. Keep global/local cell bounds and resolution controls explicit; no automatic parameter search.
- [ ] Add the separate-mesh branch to prepare-dsmc. Perform parallel meshing in a dedicated mesh-build directory, reconstruct before transfer, and keep its processor folders separate from the actual DSMC case. Reuse a matching completed mesh on retry; reject fit/mesh/cloud identity changes. Verify native flag spelling from installed tool help as needed.
- [ ] Locate actual interface faces in accepted 3D CFD cells, use bounded cell reconstruction and shared flux/suitability diagnostics, and preserve per-face provenance. Certify whole facets through a verified analytic enclosure of their vertices or an equivalent whole-facet check; centres alone are insufficient. Preserve uncovered source injection and measure disk/patch area error with an explicit profile tolerance and no density renormalization.
- [ ] Compute final particle weight, timestep, collision/resolution statistics and memory estimates from actual DSMC cells. Include both meshes, configurable analysis/MPI overhead and a safety allowance in the peak guard. Record assumptions; never alter physical density or silently raise particle weight merely to fit RAM.
- [ ] Run one small synthetic 3D native mesh/transfer fixture, using a small MPI rank count if needed to exercise the meshing lifecycle. Require checkMesh success, positive cells/states, contained facets, correct normals, closure and source coverage. No production mesh or solvers. Run only the associated new numerical checks; commit after they pass.

### Task 3: Complete PBS run, own-mesh reports and handoff

**Files:** new `run3DCase.sh`, `Allrun`, `spatial.py`, `workflow.py`, `report.py`, `README.md` and focused runner/report tests.

**Interfaces:**
- `spatial.dsmc_post(case_dir, time=None) -> dict`: existing moment readers, log/flux accounting and centreline statistics on the actual DSMC geometry; no CFD-cell identity map.
- `spatial.plots(case_dir) -> list[Path]`: thin meridional CFD/DSMC sections and fitted interface; retain .foam and VTK outputs.
- PBS launcher defaults: `PROFILE=hpc-kn1em3-ellipsoid`, 192 ranks, longq, 120-hour ceiling, `mem=1200gb`; resolves paths from the submitted checkout and invokes `Allrun ... all`.

- [ ] Add a stubbed runner check proving the launch path includes initialization, `dsmcFoam` and post-dsmc with the intended rank/profile arguments. A failed assessment must stop that path; explicit smoke permission must remain recorded. Add one synthetic own-mesh post-processing check that requires no CFD addressing map.
- [ ] Implement the dedicated job wrapper using the existing cluster environment and library-build steps. Support submission from the repository root or case script directory, avoid set-u, guard OpenFOAM initialization correctly, and keep existing run2DCase.sh/runCase.sh choices intact. Preserve resume/checkpoint behavior and explicit smoke-only authorization.
- [ ] Dispatch post-processing/plots for the new remeshed 3D path. Reuse moment/flux readers, use actual 3D coordinates, and label circular-source reference comparisons as collisionless-limit comparisons for this collisional case. Do not treat the new mesh as a CFD subset.
- [ ] Document exact qsub and restart commands, the memory units/controls, local refinement settings, actual-fit diagnostics and the fact that this is an initial HPC case without convergence validation.
- [ ] Run shell syntax, the stubbed workflow checks, the new focused tests and only existing hybrid tests affected by shared edits. Do not repeat a passing native mesh check unless mesh code changes. Run git diff --check.
- [ ] Commit the deliverable, request one focused whole-change review, address important findings with reproducing checks, and deliver the ready-to-submit script/profile. Do not submit remotely or launch local solvers.

## Self-review

The three tasks cover all approved stages and the five review concerns. Shared
interfaces use metre-based fit parameters and actual DSMC geometry throughout.
Native meshing stays separate from checkpoint ownership. Tests are bounded to the
new numerical/data-flow behavior and one small native fixture, as requested.

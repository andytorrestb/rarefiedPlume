# Portable DSMC Refinement and Cai Radial Study Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a portable utility that turns measured DSMC fields into mesh and parcel corrections, then use it in a seven-Kn radial study with fixed inlets matched to the Cai exit reference.

**Architecture:** `plumetools.refinement` owns solver-independent cell assessment, target plans, provenance, OpenFOAM readers, and a standalone hex refinement backend. `plumetools.radialstudy` supplies the Cai source mapping, expanding-streamtube mesh backend, case lifecycle, and study reports. Neither the shared decisions nor its CLI require a nozzle diameter, exit Kn, source-flow model, or study tree.

**Tech Stack:** Python >= 3.9; existing numpy, scipy, PyYAML and optional matplotlib; pytest; bash; OpenFOAM v2512; the owner's hyStrath v2512 `dsmcFoam+` / `dsmcInitialise+`. No new runtime dependency or C++ change.

**Spec:** [Approved design](../specs/2026-11-07-cai-exit-radial-resolution-design.md), including its portability requirement. Read both documents before implementation.

## Global Constraints

- Physical reference: `Kn_e = lambda_e / D_e`; argon at `[100, 10, 1, 0.1, 0.01, 0.001, 0.0001]` with the existing Cai exit inputs.
- Final spatial requirement: `Delta_edge_max(i) <= 0.8 lambda_measured(i)`; occupancy: `P_mean >= 20`. Adequate measured coverage is required. There is no mandatory `0.8 -> 0.4 -> 0.2` sweep.
- Predicted MFP sizes the initial case; sampled local density and temperature determine measured MFP through existing VHS code. An initially empty cell supplies no physical MFP.
- Initial radial target: `min(1.0 lambda_predicted(r), 0.05 D_e)`; fixed default inlet `r_a`; outer radius `25 D_e`; angle fixed before numerical revisions.
- Default startup is two transits, followed by six sampling transits divided into six blocks. Preserve physical domain, inlet, gas and sampling definition across revisions.
- Preserve all existing cases and results, `plumetools/cai2012`, `plumetools/markelov1999`, and `/home/andy/hyStrath-port`. Import their physics; do not copy formulas or edit frozen cases.
- Portable readers initially support reconstructed static single-species ASCII OpenFOAM cases. Explain unsupported input without modifying it. The pure API accepts normalized arrays from other readers.
- Generic mesh materialization writes a fresh artifact, never a parent mesh or parcel cloud. Inherited planning fields are not child measurements; a new flow run is required.
- Python code runs on Windows without OpenFOAM for unit tests; tool execution runs under Linux/WSL. Use `from __future__ import annotations`, no `match`, no `dataclass(slots=True)`, and no new runtime dependency.
- Use argument lists for subprocesses, preserve paths containing spaces, LF for new text files, and optional plot dependencies with useful CSV/JSON output when plots are unavailable.
- Do not run the scientific matrix or submit scheduler jobs. Only bounded mechanics checks are part of implementation verification.

## Review Focus

1. A case outside this repository with renamed fields and no Cai metadata must be assessable and locally refinable. Tests: Task 2 `test_custom_fields_outside_study`; Task 4 `test_refine_unrelated_cartesian_case`.
2. Empty cells, hyStrath's `0.001` placeholder, non-finite fields, and unknown/mismatched averaging windows must not yield a passing mesh or a guessed mesh correction. Tests: Tasks 1 and 2 `test_unreliable_measurements_remain_unresolved`, `test_sampling_window_required_and_consistent`.
3. Same-sized meshes or fields changed after planning must invalidate a stored plan; moving an unchanged case must not. Tests: Task 3 `test_changed_source_invalidates_plan`, `test_plan_is_relocatable`.
4. Binary, decomposed-only, dynamic or unsupported-topology cases must receive an actionable result, without touching the source. Tests: Tasks 2 and 4 `test_unsupported_input_is_read_only`, `test_failed_backend_keeps_source_and_destination_intact`.
5. A restart must neither reinitialize parcels nor mix cumulative/overlapping sampling windows into independent blocks. Tests: Task 8 `test_resume_preserves_cloud_and_sampling`, `test_overlapping_windows_are_rejected`.

---

## File Map

| Area | Files and responsibility |
|---|---|
| Portable data and decisions | `plumetools/refinement/{__init__,model,assess}.py`: normalized inputs, targets, per-cell statuses, MFP/occupancy assessment |
| Portable OpenFOAM input | `plumetools/refinement/{config,geometry,foam}.py`: standalone settings, actual mesh edges, standard/direct-field adapters |
| Portable plans and output | `plumetools/refinement/{plan,report,__main__}.py`: immutable plans, hashes, JSON/CSV, CLI |
| Generic mesh backend | `plumetools/refinement/openfoam.py`: isolated v2512 hex refinement and inherited target-field propagation |
| Study physics | `plumetools/radialstudy/{__init__,config,physics}.py`: seven-Kn settings and wrappers around existing Cai/interface physics |
| Radial mesh backend | `plumetools/radialstudy/mesh.py`: fixed streamtube, conservative graded intervals, measured-target remeshing |
| Study generation | `plumetools/radialstudy/{case,generation,__main__}.py`: plus dictionaries, immutable revisions, family CLI |
| Lifecycle and measurements | `plumetools/radialstudy/run.py`: environment, meshing, initialization, restart, disjoint sampling records |
| Study analysis | `plumetools/radialstudy/{post,report}.py`: native profiles, optional differences/plots, breakdown and required corrections |
| Family | `cases/cai2012-radial-resolution/{study.yaml,README.md,generate_cases.py,refine_case.py,analyse.py,baseCase/Allmesh,baseCase/Allrun,baseCase/Allpost}` |
| Tests | `tests/unit/_refinement.py`, focused `test_refinement_*` / `test_radialstudy_*`; `tests/openfoam/test_refinement_openfoam.py`, `test_radialstudy_openfoam.py` |
| Documentation | Portable guide `docs/mesh-refinement.md`; radial family guide; API indexes; existing CLI/configuration/output/catalog/nav registration |

These are responsibility boundaries, not permission for empty scaffolding. Add each file with its owning task's working behavior.

## Contracts Used Across Tasks

Define these dataclasses in `refinement/model.py` in Task 1. Arrays have one entry per cell, in mesh order; units are SI.

- `Geometry(centres_m, volumes_m3, max_edge_m, mesh_sha256)`.
- `Measurements(number_density, temperature_K, mean_parcels, window, source_sha256, provenance)`. `window` records start/end, sample basis, and whether its metadata was verified or explicitly declared.
- `Targets(edge_over_mfp=0.8, min_mean_parcels=20.0)`.
- `Assessment(targets, lambda_m, target_edge_m, mesh_status, occupancy_status, valid, summary)`. Status arrays use `-1` unavailable/untrusted, `0` failing, `1` passing. Summary action is `criteria_satisfied`, `needs_mesh`, `needs_parcels`, or `unresolved`.
- `RefinementPlan(action, mesh_sha256, source_sha256, source_geometry, target_edge_m, density_estimate, selected_cells, target_parcels, recommended_weight, provenance)`. `source_geometry` is the source `Geometry`; unknown global weight stays unavailable.
- `RefinementResult(destination, backend, mesh_sha256, n_cells, target_status, requires_flow_rerun, logs)`. `destination` is optional for a no-op; a no-op requires no new flow run.
- `MeshRefiner` protocol: `refine(source: Path, plan: RefinementPlan, destination: Path) -> RefinementResult`.

Window/provenance records are JSON-compatible dictionaries with named fields, not arbitrary Python objects. A stored plan is schema-versioned JSON plus `targets.npz` (loaded with `allow_pickle=False`) and a human-readable per-cell CSV. File references inside it are relative to the source case or artifact root, not the checkout location.

### Task 1: Portable Cell Assessment

**Files:** Create `plumetools/refinement/{__init__,model,assess}.py`; tests `tests/unit/{_refinement,test_refinement_assess}.py`.

**Interfaces:** Produce the contracts above and `assess(geometry: Geometry, measured: Measurements, species, *, convention: str = "vhs", targets: Targets = Targets()) -> Assessment`. `species` follows the existing MFP property protocol; no Cai configuration is accepted.

- [ ] **1. Write failing tests.** Add a normalized-array fixture `unit_inputs(*, edges, lambdas, parcels)` returning geometry, measurements, and a hard-sphere test species. Use existing gas conversion helpers to construct the requested MFPs. Pin decisions independently of any study:
  ```python
  g, m, gas = unit_inputs(edges=[0.8, 0.81, 0.4], lambdas=[1, 1, 1], parcels=[20, 40, 20])
  result = assess(g, m, gas)
  assert result.mesh_status.tolist() == [1, 0, 1]
  assert result.occupancy_status.tolist() == [1, 1, 1]
  assert result.summary["action"] == "needs_mesh"
  ```
  Also assert that a fully passing case requests no finer target; zero/NaN density, invalid temperature, unknown window, and the `0.001` empty placeholder leave affected mesh statuses unresolved. Counts below 20 request population/sampling work rather than trusted mesh refinement.
- [ ] **2. Run red:** `python -m pytest tests/unit/test_refinement_assess.py -q`; expect missing module/API failure.
- [ ] **3. Implement the contracts and assessment.** Reuse `rarefaction.core.mean_free_path`; keep the study's `>=20` mask distinct from the existing generic `>20` mask. Document a floating-point comparison tolerance of `1e-12` relative to the edge target, not a tunable relaxation. Do not silently mask uncovered cells out of the overall verdict.
- [ ] **4. Run green:** the same module passes, including exact-threshold and non-finite cases; run existing `tests/unit/test_rarefaction.py` to establish no shared-physics regression.
- [ ] **5. Commit:** only the portable contracts, assessment, and their tests; message `feat: assess DSMC mesh and occupancy from measured cells`.

### Task 2: Standalone Settings and OpenFOAM Readers

**Files:** Create `plumetools/refinement/{config,geometry,foam}.py`; extend `tests/unit/_refinement.py`; tests `tests/unit/{test_refinement_config,test_refinement_foam,test_refinement_geometry}.py`.

**Interfaces:** Produce `RefinementConfig`, `load_settings(path: Path) -> RefinementConfig`, `read_geometry(case: Path, *, time: str) -> Geometry`, and `read_snapshot(case: Path, config: RefinementConfig, *, time: Optional[str] = None) -> Snapshot`. Define `Snapshot(geometry, measured, species, convention, particle_weight, time, provenance)` in `foam.py`; it is an I/O wrapper, not required by Task 1's API.

- [ ] **1. Write failing tests.** Define `tree_hash(root: Path) -> dict[str, str]` (relative file paths to streamed SHA-256) in the test helper. Its `write_cartesian_case(root, *, profile="custom", lambdas=(0.75, 2.0), parcels=(40, 40)) -> Path` creates a valid two-cell ASCII mesh, `C`, `V`, direct fields with custom names, and sampling metadata. Pin:
  ```python
  case = write_cartesian_case(tmp_path / "unrelated case")
  snap = read_snapshot(case, load_settings(case / "refinement.yaml"))
  assert snap.geometry.max_edge_m.tolist() == [1.0, 1.0]
  assert assess(snap.geometry, snap.measured, snap.species).mesh_status.tolist() == [0, 1]
  assert not (case / "case.yaml").exists()
  ```
  Cover standard averaged moments, direct hyStrath fields, custom species suffixes, uniform fields, unequal cell volumes, unknown keys, absent/mismatched window metadata, invalid field dimensions/counts, and binary/decomposed/dynamic input explanations. On a rotated stretched hex assert actual maximum edge, rather than face diagonal, bounding-box extent, or `V^(1/3)`.
- [ ] **2. Run red:** `python -m pytest tests/unit/test_refinement_config.py tests/unit/test_refinement_foam.py tests/unit/test_refinement_geometry.py -q`.
- [ ] **3. Implement readers.** Settings have `gas` (the existing generic molecular mass/diameter/omega/reference-temperature/convention keys), `fields`, `sampling`, and `targets`; reject unknown keys. `fields.profile` is `standard`, `hystrath`, or `custom`; role overrides are `number_density`, `temperature`, `mean_parcels`, `mass_density`, `momentum_density`, `kinetic_energy_density`, `cell_centres`, and `cell_volumes`. Standard defaults are `rhoNMean`, `rhoMMean`, `momentumMean`, `linearKEMean`, and `dsmcRhoNMean`; hyStrath defaults are `rhoN_<suffix>`, `Ttra_<suffix>`, `dsmcNMean_<suffix>` with suffix defaulting to the gas species label. Use existing `foamio.volfield`, polyMesh CSR readers, and `rarefaction.foam.read_window` / `core.translational_temperature`. Direct fields supply direct temperature. Traverse face-adjacent vertex pairs to accumulate each owner's/neighbour's longest edge. Missing `C`/`V` explains the explicit `postProcess` commands needed, rather than modifying the case. Detect unsupported formats before parsing payloads. Sampling metadata for unrelated direct-field cases may be explicitly declared in settings and labeled as such. Stream hashes of mesh connectivity/points/boundary plus existing refinement-level/history files, and source fields, independently of absolute case location.
- [ ] **4. Run green:** the three modules and existing `test_foamio_shared.py`, `test_mesh_parsing.py`, `test_rarefaction.py` pass. Assert every read-only test leaves its complete source tree unchanged.
- [ ] **5. Commit:** `feat: read portable DSMC refinement inputs from OpenFOAM cases`.

### Task 3: Plans, Target Export, and Portable CLI

**Files:** Create `plumetools/refinement/{plan,report,__main__}.py`; tests `tests/unit/{test_refinement_plan,test_refinement_cli}.py`.

**Interfaces:** Produce `plan_refinement(geometry, measured, assessment, *, particle_weight: Optional[float] = None) -> RefinementPlan`, `write_plan(plan, directory: Path) -> Path`, `read_plan(path: Path) -> RefinementPlan`, `verify_sources(case: Path, plan) -> None`, `write_assessment(geometry, measured, assessment, directory: Path) -> Path`, and `main(argv=None) -> int`.

- [ ] **1. Write failing tests.** Assert:
  ```python
  assert passing_plan.action == "criteria_satisfied"
  assert passing_plan.selected_cells.size == 0
  assert failing_plan.selected_cells.tolist() == [0]
  assert failing_plan.target_edge_m[0] == pytest.approx(0.6)  # lambda=0.75
  assert read_plan(write_plan(failing_plan, tmp_path)).mesh_sha256 == failing_plan.mesh_sha256
  ```
  Pin field/mesh changes with unchanged array sizes invalidating a plan, relocation preserving validity, unknown/variable particle weight remaining unavailable, invalid measurements preventing an invented target, and no output created by help or dry-run. A same-mesh parcel correction from mean count 10 to target 20 recommends at most half the known weight and remains an estimate.
- [ ] **2. Run red:** `python -m pytest tests/unit/test_refinement_plan.py tests/unit/test_refinement_cli.py -q`.
- [ ] **3. Implement portable plans and CLI.** Provide `assess CASE --settings FILE [--time T] [--output DIR]`, `plan CASE --settings FILE [--time T] --output DIR [--dry-run]`, and reserve `apply CASE --plan FILE --backend hex --output DIR` for Task 4. Read-only mode prints a concise summary; writing is explicit. Plans retain target arrays, measured-source identity/window, conservative density estimates, and source cell IDs/centres/volumes/edge lengths for other meshers. Assess all coverage; `needs_parcels` / unresolved coverage does not become an executable full mesh plan. Do not infer constant weight from a variable-weight case. Export JSON, NPZ, and CSV; no original time-directory fields are overwritten.
- [ ] **4. Run green:** both test modules pass; run CLI on custom-field fixture without a Cai reference or repository cwd.
- [ ] **5. Commit:** `feat: export measured refinement plans with provenance`.

### Task 4: Portable OpenFOAM Hex Backend

**Files:** Create `plumetools/refinement/openfoam.py`; finish CLI apply; tests `tests/unit/test_refinement_backend.py`, `tests/openfoam/test_refinement_openfoam.py`.

**Interfaces:** Produce `HexRefiner(max_stages: int = 20, max_cells: Optional[int] = None).refine(source: Path, plan: RefinementPlan, destination: Path) -> RefinementResult`; implement the Task 1 protocol. Export a ready-to-use new mesh artifact, not a runnable solver case.

- [ ] **1. Write failing unit and integration tests.** Stub command execution only for failure/no-write paths. The real integration test builds two adjacent unit cubes independently of the study and selects only the first (`lambda=0.75`, target edge 0.6):
  ```python
  result = HexRefiner().refine(source, plan, output)
  assert result.n_cells == 9
  assert result.requires_flow_rerun is True
  assert result.target_status == "inherited_targets_met"
  assert read_geometry(output, time="0").volumes_m3.sum() == pytest.approx(2.0)
  assert tree_hash(source) == original_hash
  ```
  Add a multi-pass case, conformity extension, stale plan refusal, unsupported topology failure, failed tool process, budget exhaustion, output-path collision, and a passing-plan no-op. Neither a mapped target nor a mapped density is labeled as a new measured result.
- [ ] **2. Run red:** unit module and `python -m pytest tests/openfoam/test_refinement_openfoam.py -m needs_openfoam -q` under the sourced v2512 environment. Integration fixtures check `refineHexMesh`, `blockMesh`, `postProcess`, and `checkMesh` directly; they must not require the unrelated standard DSMC solver/library fixture.
- [ ] **3. Implement isolated materialization.** Verify hashes before creating a destination. Stage a minimal fresh case containing a copy of the source mesh, controlled `controlDict`/finite-volume dictionaries with ASCII output, and scalar target-edge / inherited-density fields. Write a cellSet of current violations. Invoke `refineHexMesh <set> -case <stage> -overwrite` with its default conformity extension, never `-minSet` or global refinement. Its field mapping propagates the fixed inherited targets to children. Regenerate `C` and `V` with `postProcess` after each pass; mapped old centroids/volumes are not new geometry. Recompute actual edges, repeat only while targets are missed, enforce stage/resource limits, run `checkMesh`, and atomically publish the fresh artifact. Keep logs and provenance. Exclude source result histories/clouds but retain mesh refinement history needed by the mesher. A passing plan returns a no-op result with `destination=None`. On failure leave the source and pre-existing output untouched. All staging cleanup must resolve inside its owned temporary root.
- [ ] **4. Run green:** real v2512 local and multi-pass tests pass; inspect actual cell counts, boundaries, volume, inherited fields, and source hashes. This proves portability through actual mesh construction, rather than only a second reader fixture.
- [ ] **5. Commit:** `feat: materialize portable hex refinement in new mesh artifacts`.

### Task 5: Cai Radial Reference and Breakdown

**Files:** Create `plumetools/radialstudy/{__init__,config,physics}.py`, `cases/cai2012-radial-resolution/study.yaml`; tests `tests/unit/{test_radialstudy_config,test_radialstudy_physics}.py`.

**Interfaces:** Produce `load_study(path: Path) -> RadialStudyConfig`, `source_for_kn(config, kn: float) -> RadialModel`, and `predict_breakdown(model, *, definition: str = "general", threshold: float = 0.05) -> dict`. `RadialModel` holds `config`, the existing equivalent `source`, resolved domain, angle policy, and source/domain identity; expose its axis state as `model.state(r)`. Case physical identity later includes the once-resolved numeric angle from Task 6.

- [ ] **1. Write failing tests.** Check exact seven-Kn ordering, inherited D/T/S/gas, density's reciprocal Kn scaling, unchanged physical hashes under numerical changes, and:
  ```python
  assert model.source.r_a / model.source.diameter_m == pytest.approx(0.7172871)
  assert model.source.T_a == pytest.approx(299.845, abs=0.001)
  assert low_kn_breakdown["r_over_D"] == pytest.approx(14.992, abs=0.002)
  assert kn001_breakdown["r_over_D"] == pytest.approx(2.788, abs=0.002)
  assert kn01_breakdown["status"] == "ALREADY_NONCONTINUUM_AT_INLET"
  ```
  Pin a shorter domain returning `NOT_REACHED_IN_DOMAIN`, invalid states, radial/general definition distinction, incompatible reference configuration, and forbidden inlet movement in a numerical revision.
- [ ] **2. Run red:** `python -m pytest tests/unit/test_radialstudy_config.py tests/unit/test_radialstudy_physics.py -q`.
- [ ] **3. Implement wrappers.** Load Cai physics through its current config/inflow APIs, snapshot only its physical inputs, and call the existing equivalent-source/state/Bird functions. Adapt axis root-search semantics; a domain-limited interface result is not a breakdown crossing. Keep initial predicted ratio 1.0 separate from final measured ratio 0.8. Defaults, limits, and mapping assumptions are the spec's; no legacy numerical constants are copied.
- [ ] **4. Run green:** these modules plus existing `test_cai2012ext_interface.py` pass. Confirm the shared `refinement` package contains no radialstudy/Cai imports or case paths.
- [ ] **5. Commit:** `feat: define the seven-Kn radial reference study`.

### Task 6: Graded Streamtube Mesh and Radial Refinement Backend

**Files:** Create `plumetools/radialstudy/mesh.py`; tests `tests/unit/test_radialstudy_mesh.py`, radial mesh tests in `tests/openfoam/test_radialstudy_openfoam.py`.

**Interfaces:** Produce `initial_mesh(model: RadialModel) -> RadialMesh`, `render_block_mesh(mesh: RadialMesh) -> str`, `read_radial_mesh(artifact: Path) -> RadialMesh`, and `RadialRefiner.refine(source: Path, plan: RefinementPlan, destination: Path) -> RefinementResult`. `RadialMesh` supplies radial boundaries, transverse counts, fixed angle/domain identity, exact cell vertices/centres/volumes/longest edges, minimum conservative transport length, and density estimates in cell order. Estimates are tagged `predicted` initially or `parent-measured` on a revision; they do not change the physical inlet. Artifacts include enough radial mesh metadata for `read_radial_mesh` without recomputing the frozen angle.

- [ ] **1. Write failing tests.** Check closed positive volume against the analytic frustum volume, orientation, actual slanted/transverse edges, finite coordinate precision at Kn `0.0001`, geometric-cap behavior at Kn 100, and fixed angle/inlet/outlet across revisions. For a parent's local measured failure, assert every revised edge respects its inherited 0.8-MFP bound and that already-sufficient radial intervals are retained when feasible. Assert no fixed eightfold or factor-of-two global rule. A target plan with no violations requests no new mesh.
- [ ] **2. Run red:** unit mesh module and a small real `blockMesh` integration fixture.
- [ ] **3. Implement the radial backend.** The baseline follows `s_0(r)`, with conservative checks over intervals and the spec's once-resolved angle. Revised bounds come from the shared plan's minimum valid measured target among transverse cells per parent radial interval. Preserve practical existing radial boundaries and insert sufficient boundaries; choose shared integer transverse counts for conformity. Use existing formatting/readers. Record required cells versus optional budgets and any extra transverse resolution. Output a mesh artifact with the same `RefinementResult` contract as the generic backend.
- [ ] **4. Run green:** unit tests and a real built-mesh comparison of `C`, `V`, edge lengths, boundaries, and total volume against `RadialMesh`. The portable Task 4 tests still pass independently.
- [ ] **5. Commit:** `feat: build fixed-domain radial meshes from measured targets`.

### Task 7: Plus Dictionaries and Non-Destructive Case Revisions

**Files:** Create `plumetools/radialstudy/{case,generation,__main__}.py`, `cases/cai2012-radial-resolution/{generate_cases.py,refine_case.py,analyse.py,baseCase/Allmesh,baseCase/Allrun,baseCase/Allpost}`; tests `tests/unit/{test_radialstudy_case,test_radialstudy_generation}.py`.

**Interfaces:** Produce `write_case(destination: Path, model: RadialModel, mesh: RadialMesh, *, target_parcels: float = 20.0, particle_weight: Optional[float] = None, delta_t_s: Optional[float] = None, parent: Optional[Path] = None) -> dict`, `generate_initial(config, output: Path, *, kn_values=None, dry_run: bool = False) -> list`, and `revise_case(parent: Path, plan: RefinementPlan, *, correction: str = "needed", delta_t_s: Optional[float] = None, dry_run: bool = False) -> dict`. `write_case` derives weight from the mesh's tagged density estimates unless an explicit measured correction is supplied. `revise_case` consumes Task 6's artifact/read-back contract; a timestep-only correction retains mesh and weight. The study CLI exposes `generate`, `refine`, and later `analyse`; wrappers delegate without duplicating decisions.

- [ ] **1. Write failing tests.** Generate seven `r0_ppc20` cases, no 63-case grid. Check inlet n/U/T against the mapping, plus application names and boundary models, planned `W = min(n_i V_i)/20`, fixed seed provenance, and all physical input hashes. Repeat generation byte-for-byte, preserve sentinel solver output, refuse changed inputs/existing child targets, create a mesh-only correction at fixed requested occupancy and a population-only correction on the identical mesh. A timestep-only revision must retain the parent's mesh and weight. A satisfied assessment creates nothing. Relocate the family and preserve relative reference/manifest behavior.
- [ ] **2. Run red:** `python -m pytest tests/unit/test_radialstudy_case.py tests/unit/test_radialstudy_generation.py -q`.
- [ ] **3. Implement generation.** Use a new family model tag and staged atomic publication; follow existing conflict behavior without invoking the Cai-specific policy/schema. Write native plus dictionaries with uniform mapped inlet, deletion at inlet/outlet, supported specular sides, Cartesian weighting, VHS gas, `seedNumber`, static mesh, and direct `dsmcVolFields`. Create stable append-only revisions, independent mesh hashes, measured source/window hashes and parent identity. A child is initialized afresh by its family lifecycle; never map an old parcel cloud to a new mesh. Optional higher ppc or stricter target variants remain explicit user requests. Add narrow family-only ignore entries when generated outputs first exist.
- [ ] **4. Run green:** both modules; `generate_cases.py --dry-run` reports seven initial cases without creating files, and an explicitly temporary generation checks all seven dictionaries/weights. Existing case trees remain unchanged.
- [ ] **5. Commit:** `feat: generate immutable radial cases and measured corrections`.

### Task 8: Solver Lifecycle and Disjoint Sampling Windows

**Files:** Create `plumetools/radialstudy/run.py`; complete `cases/cai2012-radial-resolution/baseCase/{Allmesh,Allrun,Allpost}`; tests `tests/unit/{test_radialstudy_run,test_radialstudy_sampling}.py`; extend `tests/openfoam/test_radialstudy_openfoam.py`.

**Interfaces:** Produce `mesh_case(case: Path) -> dict`, `run_case(case: Path, *, resume: bool = False) -> dict`, `read_windows(case: Path) -> list`, and `collect_case(case: Path) -> list`. Each window record has actual start/end time, step count, mesh/physical identity, field names and hashes, and disjoint sampling basis. Keep command execution behind an injectable runner for unit failures.

- [ ] **1. Write failing tests.** Pin environment/command detection, source-before-plus setup examples, timestep from both limits, startup/sampling durations invariant across revisions, six blocks with two disjoint output subwindows per block, native reset-at-write behavior, source/field identities, and restart handling. A fresh run refuses existing output; resume does not initialize or discard the cloud; incomplete checkpoints are identified using `foamio.checkpoints` and the user receives the existing concrete recovery action. Overlapping/cumulative windows cannot be accepted as independent. Include path spaces and command failure logging. If MPI is exposed, assert rank-distinct explicit seeds; otherwise clearly scope the implemented runner to serial.
- [ ] **2. Run red:** the two unit modules with stub executables, then the plus integration test under WSL.
- [ ] **3. Implement lifecycle.** Use the actual plus port's dictionaries and checkpoint layout, existing checkpoint/time helpers, and explicit argv commands. Enable `resetAtOutput` continuously; do not set the old `resetAtOutputUntilTime` that would turn later snapshots cumulative. Default write cadence is half a sampling transit, two subwindows per block; record achieved times after step rounding. Persist/recover sampling state on resume and refuse a restart that would invalidate its window definition. No AMR or automatic variable weights. `Allmesh` builds/checks geometry; `Allrun` uses fresh or explicit resume semantics; `Allpost` collects fields and windows without launching a solver.
- [ ] **4. Run green:** unit modules plus a bounded actual plus mechanics case. Use Kn 100, outer radius `0.82 D_e`, mapped inlet, one core, <=100 cells, <=2000 steps and <=5000 estimated parcels. This is a labeled mechanics configuration, not a scientific baseline. Verify mesh, initialization, inlet/outlet/model loading, cloud/checkpoint writing, direct field names and disjoint windows. Record actual time/counts; failure must remain visible.
- [ ] **5. Commit:** `feat: run radial DSMC cases with auditable sampling and restart`.

### Task 9: Study Assessment, Profiles, and Required Corrections

**Files:** Create `plumetools/radialstudy/{post,report}.py`; complete `cases/cai2012-radial-resolution/analyse.py` / CLI analysis; tests `tests/unit/{test_radialstudy_post,test_radialstudy_report}.py`.

**Interfaces:** Produce `assess_case(case: Path) -> dict`, `analyse_study(root: Path, *, kn_values=None) -> dict`, and `write_study_report(root: Path, analysis: dict, *, plots: bool = True) -> list[Path]`. Case results reference the shared `Assessment` and `RefinementPlan` artifacts, plus independent timestep, sampling and breakdown evidence. A needed timestep correction includes `recommended_delta_t_s` for Task 7's override; it must not disappear merely because mesh and occupancy passed.

- [ ] **1. Write failing tests.** A changed measured MFP, with fixed inlet Kn, changes local required edge limits. All reliable cells passing 0.8 and mean occupancy 20 requests no additional case. Missing/low-signal windows remain unresolved. Check measured/predicted MFP side by side, required population correction after changing volume, no optional-comparison requirement for `criteria_satisfied`, and honest missing-run inventory. For unequal radial meshes conserve bin volumes and moments while projecting fine profiles onto the coarse intervals. Test drift detection, disjoint-block errors, optional metric verdicts, unavailable measured gradients, and CSV output without matplotlib.
- [ ] **2. Run red:** both post/report modules.
- [ ] **3. Implement assessment through the portable API.** Process windows one at a time; retain block profiles and necessary moment accumulators rather than every full field snapshot. Aggregate through conservative density/momentum/energy moments, reusing existing moment/temperature helpers rather than naively averaging temperatures across different bulk velocities. Build native radial layers; compare only identical physical/domain/sampling identities using existing normalized L2/uncertainty verdicts. Reuse existing radial/Bird and numerical-gradient helpers for predicted/optional measured breakdown. Timestep and sampling failures remain separate from mesh acceptance. Write case settings/plan/CSV/JSON, study CSV/Markdown, and optional overlays of n/U/T, both MFPs, edge/MFP with 0.8, ppc with 20, and time ratios. Print the specific remaining correction or `criteria_satisfied`; never launch it.
- [ ] **4. Run green:** unit modules and post-processing of Task 8's mechanics output. The latter proves field/format interoperability; do not require a tiny run to prove scientific stationarity or convergence. Existing overlap/convergence tests still pass.
- [ ] **5. Commit:** `feat: report measured radial sizing and remaining corrections`.

### Task 10: Documentation and Final Verification

**Files:** Create `docs/{mesh-refinement,cai2012-radial-resolution}.md`, `docs/tutorials/cai2012-radial-resolution.md`, `docs/reference/api/{refinement,radialstudy}/index.md`, family README; update `.gitignore`, `docs/_tools/hooks.py:FAMILIES`, `mkdocs.yml`, and CLI/configuration/output/case-catalog references. Extend `docs/_tools/run_examples.py` only with bounded dictionary/read-only examples.

**Interfaces:** Document exact shipped commands and settings; all examples use the actual APIs from previous tasks. No new product behavior belongs in this task.

- [ ] **1. Document two complete workflows.** The study path is predict -> run -> measure -> correct only failures -> rerun, with fixed inlet/domain and 0.8 measured-MFP completion. The unrelated-case path supplies standalone gas/field/window settings, assesses and exports a plan, materializes a new hex mesh artifact, and hands it to that case's initializer. Explain supported inputs/backends and how another mesher implements `MeshRefiner`; include custom fields without Cai metadata. Distinguish inherited planning fields from child measurements and streamtube totals from full-plume totals.
- [ ] **2. Run the final focused checks under WSL.** Run all new `test_refinement_*` and `test_radialstudy_*` unit modules, plus existing `test_rarefaction.py`, `test_foamio_shared.py`, `test_foamio_checkpoints.py`, `test_cai2012ext_interface.py`, and `test_cai2012ext_overlap.py`. Run only the two new bounded OpenFOAM integration modules with `-m needs_openfoam`. Confirm seven initial dictionary sets, idempotence, no-op for a compliant case, and each independent correction path. Use a known writable `--basetemp`; native Windows sandbox temp discovery previously failed and does not justify changing tests.
- [ ] **3. Verify portability and evidence artifacts.** Keep the unrelated Cartesian case's selected cells, actual refined mesh statistics, total volume, source hashes and logs in a temporary verification directory. Record the plus mechanics run's actual executable/version, fields/windows/counts and outcome. Check no old case, frozen physics module, or WSL port file changed. No scientific matrix or scheduler command is executed.
- [ ] **4. Run documentation checks:** `python docs/_tools/check_api.py`, `python docs/_tools/check_links.py`, `python docs/_tools/run_examples.py`, `python -m mkdocs build --strict`, and `git diff --check`. Use the existing docs environment; install/build work, if needed, remains part of the approved execution environment setup and must not alter physics or test requirements. Fix new failures; record unrelated baseline failures precisely if any.
- [ ] **5. Commit:** `docs: explain portable refinement and the Cai radial study`; finish with a whole-branch review and a concise evidence-based report. Report unrun scientific work explicitly rather than treating mechanics checks as validation.

## Self-Review Record

- Spec coverage: physics and fixed domain (Tasks 5-6), measured acceptance and population (1-3, 9), actual portable refinement (4), immutable cases (7), plus execution and sampling (8), profiles and optional convergence (9), docs and evidence (10).
- Portability is demonstrated twice: standalone custom-field ingestion outside the repository (2), and an actual refined unrelated Cartesian mesh (4). The radial backend consumes the same plan/protocol (6).
- Review Focus tests are assigned above; stale sources, no-op completion, missing coverage, source preservation, and restart windows have explicit owners.
- Interfaces and model names are shared in one contract section; only OpenFOAM I/O uses `Snapshot`. The pure assessment/target decision API needs normalized arrays, not a case family.
- The selected 0.8 measured-MFP target is a stopping condition. No task depends on a predetermined finer-mesh grid or optional profile-convergence evidence for a criteria pass.

## Implementation References

- Approved spec and the existing source modules cited there.
- Installed v2512 source inspected during planning:
  `/usr/lib/openfoam/openfoam2512/applications/utilities/mesh/advanced/refineHexMesh/refineHexMesh.C` (field mapping, consistency extension, mesh/history writes).
- OpenCFD primary source reference for the related mesh utility:
  <https://api.openfoam.com/2506/refineMesh_8C_source.html>. Use the installed
  v2512 implementation as the version-specific authority.
- The port's `dsmcField.C` documents reset-at-output behavior;
  `dsmcVolFields.C` supplies direct field/count/MFP conventions; `dsmcCloud.C`
  supplies `seedNumber` semantics. No port edit is planned.

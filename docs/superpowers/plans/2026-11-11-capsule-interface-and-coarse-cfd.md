# Fitted Capsule Interface on a Coarse, Stabilised 3D CFD — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `fitted-capsule` interface mode (cylinder + tangent spheroidal tip, conservatively fitted and flux-screened) and a coarse 30 D × 15 D CFD profile kept stable by a kinematic-viscosity cap.

**Architecture:** Part 1 is configuration plus a small transport-library extension and two breakdown-step checks. Part 2 introduces `shapes.py` (a shared shape interface; the existing spheroid wrapped bit-identically), a capsule fitter/screen, and STL-based snappy meshing, all routed through the existing `spatial.py` pipeline so every downstream utility is reused.

**Tech Stack:** Python 3.11 (numpy, scipy), pytest, OpenFOAM v2512 (C++ thermophysical templates, snappyHexMesh, surfaceFeatureExtract), bash (Allrun), PBS.

**Spec:** `docs/superpowers/specs/2026-11-11-capsule-interface-and-coarse-cfd-design.md`

## Global Constraints

- Branch `capsule-interface`. Before Task 1 the owner commits or stashes their uncommitted edits to `cases/cai2012-hybrid/Allrun`, `hybrid.yaml`, `run3DCase.sh`, `runSetup.sh`; do not commit those edits on their behalf.
- Run Python and tests with the project venv: `.venv/bin/python -m pytest <path> -q` from the repo root.
- Do not run CFD, DSMC, snappyHexMesh or any production mesh locally. No convergence or sensitivity study.
- Existing profiles and the `fitted-ellipsoid` / `fitted-ellipse` modes must behave identically; existing tests keep passing.
- New profile name: `hpc-kn1em3-capsule`. Values verbatim from the spec:
  `domain: {x_max_over_D: 30.0, half_width_over_D: 15.0}`;
  `mesh: {geometry: cartesian, core_cell_over_D: 0.05, core_x_over_D: 3.0, core_half_over_D: 1.5, outer_cells: 40, max_outer_expansion: 60.0}`;
  `cfd: {end_time_transits: 1.0, viscosity_cap: {enabled: true, factor: 1.0}}`;
  `interface: {mode: fitted-capsule}`.
- Viscosity cap: `nu_max = factor · viscous_dt_factor · h² / dt_visc` (= factor · ν_bg); `mu_eff = min(mu_tab(T), nu_max · p/(R T))`; `kappa_eff = kappa_tab · mu_eff/mu_tab`. Off by default.
- Stability acceptance: zero capped cells inside the retained region; retained volume fraction change ≤ 0.02 between the last two writes (both are errors when violated).
- Capsule screen defaults: `interface.fit.screen: {radius_samples: 12, length_samples: 24, axial_points_per_D: 40, azimuthal_points: 90}`. Selection = largest volume with screened return ≤ `interface.criteria.max_return_mass_fraction` and subsonic area ≤ `max_subsonic_normal_area_fraction`; otherwise lowest-return feasible with `screen.status = "no_passing_candidate"`.
- STL chord error ≤ 0.1 × `interface.fit.geometry_relative_tolerance`.
- Spheroid fit identity must hash the same `[xc, a, b]` list as today.

## Review Focus

1. **Capsule with L = 0** (excluded voxels right at the exit): must degrade to a pure tip (half-spheroid on the source), not divide by zero or emit an empty cylinder — Task 6 test `test_fit_handles_zero_length_cylinder`.
2. **Only one CFD write exists** when checking retained stability: must report `"previous_write": null` and skip, not crash — Task 3 test `test_retained_stability_skips_without_previous_write`.
3. **Rerun of `interface` on a case already meshed with the other shape** (ellipsoid → capsule): must refuse via the existing identity check, not silently reuse the DSMC mesh — Task 9 test `test_capsule_interface_refuses_existing_ellipsoid_mesh`.
4. **Screen samples falling outside the CFD domain or on excluded cells** (tip touching the clearance band): candidate is marked infeasible, never an exception that aborts the sweep — Task 7 test `test_screen_marks_unsampleable_candidate_infeasible`.
5. **`viscosity_cap.enabled: true` with a non-positive or non-finite factor**: config error at load time — Task 2 test `test_viscosity_cap_factor_must_be_positive`.

---

## Part 1 — Coarse CFD and stability

### Task 1: Profile and spatial-mode routing

**Files:**
- Modify: `plumetools/hybridplume/config.py` (`INTERFACE_MODES` ~L41; validation ~L136-146; add `is_spatial_fit`)
- Modify: `plumetools/hybridplume/workflow.py` (every `== "fitted-ellipsoid"` / `in ("fitted-ellipse", "fitted-ellipsoid")` test: ~L132, 274, 339, 364, 468, 553, 564, 1077)
- Modify: `cases/cai2012-hybrid/Allrun` (L184, L191, L196)
- Modify: `cases/cai2012-hybrid/hybrid.yaml` (new profile after `hpc-kn1em3-ellipsoid`)
- Test: `tests/unit/test_hybridplume_capsule_profile.py`

**Interfaces:**
- Produces: `config.SPATIAL_FIT_MODES = ("fitted-ellipsoid", "fitted-capsule")`; `config.is_spatial_fit(cfg) -> bool`; `config.is_fitted(cfg) -> bool` (ellipse or spatial). Allrun shell function `spatial_fit` (exit 0 when the mode is in `SPATIAL_FIT_MODES`).

- [ ] **Step 1: Write failing tests**

```python
def test_capsule_profile_mesh_matches_spec():
    cfg = load(STUDY, "hpc-kn1em3-capsule"); src = config.build_source(cfg)
    p = mesh.plan(cfg, src)
    assert 1.90e6 < p.n_cells < 2.02e6
    assert p.max_expansion_ratio < 1.11
    assert abs(mesh.predicted_source(p, src).area_error) < 0.01
    assert cfg.get("interface.mode") == "fitted-capsule"

def test_is_spatial_fit():
    assert config.is_spatial_fit(load(STUDY, "hpc-kn1em3-capsule"))
    assert config.is_spatial_fit(load(STUDY, "hpc-kn1em3-ellipsoid"))
    assert not config.is_spatial_fit(load(STUDY, "planar-laptop-kn0p01"))

def test_fitted_capsule_requires_cartesian():  # overlay mesh.geometry: planar -> ConfigError

def test_allrun_has_no_literal_ellipsoid_mode_tests():
    text = (ROOT/"cases/cai2012-hybrid/Allrun").read_text()
    assert "= fitted-ellipsoid ]" not in text and "spatial_fit()" in text
```

Use the same `load` helper the existing workflow tests use (see `tests/unit/test_hybridplume_workflow.py`).

- [ ] **Step 2: Run, expect FAIL** — `.venv/bin/python -m pytest tests/unit/test_hybridplume_capsule_profile.py -q` (unknown profile / missing function).
- [ ] **Step 3: Implement.** Add the mode, `SPATIAL_FIT_MODES`, `is_spatial_fit`, `is_fitted`; cartesian check generalised to `is_spatial_fit`. Replace each workflow mode test with the helper (no behavioural change for existing modes). In Allrun add `spatial_fit() { case "$(cfg_value interface.mode automatic)" in fitted-ellipsoid|fitted-capsule) return 0;; esac; return 1; }` and use it at the three sites. Add the profile with the Global Constraints values; inherit everything else from `hpc-kn1em3-ellipsoid` by repeating its `source`, `cfd`, `dsmc`, `resources` blocks (profiles overlay `base`, not each other); interface `fit` keeps `{clearance_cells: 0.25, geometry_relative_tolerance: 0.02, max_source_area_error: 0.01}`. Label: `"HPC 3D fitted capsule: Kn_e 0.001 (coarsest stable mesh, 30D x 15D)"`.
- [ ] **Step 4: Run new tests plus** `tests/unit/test_hybridplume_workflow.py tests/unit/test_hybridplume_ellipsoid.py` and `bash -n cases/cai2012-hybrid/Allrun` — all PASS.
- [ ] **Step 5: Commit** — `feat: add fitted-capsule mode routing and coarse 30D x 15D profile`.

### Task 2: Viscosity cap — Python side

**Files:**
- Modify: `plumetools/hybridplume/cfd.py` (`timing` ~L170; `render_thermophysical_properties` ~L100; `write_case` ~L298)
- Modify: `plumetools/hybridplume/config.py` (validation)
- Test: `tests/unit/test_hybridplume_viscosity_cap.py`

**Interfaces:**
- Produces: `timing(...)` gains key `"nu_max_m2_s": float | None` (None when disabled). `render_thermophysical_properties(source, table=None, nu_max=None) -> str`: when `nu_max` is set, `transport cappedTabulated;` and `nuMax <value>;` inside the `transport` subdict; otherwise byte-identical to today. `setup.json` (written by `workflow`) carries `nu_max_m2_s` via the timing dict.

- [ ] **Step 1: Failing tests**

```python
def test_nu_max_equals_background_kinematic_viscosity():
    cfg = load(STUDY, "hpc-kn1em3-capsule"); src = config.build_source(cfg)
    tm = cfd.timing(cfg, src, mesh.plan(cfg, src)); bg = cfd.background_state(cfg, src)
    assert tm["nu_max_m2_s"] == pytest.approx(src.gas.viscosity(bg["temperature_K"]) / bg["density_kg_m3"], rel=1e-12)

def test_thermo_dict_selects_capped_transport_only_when_enabled():
    # enabled: "transport       cappedTabulated;" and "nuMax" present
    # ellipsoid profile: output identical to render_thermophysical_properties(source) today, no "nuMax"

def test_viscosity_cap_factor_must_be_positive():   # factor 0, -1, nan -> ConfigError
```

- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement** per Global Constraints formula; validation in `config.validate` for `cfd.viscosity_cap.factor` (finite, > 0) when enabled.
- [ ] **Step 4: Run new tests + `tests/unit/test_hybridplume_workflow.py::test_cfd_case_uses_vhs_transport_and_caps_the_viscous_step`** — PASS.
- [ ] **Step 5: Commit** — `feat: write a kinematic-viscosity cap for the coarse CFD profile`.

### Task 3: Viscosity cap — C++ transport, and breakdown accountability

**Files:**
- Create: `applications/thermophysicalModels/hybridPlumeThermo/cappedTabulatedTransport.H`
- Modify: `applications/thermophysicalModels/hybridPlumeThermo/hybridPlumeThermos.C` (second `makeThermos` with `cappedTabulatedTransport`; update the header comment)
- Modify: `plumetools/hybridplume/breakdown.py` (two functions), `plumetools/hybridplume/workflow.py::breakdown` (~L264-349)
- Test: `tests/unit/test_hybridplume_viscosity_cap.py` (extend)

**Interfaces:**
- Consumes: `cfd.timing(...)["nu_max_m2_s"]` (Task 2).
- Produces: `bd.capped_cells(state: CfdState, gas: Gas, nu_max: float) -> np.ndarray[bool]` (`gas.viscosity(T) > nu_max·rho`); `bd.retained_change(previous: np.ndarray, current: np.ndarray, volumes: np.ndarray) -> float` (|V_prev − V_cur| / V_cur over retained volume). `breakdown.json` gains `"viscosity_cap": {"nu_max_m2_s", "capped_cells", "capped_in_retained"}` and `"retained_stability": {"previous_write": str|None, "volume_fraction_change": float|None, "limit": 0.02}`.

- [ ] **Step 1: C++ wrapper.** `template<class Thermo> class cappedTabulatedTransport : public tabulatedTransport<Thermo>` with member `scalar nuMax_` read from `dict.subDict("transport").get<scalar>("nuMax")`; `typeName()` returns `"cappedTabulated<" + Thermo::typeName() + '>'`; overrides `mu(p,T)`, `kappa(p,T)`, `alphah(p,T)` per Global Constraints (`this->rho(p,T)` from the equation of state); `write` adds `nuMax`. Mirror the constructor/clone/New set of `tabulatedTransport.H` (see `$FOAM_SRC/thermophysicalModels/specie/transport/tabulated/`).
- [ ] **Step 2: Verify compile** in the HPC env: `source cases/cai2012-kn-extension/hpc_env.sh && bash applications/thermophysicalModels/Allwmake` — builds `libhybridPlumeThermo.so` with no errors; `strings $FOAM_USER_LIBBIN/libhybridPlumeThermo.so | grep cappedTabulated` non-empty. If the env cannot be sourced on this node, record that and hand the compile to the owner's HPC session.
- [ ] **Step 3: Failing Python tests**

```python
def test_capped_cells_flags_only_low_density():   # two-cell state: rho_bg and rho_bg/100 at 300 K -> [False, True]
def test_breakdown_errors_when_capped_cell_is_retained(...)  # monkeypatch state so a retained cell is capped -> WorkflowError "viscosity cap active in retained continuum"
def test_retained_change_is_volume_weighted():
def test_retained_stability_errors_above_limit(...)  # change 0.03 -> WorkflowError
def test_retained_stability_skips_without_previous_write(...)  # one write -> previous_write None, no error
```

Reuse `_hybridmesh.box_mesh` for a tiny mesh, as `test_hybridplume_ellipsoid.py` does.
- [ ] **Step 4: Implement** the two functions and wire them into `workflow.breakdown`: accountability only when `cfd.viscosity_cap.enabled`; stability only when `breakdown.retained_stability_limit` is set (add `0.02` to the capsule profile in this task; absent → section omitted, so existing profiles are unchanged). The previous write is the predecessor of `time` among the numeric time dirs; classify it with the same thresholds/criterion. Write `viscosityCapped` into `fields` for ParaView.
- [ ] **Step 5: Run** new tests + `tests/unit/test_hybridplume_workflow.py` — PASS.
- [ ] **Step 6: Commit** — `feat: capped tabulated transport with retained-region accountability`.

**Owner checkpoint after Task 3:** submit CFD + breakdown for `hpc-kn1em3-capsule` (`STAGE` up to breakdown) so Part 2 can be fitted on real data. Part 2 does not depend on that run to be implemented.

---

## Part 2 — Fitted capsule interface

### Task 4: Shape interface with the spheroid wrapped

**Files:**
- Create: `plumetools/hybridplume/shapes.py`
- Test: `tests/unit/test_hybridplume_shapes.py`

**Interfaces:**
- Produces (both frozen dataclasses): `Spheroid(xc, a, b)`, `Capsule(R, L, a)` with
  `q(points: (n,3)) -> (n,)`, `deviation(vertices) -> float`, `outward_gradient(points) -> (n,3)`,
  `verify(lower, upper, retained, clearance: float) -> dict` (same keys as `verify_ellipsoid`),
  `grown(factor: float) -> Self`, `section(n: int) -> (m,2)` (x, y at z=0, upper and lower),
  `volume() -> float`, `base_radius() -> float`, `params() -> list[float]`,
  `to_json() -> dict` (includes `"shape": "spheroid"|"capsule"`), and module function
  `shapes.from_json(doc: dict) -> Spheroid | Capsule` (a doc without `"shape"` is a spheroid — today's `ellipsoid.json`).
- `Spheroid.verify` calls `ellipsoid.verify_ellipsoid` with `{"xc","a","b","clearance_m"}`; `Spheroid.params() == [xc, a, b]`.
- `Capsule.verify`: excluded boxes (padded by clearance) fail if they intersect the cylinder (`x`-overlap with [0, L] and min radial distance of the box's y-z rectangle to the axis < R) or the tip spheroid (`verify_ellipsoid` with centre L, axes (a, R, R), restricted to boxes with `hi_x > L`). `grown(g)` scales R and a, keeps L. `base_radius() == R`. Volume `πR²L + (2/3)πR²a`.

- [ ] **Step 1: Failing tests** — `q == 1` on sampled surface points of both shapes; `outward_gradient · (p − interior point) > 0`; capsule volume formula; `grown(1.01)` capsule fails `verify` on a voxel set it passes at 1.0; spheroid `verify` output equals `verify_ellipsoid` output on the `voxels()` fixture from `test_hybridplume_ellipsoid.py`; `from_json(Spheroid(...).to_json())` round-trips; a legacy `ellipsoid.json` dict (no `"shape"`) loads as `Spheroid`.
- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** — `feat: shared fitted-shape interface for spheroid and capsule`.

### Task 5: Route `spatial.py` and `fitted.py` through shapes (spheroid regression)

**Files:**
- Modify: `plumetools/hybridplume/spatial.py` (`certify_facets` L32-47, `interface` L50-81, `dsmc_setup` L84-130, `_geometry` L154-189, `dsmc_fields` L200+, `plots` L269+)
- Modify: `plumetools/hybridplume/fitted.py::_identity` (L60-65)
- Modify: `plumetools/hybridplume/ellipsoid_mesh.py` (factor `_background`, `_source_regions`, `_finish_dicts` helpers out of `write_mesh_inputs` for reuse; seed check via `shape.q`)
- Test: `tests/unit/test_hybridplume_ellipsoid.py` (extend)

**Interfaces:**
- Consumes: Task 4.
- Produces: `spatial.FIT_FILES = {"fitted-ellipsoid": "ellipsoid.json", "fitted-capsule": "capsule.json"}`; `spatial.load_shape(case_dir, cfg) -> tuple[dict, Spheroid|Capsule]`; `certify_facets(vertices, shape, lower, upper, retained, tolerance) -> dict` (same output keys); `fitted._identity(case_dir, cfg, shape, time)`; `ellipsoid_mesh` helpers `background(cfg, source) -> dict` (divisions, size, seed, block dict text), `source_regions(cfg, source, size, base_radius) -> list[RefinementRegion]`, `finish_dicts(system, cfg, source, block, snappy) -> list[Path]`.

- [ ] **Step 1: Regression test first** — golden values: run `fit_ellipsoid` on the `voxels()` fixture, `certify_facets` on its sampled surface, and `fitted._identity` on a tmp case; assert results equal those produced by calling the pre-refactor code paths (capture them in the test via the unchanged `verify_ellipsoid` and a hard-coded identity computed at the start of this task on the current commit). Also assert the rendered `snappyHexMeshDict` for the ellipsoid profile is byte-identical before/after (capture text in the test at the start of the task).
- [ ] **Step 2: Run — PASS on the current code** (establishes the golden).
- [ ] **Step 3: Refactor** to consume shapes; no behaviour change.
- [ ] **Step 4: Run** `tests/unit/test_hybridplume_ellipsoid.py tests/unit/test_hybridplume_workflow.py` — PASS.
- [ ] **Step 5: Commit** — `refactor: route the 3D fitted path through the shape interface`.

### Task 6: Capsule geometric fit (exact sweep)

**Files:**
- Create: `plumetools/hybridplume/capsule.py`
- Test: `tests/unit/test_hybridplume_capsule.py`

**Interfaces:**
- Consumes: `shapes.Capsule`.
- Produces: `feasible_candidates(lower, upper, retained, *, source_radius: float, clearance: float, radius_samples: int, length_samples: int) -> list[dict]`, each `{"R","L","a","volume_m3","verification"}`, all certified by `Capsule.verify`; empty list if none.

Algorithm (spec §2b): domain limits as in `fit_ellipsoid` (`domain_pad = max(clearance, 0.25·min cell)`); excluded slabs = `np.unique` over excluded boxes' x-intervals with min radial distance per slab (same reduction as `fit_ellipsoid`, padded by clearance). `R_max = min(source_radius, r_domain) − clearance`; R grid = `R_max·(1 − k/radius_samples·0.5)`, k = 0..radius_samples−1 (upper half of the range, where the probe showed the useful shapes). `L_max(R)` = min `x_lo` of slabs with `ρ < R`, minus nothing (slab already padded), clipped to `[0, x_max]`. L grid = `L_max·j/length_samples`, j = 0..length_samples. `a_max(R,L)` = min over slabs with `x_hi > L` and `ρ < R` of `max(x_lo − L, 0)/√(1 − ρ²/R²)`, capped by `x_max − L`; skip candidates with `a_max ≤ 0` unless `L > 0` (a flat-ended cylinder is not a capsule — require `a > 0`).

- [ ] **Step 1: Failing tests**

```python
def test_unbounded_column_reaches_source_and_domain_limits():  # all retained -> max-volume candidate R ≈ R_max, L + a ≈ x_max - pad
def test_step_sets_length_limit():          # exclude r < 0.3 for x in [0.6, 0.8] -> every candidate with R > 0.3+clearance has L <= 0.6
def test_tip_only_exclusion_limits_a():     # single excluded box on axis at x in [1.2,1.4] -> L + a <= 1.2 for all candidates
def test_fit_handles_zero_length_cylinder(): # exclusion at x < 0.05, r > 0.2 -> candidates with L == 0 present, a > 0, all verified
def test_every_candidate_is_certified():
```

Build voxel fixtures like `voxels()` in `test_hybridplume_ellipsoid.py`, sized so these numbers are on cell edges.
- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** — `feat: exact conservative capsule candidate sweep`.

### Task 7: Capsule flux screen and selection

**Files:**
- Modify: `plumetools/hybridplume/capsule.py`
- Modify: `plumetools/hybridplume/fitted.py` (`sample_containing_cells`: optional cached grid)
- Test: `tests/unit/test_hybridplume_capsule.py` (extend)

**Interfaces:**
- Consumes: Task 6 candidates; `halfrange.incoming_flux`, `halfrange.integrate`.
- Produces: `fitted.RectilinearLocator(mesh, axes=(0,1,2))` with `.locate(points, allowed) -> cells` (the body of `sample_containing_cells`, grid built once; `sample_containing_cells` becomes a thin wrapper — existing behaviour unchanged). `capsule.surface_samples(shape, axial_points_per_D, azimuthal_points, D) -> (points, normals, areas)` (cylinder wall + tip, no base disc). `capsule.screen(candidates, locator, state, retained, gas, D, cfg) -> list[dict]` adding `"return_mass_fraction", "subsonic_area_fraction", "screen_feasible"`. `capsule.select(screened, crit) -> dict` → `{"status": "fitted", "shape": Capsule.to_json(), "screen": {"status": "passed"|"no_passing_candidate", "candidates": [...]}}` or `{"status": "no_feasible_capsule", "reason": ...}`.

- [ ] **Step 1: Failing tests**

```python
def test_surface_area_matches_analytic():   # sum(areas) within 0.5 % of 2πRL + spheroid-cap area
def test_screen_return_is_small_for_radial_outflow():  # synthetic box mesh, U = 500 m/s radial-from-source, T 300 -> return < 1e-3
def test_screen_return_is_large_for_axial_flow_on_cylinder(): # U axial only -> cylinder wall Mn = 0 -> return > 0.3
def test_select_prefers_largest_passing_volume():
def test_select_reports_no_passing_candidate():  # all above limit -> lowest-return candidate, status no_passing_candidate
def test_screen_marks_unsampleable_candidate_infeasible():  # sample on excluded cell -> screen_feasible False, sweep continues
def test_sample_containing_cells_unchanged():  # existing ellipsoid sampler test still passes via wrapper
```

- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement.** Catch the locator's `WorkflowError` per candidate → `screen_feasible=False`.
- [ ] **Step 4: Run** capsule + ellipsoid tests — PASS.
- [ ] **Step 5: Commit** — `feat: pre-mesh half-range flux screen for capsule candidates`.

### Task 8: Capsule STL and snappy inputs

**Files:**
- Create: `plumetools/hybridplume/capsule_mesh.py`
- Test: `tests/unit/test_hybridplume_capsule.py` (extend)

**Interfaces:**
- Consumes: `ellipsoid_mesh.background`, `source_regions`, `finish_dicts` (Task 5); `shapes.Capsule`.
- Produces: `capsule_stl(shape, chord_tol: float, back: float) -> (vertices (n,3), triangles (m,3))`; `write_stl(path, vertices, triangles)` (ASCII, solid `interface`); `write_mesh_inputs(case_dir, cfg, source, shape) -> list[Path]` writing `mesh-build/constant/triSurface/interface.stl`, `system/surfaceFeatureExtractDict` (`interface.stl { extractionMethod extractFromSurface; includedAngle 150; writeObj no; }`), and a snappyHexMeshDict whose geometry is `interface.stl { type triSurfaceMesh; name interface; }`, `refinementSurfaces interface { level (L L); patchInfo { type patch; } }` at `dsmc.mesh.surface_level`, and `features ( { file "interface.eMesh"; level L; } )`.

STL: cylinder from x = −back to L, tip as rings of the spheroid; back cap disc at x = −back. Ring count/segments chosen so chord sagitta / local radius ≤ chord_tol (`chord_tol = 0.1·geometry_relative_tolerance`, `back` = one background cell).

- [ ] **Step 1: Failing tests**

```python
def test_stl_is_closed_and_outward():  # every edge shared by exactly 2 triangles; signed volume > 0 and ≈ capsule volume + πR²·back within 1 %
def test_stl_chord_error_within_tolerance():  # max |sqrt(q)-1| over triangle centroids ≤ chord_tol (centroids with x ≥ 0)
def test_snappy_dict_uses_trisurface_and_features(tmp_path):  # strings present; seed outside capsule (shape.q(seed) > 1)
```

- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run — PASS.**
- [ ] **Step 5: Commit** — `feat: watertight capsule STL and snappy inputs with feature edges`.

### Task 9: Wire the capsule into `spatial.interface` and Allrun

**Files:**
- Modify: `plumetools/hybridplume/spatial.py::interface` (dispatch by mode), `plots` (section from shape; output name `breakdown_capsule.png` for capsule)
- Modify: `cases/cai2012-hybrid/Allrun` (in the `spatial_fit` mesh block, run `surfaceFeatureExtract` after `blockMesh` when `system/surfaceFeatureExtractDict` exists)
- Modify: `docs/hybridplume.md` (short section: capsule mode, viscosity cap, new profile; link the spec)
- Test: `tests/unit/test_hybridplume_capsule.py` (extend)

**Interfaces:**
- Consumes: Tasks 4–8.
- Produces: for `fitted-capsule`, `spatial.interface` writes `results/capsule.json` (output of `capsule.select`) and `interface.json` with `"mode": "fitted-capsule"`, `"fit": {...}`, `"screen": {"status", "selected": {...}}`, `status "mesh_pending"`; then `capsule_mesh.write_mesh_inputs`. Uses `cfd.read_state` and the cached locator for the screen.

- [ ] **Step 1: Failing tests**

```python
def test_capsule_interface_writes_fit_screen_and_mesh_inputs(tmp_path, monkeypatch):  # synthetic case via _hybridmesh.box_mesh + breakdown.npz; asserts files and statuses
def test_capsule_interface_refuses_existing_ellipsoid_mesh(tmp_path, monkeypatch):    # dsmc/constant/polyMesh/owner present with ellipsoid identity -> WorkflowError
def test_allrun_runs_surface_feature_extract_for_stl_meshes():  # text check: "surfaceFeatureExtract" inside the spatial_fit mesh block
```

Follow the monkeypatch pattern in `test_retry_refuses_a_changed_dsmc_mesh_before_reusing_cloud` (`test_hybridplume_ellipsoid.py` L185).
- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** the full hybridplume unit set: `.venv/bin/python -m pytest tests/unit -q -k hybridplume` and `bash -n cases/cai2012-hybrid/Allrun` — PASS.
- [ ] **Step 5: Commit** — `feat: fit, screen and mesh the capsule interface in the 3D pipeline`.

**Owner checkpoint after Task 9:** `PROFILE=hpc-kn1em3-capsule msub cases/cai2012-hybrid/run3DCase.sh` (VEGA uses `msub`). Expected: CFD stable per spec §1c; `capsule.json` screen status; actual-mesh `interface.json` return fraction ≈ 0.05 (borderline); resource guard verdict for the 30 D DSMC box.

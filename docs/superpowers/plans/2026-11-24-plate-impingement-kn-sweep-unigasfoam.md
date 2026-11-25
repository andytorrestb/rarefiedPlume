# Plate-impingement Kn sweep on uniGasFoam — Implementation Plan (phase A: to submission)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate the 21 plate-impingement cases (7 Kn_e × z/D 25/50/100) for uniGasFoam, with cluster job scripts, and pass an under-resolved laptop smoke run — everything the owner needs before submitting.

**Architecture:** A new package `plumetools/plateugf/` glues two existing, unedited packages: `plumetools.nozzleplume` (case.yaml schema, plate/symmetry geometry, snappyHexMesh octree, nozzle patch, run window, decomposeParDict) and `plumetools.cai2012ugf` (uniGasFoam dictionary renderers, `unigas.yaml` settings, restart helpers). Mean-free-path sizing is computed by the generator and written into each case's `case.yaml` as nozzleplume `mesh.refinement_regions` cones, so nozzleplume meshes it unchanged. A small adapter exposes the attributes the cai2012ugf renderers read.

**Tech Stack:** Python 3 (numpy, PyYAML, pytest), OpenFOAM v2512 (blockMesh, snappyHexMesh, topoSet, createPatch), uniGasFoam fork `dsmc-time-step-limit` ≥ bd8ba5c, bash, PBS.

**Spec:** `docs/superpowers/specs/2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md`

**Phase B (separate plan, after submission):** conversion of plate fields, plate-load sampling, visualisation, Cai collisionless plate solution, comparison, the three shortened Kn_e = 100 prototype runs. Nothing in phase B changes what the cluster jobs write; Task 6's smoke check guards that.

## Global Constraints

- Never edit `plumetools/nozzleplume`, `plumetools/cai2012ugf`, `plumetools/cai2012`, `plumetools/cai2012ext`, or any existing case family.
- Kn_e ∈ {100, 10, 1, 0.1, 0.01, 0.001, 0.0001}; z/D ∈ {25, 50, 100}; D = 0.2 m; S0 = 2; T0 = 300 K; argon (mass 6.63e-26 kg, d 4.17e-10 m, ω 0.74, Tref 273 K, 0 internal DOF).
- Plate 34.14 (y) × 11.88 (z) × 0.25 (x) m, centred, upstream face at x = z/D·D; diffuse, 300 K.
- Quarter domain: `symmetry_planes: [y, z]`; margins 1.75 m behind the plate, 2 m beside it.
- Sizing: s = clamp(min(λ0·n0/n, r/20), 0.01 m, 0.64 m); background 0.64 m; levels 0–6 (0.01 m at level 6); plate level 3.
- Weight: `target_particles_per_exit_cell: 40` (weight null → derived).
- Run window: `transient_domain_transits: 6`, `sampling_domain_transits: 3`; startup writes every transit, sampling writes 4 per transit.
- uniGasFoam settings: `cases/cai2012-unigas/baseCase/unigas.yaml`, copied verbatim; `meta.unigas_commit: bd8ba5c3d9ffbbf95a30629eae7012b913188f4a`.
- No cell budget: `mesh.max_cells: 1000000000` (typo guard only). Nothing pruned on the cluster.
- PBS headers as the first family's (`longq`, `nodes=1:ppn=192`, 120 h, the owner's `-M`).
- Case names `Kn<label>_zD<25|50|100>`, Kn labels as the first family (`Kn100`, `Kn0p1`, `Kn0p0001` …).

## Review Focus

1. A Kn_e ≤ 0.001 case whose λ-region extends past the plate or the domain: cones must be clipped to the domain, not refused.
2. The exit disk must sit wholly in the finest region (nozzleplume refuses otherwise): every case, including the smoke override.
3. A resumed run must pick the right leg from the newest complete write (restart helpers reused, leg from plate run window).
4. Regenerating over a case that already has solver output must refuse, not overwrite.
5. The symmetry patches must get no uniGas boundary entry, and every other non-nozzle, non-plate patch must get a deletion entry.

---

### Task 1: Mean-free-path sizing → refinement cones

**Files:**
- Create: `plumetools/plateugf/__init__.py`, `plumetools/plateugf/sizing.py`
- Test: `tests/unit/test_plateugf_sizing.py`

**Interfaces:**
- Produces:
  - `density_table(x_max_m: float, rho_max_m: float, radius_m: float, speed_ratio: float) -> DensityTable` — collisionless n/n0 on a log-spaced (x, ρ) grid via `plumetools.cai2012.analytical.density_ratio` (n_radial=32, n_azimuthal=64); `DensityTable.n_ratio(x, rho) -> ndarray` (log-interpolated, 0 for x ≤ 0).
  - `target_cell_size(table, x, rho, *, knudsen: float, diameter_m: float, floor_m=0.01, ceiling_m=0.64, geo_fraction=0.05) -> ndarray`.
  - `refinement_regions(table, *, knudsen, diameter_m, background_m=0.64, max_level=6, floor_m=0.01, geo_fraction=0.05, x_limit_m, rho_limit_m) -> list[dict]` — one dict per level 1..max_level, keys exactly nozzleplume `RegionConfig` fields: `{name: "mfpL<level>", level, shape: "cone", x_min_m: 0.0, x_max_m, radius1_m, radius2_m}`, rounded up to 1 mm.

Region for level L: points where s < background/2^(L−1). Cone fit: `x_max` = largest sampled x in the region plus one level-L cell, clipped to `x_limit_m`; `radius1` = max(envelope at the smallest x, D/2) plus one level-L cell; slope = max over sampled x of (envelope(x) − radius1)/x (≥ 0); `radius2 = radius1 + slope·x_max`, clipped to `rho_limit_m`. Envelope(x) = largest sampled ρ in the region at that x.

- [ ] **Step 1: Write the failing tests**

```python
def test_rarefied_regions_are_geometric():           # Kn 100 and Kn 0.1 give identical cones
def test_every_region_point_is_inside_its_cone():     # sample 20k (x, rho) in the region; all inside
def test_regions_nest_and_finest_holds_exit_disk():   # level 6 cone radius1 >= 0.1 + 0.01; x_max >= 0.01
def test_dense_case_reaches_further_than_rarefied():  # Kn 1e-4 level-1 x_max > Kn 100 level-1 x_max
def test_cones_clip_to_domain():                      # x_limit 7.0 → every x_max <= 7.0
def test_floor_and_ceiling():                         # target_cell_size within [0.01, 0.64] everywhere
```

- [ ] **Step 2: Run** `python -m pytest tests/unit/test_plateugf_sizing.py -q` — expect ImportError.
- [ ] **Step 3: Implement** `sizing.py` as above.
- [ ] **Step 4: Run** the tests — expect all pass.
- [ ] **Step 5: Commit** `feat(plateugf): mean-free-path sizing as refinement cones`.

### Task 2: uniGasFoam dictionaries through an adapter, with the plate

**Files:**
- Create: `plumetools/plateugf/dictionaries.py`
- Test: `tests/unit/test_plateugf_dictionaries.py`

**Interfaces:**
- Consumes: `plumetools.nozzleplume.workflow.Case` (`cfg`, `geom`, `exit_state`, `run`, `plan`); `plumetools.cai2012ugf.settings.UniGasSettings`; cai2012ugf renderers `render_uni_gas_properties`, `render_field_properties_dict`, `render_initialisation_dict`, `render_hybrid_decomposition_dict`, `render_control_dict`, `render_fv_schemes`, `render_fv_solution`, `render_chem_react_dict`; nozzleplume `render_decompose_par_dict`, `mesh.write_mesh_setup`.
- Produces:
  - `ugf_views(case) -> tuple[SimpleNamespace, SimpleNamespace, SimpleNamespace]` — (cfg, exit_state, run) carrying exactly: `cfg.gas.{species_name, t_ref_K, mass_kg, diameter_m, internal_degrees_of_freedom, omega}`, `cfg.geometry.plume_axis = (1, 0, 0)`, `cfg.mesh.patch_names`, `cfg.dsmc.{initial_number_density_per_m3, initial_temperature_K}`; `exit_state.{T0_K, number_density_per_m3, velocity_m_per_s}`; `run.{n_equivalent_particles, average_start_s, end_time_s, delta_t_s, domain_transit_s (= run.transit_s), write_interval_s (= transit_s / 4)}`.
  - `render_boundaries_dict(case, settings) -> str` — nozzle: `uniGasDeletionPatch` + `uniGasFreeStreamInflowPatch` (n0, T0, U0 along +x); plate: `uniGasDiffuseWallPatch` {velocity (0 0 0); temperature <plate T>;}; upstreamVacuum, vacuum: `uniGasDeletionPatch`; symmetry: no entry.
  - `write_case_dictionaries(case, settings) -> list[Path]` — the same file set as `cai2012ugf.dictionaries.write_case_dictionaries`, plus nozzleplume's mesh setup.

- [ ] **Step 1: Write the failing tests** (build a Case from a tmp `case.yaml` copied from Task 5's template, Kn 10, z/D 25):

```python
def test_exit_state_matches_first_family():  # n0, U0 equal cai2012.inflow.from_config(cases/cai2012-unigas/baseCase) at Kn 10 to 1e-12 rel
def test_plate_is_diffuse_wall_at_300K():
def test_symmetry_has_no_boundary_entry_and_open_patches_delete():
def test_other_dictionaries_equal_cai2012ugf_for_same_views():  # render_uni_gas_properties(views) == file written
def test_sampling_writes_four_per_transit():  # controlDict.sampling writeInterval ≈ transit/4
```

- [ ] **Step 2: Run** — expect ImportError.
- [ ] **Step 3: Implement** `dictionaries.py`.
- [ ] **Step 4: Run** — expect pass.
- [ ] **Step 5: Commit** `feat(plateugf): uniGasFoam dictionaries for nozzleplume cases, with the plate`.

### Task 3: Case workflow CLI and restart helpers

**Files:**
- Create: `plumetools/plateugf/workflow.py`, `plumetools/plateugf/__main__.py`
- Test: `tests/unit/test_plateugf_workflow.py`

**Interfaces:**
- Consumes: Task 2; `nozzleplume.workflow.load/report`; `cai2012ugf.hpc.{leg_for, lost_clouds, mark_restart, UgfRestartPlan, _slack}`; `foamio.checkpoints.plan_restart`.
- Produces: `python -m plumetools.plateugf <step> <case_dir>`:
  - `describe` — nozzleplume report + memory (4.3 KB/cell + 240 B/particle) and disk estimate; writes nothing.
  - `setup` — describe, checks (`raise_if_failed`), `write_case_dictionaries`.
  - `resume-plan --np N` — `UgfRestartPlan.shell()` with the leg from the plate run window; exit 3 on refusal.
  - `mark-restart --time T --np N`, `finished`, `check-build` (expects `meta.unigas_commit`, compares `$FOAM_USER_APPBIN/uniGasFoam.commit`).
  - `smoke-check [--time T]` — Task 6's pass criteria on a reconstructed or decomposed write; exit 0/1 with one line per criterion.

- [ ] **Step 1: Failing tests:** `test_resume_plan_fresh_case_is_startup`, `test_leg_after_average_start_is_sampling`, `test_check_build_mismatch_reports`, `test_smoke_check_fails_without_plate_fields`.
- [ ] **Step 2–4:** run (fail) → implement → run (pass).
- [ ] **Step 5: Commit** `feat(plateugf): case workflow, restart plan, build and smoke checks`.

### Task 4: Case scripts

**Files:**
- Create: `cases/plate-unigas/baseCase/{Allmesh, Allrun, Allclean, runMesh.sh, runCase.sh}`

`Allmesh` = r-4d's `Allmesh` with `python -m plumetools.plateugf setup .` in place of the nozzleplume setup, and its results guard; meshing serial (`blockMesh`, `snappyHexMesh -overwrite`, `topoSet`, `createPatch -overwrite`, `checkMesh`, `python -m plumetools.nozzleplume verify .`). `Allrun` = the first family's `Allrun` with `plumetools.cai2012ugf.hpc` → `plumetools.plateugf`, the `simple_split` override from `nozzleplume.dictionaries.simple_split`, `--stop-at-transits` computed from the plate run, and the final convert step removed (phase B). `runMesh.sh`: PBS, `longq`, 1 node, 48 h (serial snappyHexMesh is r-4d's proven path, at up to 30M cells). `runCase.sh`: the first family's, `-np 192`, resubmit-until-done message. `Allclean` as r-4d's guarded one.

- [ ] **Step 1:** write the scripts; `bash -n` each — expect no output.
- [ ] **Step 2: Commit** `feat(plate-unigas): case scripts`.

### Task 5: Study, template and generator

**Files:**
- Create: `cases/plate-unigas/{study.yaml, generate_cases.py, README.md}`, `cases/plate-unigas/baseCase/{case.yaml, unigas.yaml}`, `plumetools/plateugf/generation.py`
- Test: `tests/unit/test_plateugf_generation.py`

**Interfaces:**
- Consumes: Task 1 (`refinement_regions`), Task 3 (`describe` numbers).
- Produces: `generate(family_root: Path, *, dry_run=False, only: list[str] | None = None) -> list[CaseEntry]`; `stage_smoke(family_root: Path) -> Path` (writes `Smoke/Kn100_zD25`: background 1.28 m, floor 0.02 m, geo_fraction 0.2, max_level 6, transits 1 + 1, `n_subdomains: 12`, `decomposition_method: scotch`); `CaseEntry(name, knudsen, standoff_zD, cells_estimate, particles_estimate, memory_gb, steps)`. `generate_cases.py [--dry-run] [--smoke]` prints the table. Each case directory: `case.yaml` (template + `exit.knudsen`, `plate.standoff_m`, `domain.x_max_m`, `mesh.refinement_regions`, `meta.case_name`), `unigas.yaml`, the scripts. Refuses to overwrite a case with solver output; rewrites otherwise.

- [ ] **Step 1: Failing tests:** `test_21_cases_named_and_parameterised`, `test_dry_run_writes_nothing`, `test_refuses_case_with_output`, `test_every_case_loads_with_nozzleplume_and_exit_in_finest_region`, `test_smoke_case_is_small` (estimated cells < 200k).
- [ ] **Step 2–4:** fail → implement → pass.
- [ ] **Step 5:** `python cases/plate-unigas/generate_cases.py --dry-run` — table of 21 rows, no file written.
- [ ] **Step 6:** `python cases/plate-unigas/generate_cases.py` — 21 case directories; then `python -m plumetools.plateugf describe cases/plate-unigas/Cases/Kn0p0001_zD100` — no failed checks.
- [ ] **Step 7: Commit** `feat(plate-unigas): study, template and generator for the 21 cases` (generated `Cases/` included, as the first family does).

### Task 6: Laptop smoke run

Not code: the run the owner asked to find started.

- [ ] **Step 1:** In WSL (OpenFOAM v2512 + fork sourced; `df -h /mnt/c` > 20 GB free), `python cases/plate-unigas/generate_cases.py --smoke`, then copy `Smoke/Kn100_zD25` to `~/ugf-laptop/plate/Kn100_zD25`.
- [ ] **Step 2:** `./Allmesh` — expect `Mesh OK`, patches nozzle/upstreamVacuum/vacuum/plate/symmetry*, verify passes.
- [ ] **Step 3:** `./Allrun -np 12 --skip-build-check` detached (`nohup … > log.Allrun 2>&1 &`) — expect it past step 20 within minutes.
- [ ] **Step 4:** when the startup leg's first write exists: `python -m plumetools.plateugf smoke-check .` — pass criteria (spec §8): plate patch carries `p_Ar`, `surfaceShearStress_Ar`, `fD_Ar`, `surfaceHeatTransfer_Ar`, finite, `p_Ar > 0` on faces under the plume once the plume has reached the plate; mean volume fields present; no particle count lost through symmetry planes (gas present in cells adjacent to y = 0 and z = 0).
- [ ] **Step 5:** restart check: stop the run, `./Allrun -np 12 --skip-build-check` again — expect "leg … from t = <newest write>".
- [ ] **Step 6:** report to the owner: what passed, what didn't, and the commands to submit.

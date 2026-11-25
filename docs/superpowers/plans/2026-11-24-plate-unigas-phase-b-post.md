# Plate-impingement Kn sweep — Implementation Plan (phase B: post-processing and the Cai comparison)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Plate-load sampling, flow visualisation and the collisionless Cai plate comparison for any plate-unigas case, proven end to end on three under-resolved Kn_e = 100 laptop runs (z/D = 25, 50, 100).

**Architecture:** `plumetools/plateugf/` gains `caiplate.py` (collisionless plate loads by quadrature over the exit disk), `post.py` (uniGasFoam plate fields into nozzleplume's `PlateLoads`, maps, profiles, the comparison) and `viz.py` (matplotlib slices on the symmetry planes from `writeCellCentres`/`writeCellVolumes`). A case `Allpost` runs them. The cluster jobs are unchanged.

**Tech Stack:** Python 3, numpy, matplotlib, OpenFOAM v2512 `postProcess`/`reconstructPar`.

**Spec:** `docs/superpowers/specs/2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md` §9 (owner, 2026-11-24: "yes. just test across Kn_e = 100").

## Global Constraints

- Never edit `plumetools/nozzleplume`, `plumetools/cai2012ugf`, `plumetools/cai2012`, `plumetools/cai2012ext`; import only.
- uniGasFoam plate fields: `p_Ar`, `surfaceShearStress_Ar`, `fD_Ar` (force per area, momentum the gas hands the wall), `surfaceHeatTransfer_Ar`; mean field name from `unigas.yaml` `fields.mean_field`.
- Diffuse plate at `plate.temperature_K` (300 K): p = Π_nn + m Γ √(π k T_w / 2m), τ = Π_nt, q = E − 2 k T_w Γ.
- Differences are reported, never judged against a tolerance; no claim of validation.
- Prototype runs are laptop-only, coarse, 2 startup + 1 averaging transits; never results.

## Review Focus

1. Lee and edge faces of the plate must not enter the upstream-face maps or the Cai comparison (Cai's solution is for the upstream face).
2. Totals must be whole-plate (quarter × 4), consistently for simulation and Cai.
3. A time directory without the averaged plate fields (run stopped before averaging) must fail with a clear message, not plot zeros.
4. Face centres on the symmetry planes' edges (y = 0 or z = 0 rows) must be inside the centreline profiles.
5. Colour ranges fixed across cases, so Kn_e = 100's three stand-offs are comparable.

---

### Task B1: Prototype staging and launch

**Files:** Modify `plumetools/plateugf/generation.py`, `cases/plate-unigas/study.yaml`; Test `tests/unit/test_plateugf_generation.py`.

- Produces: `stage_prototype(family_root: Path, case: str) -> Path` → `Prototype/<case>`: the smoke overrides plus `study.yaml` `prototype:` (`transient_domain_transits: 2`, `sampling_domain_transits: 1`, `target_particles_per_exit_cell: 10`, `n_subdomains: 4`).
- [ ] Failing test `test_prototype_cases_are_short_and_coarse` (2 + 1 transits, weight target 10, 4 ranks, < 250k cells, exit level = max) → implement → pass → commit.
- [ ] Stage Kn100_zD25/50/100; laptop runner runs all three concurrently, 4 ranks each, `simple (4 1 1)` (slabs along the plume axis), then `Allpost`.

### Task B2: Collisionless plate loads (`caiplate.py`)

- Produces: `plate_loads(points: ndarray (n,3), *, radius_m, number_density, speed, temperature_K, mass_kg, wall_temperature_K, n_radial=48, n_azimuthal=96) -> CaiPlate` with `.n, .flux, .pressure, .shear (n,2: y,z), .heat_flux` per point; also `moments(points, ...)` for tests.
- [ ] Failing tests: `test_density_matches_cai_free_plume` (n at 6 off-axis points vs `cai2012.analytical.density_ratio`, rel 2e-3); `test_axial_speed_matches_cai_centreline` (Γ/n on the axis vs `centerline_axial_speed_ratio`); `test_far_field_pressure_falls_as_inverse_square` (p(2z)/p(z) → 1/4 within 3% at z = 20, 40 m); `test_shear_is_zero_on_axis_and_points_outward`.
- [ ] Implement (radial moments M_k(μ) tabulated on μ ∈ [0, 1], interpolated per disk node) → pass → commit.

### Task B3: Plate post-processing and the comparison (`post.py`)

- Produces: `read_plate(case_dir, time=None) -> PlateLoads` (nozzleplume's dataclass, from `fD_<field>`, `surfaceHeatTransfer_<field>`); `upstream(loads) -> mask`; `profiles(loads, values) -> dict` (rows nearest y = 0 and z = 0); `compare(case) -> dict`; `run_post(case_dir, time=None) -> list[Path]` writing `results/plate_loads.csv`, `plate_summary.yaml`, `plate_profiles.csv`, `cai_comparison.csv`, `cai_comparison.yaml`, figures `plate-pressure.png`, `plate-shear.png`, `plate-heat-flux.png`, `plate-profiles.png`, `cai-pressure-shear-profiles.png`, `cai-ratio-maps.png`.
- CLI: `python -m plumetools.plateugf post <case> [--time T]`; `python -m plumetools.plateugf family-table <family>`.
- [ ] Failing tests on a synthetic case (tiny polyMesh with a plate patch of 6 front faces + 2 edge faces, synthetic fields): `test_read_plate_uses_unigas_fields`, `test_edge_faces_are_left_out_of_the_comparison`, `test_missing_mean_fields_are_refused`, `test_profiles_take_the_symmetry_rows` → implement → pass → commit.

### Task B4: Flow slices (`viz.py`)

- Produces: `slices(case_dir, time, out_dir) -> list[Path]` — n/n0 (log, 1e-6..1), |U| (0..1200 m/s), T (0..900 K) on y = 0 (x–z) and z = 0 (x–y) from cells touching each plane, plate outline drawn; fixed ranges.
- [ ] Failing test `test_slice_cells_touch_the_plane` (synthetic C/V) → implement → pass → commit.

### Task B5: `Allpost`, `runPost.sh`, regenerate

- [ ] `baseCase/Allpost`: `reconstructPar -newTimes` if decomposed, `postProcess -func writeCellCentres -latestTime`, `writeCellVolumes`, `python -m plumetools.plateugf post .`, `python -m plumetools.plateugf viz .`. `runPost.sh` (PBS, shortq, 1 node). `bash -n`; regenerate Cases/ (no output yet) → commit.

### Task B6: The three Kn_e = 100 prototypes, end to end

- [ ] When each prototype finishes: `./Allpost`; check figures exist and the Cai comparison runs; `family-table` over the three; report numbers (centre p sim vs Cai, total force ratio) without judging them.

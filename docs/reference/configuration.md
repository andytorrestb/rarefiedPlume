# Configuration files

Each family has its own loader and its own schema; there is no universal
`case.yaml`. The loaders are authoritative: defaults below are the dataclass
defaults in the named module, and every loader **rejects unknown keys** (a typo
is an error, not a silently ignored setting). A key shown as *required* has no
default. Units are in the key names. Values a committed case sets are given
where they matter; the committed files themselves annotate every value
(`[PAPER]`, `[ASSUMPTION]`, `[DERIVED]`, `LEGACY`, …).

| File | Family | Loader | Role |
|---|---|---|---|
| `case.yaml` (`model: source_flow`) | 3d-inflow | `plumetools.config.load_case_config` | [legacy source flow](#legacy-source-flow-caseyaml-3d-inflow) |
| `case.yaml` (`model: markelov1999_axisymmetric`) | markelov1999 | `plumetools.config.load_case_config` | [same loader, extra sections](#markelov1999-sections) |
| `study.yaml` | markelov1999 | `plumetools.markelov1999.study` | [pressure matrix](#markelov1999-studyyaml) |
| `case.yaml` (`model: cai2012_circular_plume`) | cai2012, health, Kn extension | `plumetools.cai2012.config.load_case_config` | [Cai case](#cai2012) |
| `study.yaml` | cai2012, Kn extension | `plumetools.cai2012.study.load_study` | [Kn × particle matrix](#cai2012-studyyaml) |
| `study.yaml` | cai2012-health | `plumetools.cai2012.health.load_study` | [weight × sampling matrix](#cai2012-health) |
| `interfaces.yaml`, `<case>/interface.yaml` | Kn extension | `plumetools.cai2012ext.interface_study`, `interface_config` | [interface track](#interfaceyaml-and-interfacesyaml) |
| `scaling.yaml`, `scaling-registration.yaml` | Kn extension | `plumetools.cai2012ext.scaling_workflow`, `scaling_inventory` | [field scaling](#scalingyaml) |
| `diagnostics.yaml`, `diagnostics:` section | Kn extension, nozzle plume | `plumetools.rarefaction.settings` | [rarefaction diagnostics](#rarefaction-settings) |
| `assessment.yaml` | virtual collisions | `plumetools.virtualcollision.assessment.load_yaml` | [assessment](#assessmentyaml) |
| `hybrid.yaml` | cai2012-hybrid | `plumetools.hybridplume.config` | [hybrid](#hybridyaml) |
| `case.yaml` (`model: nozzle_exit_plume`) | r-4d | `plumetools.nozzleplume.config` | [nozzle plume](#nozzleplume-caseyaml-r-4d) |
| `viz.yaml`, `plumetools/viz/slices.yaml` | every imaging family | `plumetools.viz.spec` | [slice imagery](#slice-imagery-specs) |

The `model:` key is how loaders refuse each other's files: a Cai loader given a
source-flow `case.yaml` stops with a message naming the right loader.

## Legacy source-flow `case.yaml` (3d-inflow)

Top-level keys: `model` (`source_flow`), `execution_profile` (`smoke` or
`production`, default `production`), `meta` (free-form), and the sections
below. `load_case_config` also cross-checks the file against the case: the
inflow patch must exist in `constant/polyMesh/boundary`, `gas.species_name`
must be in `constant/dsmcProperties` `typeIdList`, and `output.patches` names
must exist unless `legacy.emit_absent_patches` is true (then a warning).

| Key | Default | Allowed / constraint |
|---|---|---|
| `mesh.type` | `block_mesh_ogrid` | `block_mesh_ogrid`, `snappy_hex_sphere`, `snappy_markelov` |
| `mesh.sphere_radius_m` | 0.5 | must equal `geometry.sphere_radius_m` |
| `mesh.box_half_width_m`, `box_length_m` | 2.5, 5.0 | m |
| `mesh.outer_patch_type` | `patch` | `wall` makes the outer boundary reflect (needed by stock `FreeStream`) |
| `mesh.patch_names` | `{inflow, outer, symmetry}` roles | |
| `mesh.background_cell_size_m`, `refinement_level`, `n_cells_between_levels` | 0.125, 3, 3 | snappy types; surface cell = background / 2^level |
| `mesh.n_tangential`, `n_radial`, `radial_grading`, `projection` | 20, 24, 10.0, `arc` | O-grid only; `projection`: `arc`, `searchable_sphere`, `none` (`arc` warns: the surface is not a true hemisphere) |
| `geometry.patch` | `inflow` | |
| `geometry.sphere_radius_m` | 0.5 | |
| `geometry.plume_axis`, `polar_axis` | `[1,0,0]`, `[0,0,1]` | |
| `gas.gamma`, `molar_mass_g_per_mol`, `species_name` | 1.4, 28.0134, `Ar` | the model gas; `species_name` is the solver species |
| `stagnation.p0_pa`, `T0_K`, `throat_radius_m` | —, —, 0.0041275 | `throat_radius_m` < `sphere_radius_m` |
| `angular.exponent_offset`, `form`, `normalization_gamma`, `quadrature_points` | 0.41, `legacy_separable`, 1.4, 500 | |
| `output.dialect` | `standard` | `standard` (OpenFOAM `dsmcFoam`), `mnf` (`dsmcFoam+`) |
| `output.inflow_model` | `FreeStream` | `FreeStream`, `plumeFieldInflow`, `dsmcFreeStreamInflowFieldPatch` |
| `output.patches` | `{}` | per patch: `type` and optional `value` |
| `sampling.*` | radial lines 0/45/90°, 0.5–5 m, 1000 points, `U_Ar`/`Ttra_Ar`/`rhoN_Ar` | validated, **not consumed** (RP-17) |
| `legacy.rhoN_T0_K` | 300.0 | `null` = use `stagnation.T0_K` |
| `legacy.clip_beyond_limiting_angle`, `emit_absent_patches` | false, true | frozen-behaviour switches |

### markelov1999 sections

The same loader, with `mesh.type: snappy_markelov` and `model:
markelov1999_axisymmetric`, reads further sections (defaults from
`plumetools.config`):

| Section | Keys (default) |
|---|---|
| `mesh` (extra) | `x_min_m` … `z_max_m` (all six required for `snappy_markelov`, each pair ordered), `inflow_refinement_level`, `cylinder_refinement_level`, `plate_refinement_level`; `patch_names` must fill `inflow`, `cylinder`, `plate`, `outer`, `symmetry`, `upstream_vacuum` |
| `bodies` | `gap_m` (0.1524), `symmetry_plane_y_m` (0.0); `cylinder.{radius_m 0.0762, length_m 0.4572, centre_x_m 0.29845, centre_z_m 0.0}`; `plate.{width_y_m 0.1524, height_z_m 0.381, thickness_m 0.0127}` |
| `dsmc` | `species.{name N2, molar_mass_g_per_mol 28.0134, mass_kg, gamma 1.4, diameter_m 4.17e-10, omega 0.74, internal_degrees_of_freedom 2, source}` (the implied γ must match `species.gamma`); `n_equivalent_particles` (null = derived), `delta_t_s` 2e-7, `end_time_s` 4e-3, `write_interval_s` 5e-4, `average_start_s` 2e-3, `binary_collision_model` `LarsenBorgnakkeVariableHardSphere`, `rotational_collision_number` 5.0, `wall_interaction_model` `MaxwellianThermal`, `n_subdomains` 4 |
| `resolution` | `target_particles_per_cell` 20, `sizing_region` `cylinder` (`cylinder`, `plate`, `wake`), `wake_offset_m` 0.0254, `report_only` false |
| `checks` | `max_cell_over_mfp` 1.0, `max_courant` 1.0, `warn_courant` 0.5, `min_particles_per_cell` 5.0, `enabled` true |
| `post` | `surface_pressure_field` `fDMean`, `cylinder_patch`, `plate_patch`, `windows` (name, `azimuth_deg`, `half_angle_deg`, `axial_half_height_m`; default windward at 180°) |

### markelov1999 `study.yaml`

`base_case`, `cases_dir`, `manifest_name`, `gap_in` (6.0), `gap_dir`
(`gap06in`), `pressure_cases` (list of `name`, `pressure_psi`, `enabled`,
`note`), `meta`. The generator writes the pressure into each clone's
`case.yaml` structurally.

## cai2012

`plumetools.cai2012.config`; `model: cai2012_circular_plume`; top-level
sections below plus `meta`. `exit.knudsen` is the only per-case physical input;
**null** in a `null`-able key means *derive it*.

| Key | Default | Constraint / note |
|---|---|---|
| `nozzle.diameter_m` | 0.2 | R₀ = D/2 is derived; there is no radius key |
| `exit.speed_ratio` | 2.0 | S₀ = U₀/√(2RT₀) |
| `exit.T0_K` | 300.0 | an assumption: normalised outputs do not depend on it |
| `exit.knudsen` | 100.0 | written per case by the generator; may not be overridden |
| `exit.characteristic_length` | `diameter` | `diameter` or `radius` (radius halves λ₀ and doubles n₀) |
| `gas.species_name`, `mass_kg`, `diameter_m`, `omega`, `t_ref_K` | `Ar`, 6.63e-26, 4.17e-10, 0.74, 273.0 | `diameter_m` and `t_ref_K` are a pair |
| `gas.internal_degrees_of_freedom` | 0 | argon: must be 0 |
| `gas.mean_free_path_convention` | `vhs` | `vhs` (carries (T/T_ref)^(ω−½)) or `hard_sphere` |
| `geometry.symmetry_mode` | `none` | `none` or `half_y` (opt-in only) |
| `geometry.x_max_over_D`, `y_half_over_D`, `z_half_over_D`, `plume_axis` | 10, 10, 10, `[1,0,0]` | |
| `mesh.type` | `graded_cartesian` | the only generator |
| `mesh.core_cell_size_m` | null | null → min(`target_cell_over_mfp`·λ₀, `max_core_cell_over_D`·D) |
| `mesh.target_cell_over_mfp`, `max_core_cell_over_D` | 1.0, 0.05 | `max_core_cell_over_D` ≤ 0.25 |
| `mesh.core_x_over_D`, `core_half_over_D` | 2.0, 1.5 | `core_half_over_D` > 0.5 (the exit lies in the core); both inside the domain |
| `mesh.outer_x_cells`, `outer_lateral_cells`, `outer_expansion` | 40, 30, 40.0 | ≥ 1, ≥ 1, ≥ 1.0 (a cap; the actual ratio is solved for) |
| `mesh.max_cells` | 4 000 000 | hard budget; the core is coarsened to fit, and reports say so |
| `mesh.min_core_cell_size_m` | null | |
| `mesh.patch_names` | `nozzle`, `upstream_vacuum`, `outer`, `symmetry` roles | the first three required and distinct |
| `dsmc.binary_collision_model` | `VariableHardSphere` | the only value accepted |
| `dsmc.collisions_enabled` | true | on at every Kn |
| `dsmc.n_equivalent_particles` | null | null → derived for `resolution.target_particles_per_cell` in the exit cell |
| `dsmc.numerical_particle_multiplier` | 1.0 | weight = baseline / multiplier; set by the generator only |
| `dsmc.run_time_multiplier`, `output_frequency_multiplier` | 1.0, 1.0 | set study-wide by `study.yaml` only |
| `dsmc.delta_t_s`, `courant_target` | null, 0.2 | Δt from the smallest cell and U₀ + 3σ |
| `dsmc.end_time_s`, `average_start_s`, `write_interval_s` | null | derived from the transient basis |
| `dsmc.transient_basis` | `cai` | `cai` (`transient_collision_times`, 10 000) or `transits` (`transient_domain_transits`, 3.0) |
| `dsmc.sampling_domain_transits` | 3.0 | averaging length |
| `dsmc.initial_number_density_per_m3`, `initial_temperature_K` | 1e10, 300 | near-vacuum seed |
| `dsmc.n_subdomains` | 4 | ≥ 1 |
| `resolution.target_particles_per_cell`, `reference_knudsen`, `report_only` | 20.0, 0.01, false | Cai's Δx/λ₀ and Δt/t₀ are referred to Kn = 0.01 |
| `checks.enabled`, `max_cell_over_mfp`, `max_courant`, `warn_courant`, `min_particles_per_cell` | true, 1.0, 1.0, 0.5, 5.0 | `max_courant` is a hard error; 0 < `warn_courant` ≤ `max_courant` |
| `checks.min_nozzle_faces`, `max_nozzle_area_error` | 12, 0.05 | hard errors |
| `output.dialect`, `inflow_model` | `standard`, `plumeFieldInflow` | the only values accepted |
| `post.centerline_x_over_D_max`, `centerline_points`, `centerline_radius_over_D` | 10.0, 200, 0.25 | the committed template uses 0.10 |
| `post.plane`, `plane_half_over_D`, `contour_levels`, `plane_points_x`, `plane_points_lateral` | `xz`, 5.0, (0.1, 0.01, 0.001), 240, 240 | `plane`: `xz` or `xy` |

### cai2012 `study.yaml`

| Key | Default | Note |
|---|---|---|
| `base_case`, `cases_dir`, `manifest_name` | `baseCase`, `Cases`, `manifest.yaml` | |
| `kn_cases` | — | list of `name`, `Kn`, `enabled` (true), `overrides` (a partial `case.yaml`), `note` |
| `particle_levels` | — | list of `name`, `multiplier`, `enabled`, `note`; a 1× level is **required** |
| `run_time_multiplier`, `output_frequency_multiplier` | 1.0, 1.0 | study-wide; scale `endTime` (not the averaging start) and divide the write interval |
| `meta` | — | free-form; in the Kn extension, `meta.knudsen_suite` selects the suite |

`overrides` may **not** set `exit.knudsen`, `model`,
`dsmc.numerical_particle_multiplier`, `dsmc.run_time_multiplier` or
`dsmc.output_frequency_multiplier` (`study.FORBIDDEN_OVERRIDES`). Cases are named
`<KnCase>_<level>`.

## cai2012-health

`study.yaml` (`plumetools.cai2012.health`): `base_case`, `cases_dir`,
`manifest_name`, `weights` (list of `name`, `target_particles_per_cell`,
`enabled`, `note`), `samplings` (list of `name`, `sampling_domain_transits`,
`enabled`, `note`), `meta`. There is **no** `overrides` mechanism; the two keys
written into each case are `resolution.target_particles_per_cell` and
`dsmc.sampling_domain_transits`. The rest is `baseCase/case.yaml`, a Cai
`case.yaml` at `exit.knudsen: 100`. `results/cost-model.yaml` holds measured
wall clock per case (`per_parcel_s`, `per_cell_s`, the measurements); it is
read, not configured.

## cai2012-kn-extension

`study.yaml` uses the Cai study schema above. Additional policy
(`plumetools.cai2012ext.policy`):

- `meta.knudsen_suite: full_decades` requires exactly the seven decade cases
  100 … 0.0001; **no declaration** keeps the historical four (1e-4, 1e-3, 1,
  10) and refuses Cai's 100, 0.1 and 0.01 even when disabled. Unknown profiles
  are refused.
- One particle level (1×) only.
- Per-Kn `overrides` may change the mesh, time step and parcel target; gas,
  geometry, speed ratio, startup and averaging stay fixed.
- `baseCase/case.yaml` must equal `cases/cai2012-health/baseCase/case.yaml`
  key by key except the keys in `policy.DELIBERATE_OVERRIDES`, each with its
  reason in the code (for example `exit.knudsen`,
  `resolution.target_particles_per_cell` = 40,
  `dsmc.sampling_domain_transits` = 4.5).

### `interface.yaml` and `interfaces.yaml`

`interfaces.yaml` (the track): `overlay` (`baseCaseInterface`), `tree`
(`Interface`), `manifest_name`, `resources` (`bytes_per_cell`,
`bytes_per_parcel`, `memory_ceiling_bytes`, `memory_ceiling_pbs`), and the
groups `production` (list: `name`, `kn_case`, `interface`, `overrides`,
`note`), `validation` (`knudsen`, `name_prefix`, `bird_thresholds`, `interface`,
`overrides`, `overlap`) and `smoke` (list: `name`, `knudsen`, `interface`,
`overrides`).

`interface:` blocks — the per-case `interface.yaml` — are
`cai2012ext.interface_config.InterfaceConfig`:

| Key | Default | Allowed |
|---|---|---|
| `enabled` | false | (true in generated interface cases) |
| `purpose` | `production` | `production`, `interface-validation`, `development-smoke` |
| `run_mode` | `production` | `production`, `exploratory` |
| `surface` | `sphere` | `sphere` |
| `bird_definition` | `general` | `general`, `radial_approx` |
| `bird_threshold` | 0.05 | |
| `core_half_angle_deg` | null | a decision the track leaves to you |
| `core_angle_provisional`, `source_validity_provisional` | false | mark a value as provisional |
| `source_validity_min_radius_over_D` | null | a decision the track leaves to you |
| `source_offset_over_D` | (0, 0, 0) | |
| `angular_exponent_offset` | 0.41 | |
| `tails` | `exclude` | `exclude`, `inject_exploratory` |
| `radius_scan_points`, `theta_scan_points` | 400, 91 | |
| `domain_clearance_over_D`, `diagnostic_radius_over_D`, `diagnostic_clearance_over_D` | 0.5, 10.0, 0.5 | |
| `meshed_flux_tolerance` | 0.02 | |
| `validation_evidence` | null | path to a passing `overlap.yaml` — provenance, so a change is *stale*, not a conflict |
| `label` | `''` | provenance |

The `overlap:` block (validation group): `baseline_case`, `region`
(`r_min_over_D`, `r_max_over_D`, `theta_max_deg`, `r_bin_over_D`,
`theta_bin_deg`, `interface_buffer_over_D`), `targets` (`number_density`,
`speed`, `temperature` — project targets, not published ones), `uncertainty`
(`blocks`, `confidence`), `low_signal_relative_se`.

### `scaling.yaml`

`schema_version` (1), `registration` (optional path), `output`
(`results/scaling/full-suite`), `resources` (`cores` ≤ 8, `max_input_bytes`,
`max_total_input_bytes`, `chunk_size`, `smoke_max_points_per_case`,
`smoke_max_nfev`), `fit` (`tracks` F/C, `subsets` exploratory/suitable,
`families` baseline/axial/axial_transverse/radial_angular, `degree`,
`angular_order`, `feature` `inverse_kn`, bounds, `quadrature_order`,
`regularization`, `floor`, angle edges, `radial_bins`, `max_nfev`,
`holdout_kn`, `block_resamples`), `evaluation` (`n_radial`, `n_azimuthal`,
`temperature_transfer`, `angular_radii_D`, `sensitivity`).
`scaling-registration.yaml`: `scope`, `purpose`
(`analysis-and-specification-only`), `enabled`, `allow_paper_decades`,
`prepared_decades`, `cases` (list: `name`, `path` inside the study, `Kn`,
optional `field_npz`, `note`). The portable NPZ schema is in
[the scaling explanation](../cai2012-scaling.md#inventory-registration-and-source-compatibility).

## Rarefaction settings

The diagnostics' thresholds (`plumetools.rarefaction.settings.DiagnosticSettings`)
appear as `diagnostics:` in a nozzle-plume `case.yaml`, as the Kn extension's
`baseCase/diagnostics.yaml`, and inside the file given to
`python -m plumetools.rarefaction --settings`:

| Key | Default | Meaning |
|---|---|---|
| `enabled` | true | |
| `min_sampled_particles` | 20 | a cell is valid above this |
| `particle_count_basis` | `per_step` | `per_step` or `accumulated` (a much weaker test) |
| `cells_per_mean_free_path` | 3.0 | the resolution criterion λ/Δ ≥ 3 |
| `practical_min_cell_size_m` | 0.01 | the finest cell considered buildable |
| `mesh_budget_cells` | 10 000 000 | budget for the "cells needed" estimate |
| `reference_lengths_m` | `{}` | extra Kn reference lengths; `D` is reserved |
| `select` | null | which diagnostics to compute |
| `field_names` | `{}` | overrides of the input field names (`rhoNMean`, `rhoMMean`, `momentumMean`, `linearKEMean`, `dsmcRhoNMean`, `C`, `V`) |
| `invalid_value` | −1.0 | written where a diagnostic is undefined |
| `n_worst_cells`, `worst_cells_separation` | 10, 4.0 | the summary's worst-violation list |

The stand-alone settings file has exactly three top-level keys: `gas`
(`molecular_mass_kg` and `molecular_diameter_m` required; `omega` 0.5,
`t_ref_K` 273, `mean_free_path_convention` `vhs`, `species_name`),
`diagnostics` (the table above) and `background_cell_size_m`. Details:
[Rarefaction diagnostics](../rarefaction-diagnostics.md).

## `assessment.yaml`

`families` (`cai2012`, `cai2012-kn-extension`), `kn` (the seven values),
`source_case` (supplies D, S₀, T₀ and the gas), `smoke` and `full` (`kn`,
`samples`, `hazard_steps`, `partner_radial`, `partner_azimuthal`,
`partner_speed`, `batches`, `seed`, `domain_length_D`, `domain_radius_D`,
`chunk_size`), `assessment_bins` (`x_count`, `r_count`, `x_max_over_D`,
`r_max_over_D`, `density_floor_over_n0`).

## `hybrid.yaml`

`schema: hybridplume/1`, `default_profile`, `base`, `profiles`. A profile is a
partial mapping deep-merged onto `base`; `setup` writes the merge to
`<case>/hybrid.resolved.yaml`. Required after merging: `source.diameter_m`,
`speed_ratio`, `temperature_K`, `knudsen`; `gas.mass_kg`, `diameter_m`,
`t_ref_K`, `omega`; `domain.x_max_over_D`, `half_width_over_D`;
`mesh.core_cell_over_D`; `dsmc.backend`.

| Constraint | Rule |
|---|---|
| `dsmc.backend` | `dsmcFoam` only |
| `interface.mode` | `automatic`, `manual`, `none` (full DSMC) |
| `interface.manual.type` | `sphere`, `cylinder`, `box` |
| forbidden keys | `mesh.axisymmetric`, `mesh.wedge`, `dsmc.axisymmetric`, `dsmc.radial_weighting`, `domain.symmetry` — the family is full 3-D |
| `domain.half_width_over_D`, `mesh.core_half_over_D` | > 0.5 |
| `cfd.background_pressure_ratio` | in (0, 0.1) — a numerical floor, not an ambient gas |
| `source.knudsen` outside Cai's 100, 0.1, 0.01 | allowed, reported as an extension |

Section contents are listed in the [hybrid tutorial](../tutorials/cai2012-hybrid.md#inputs-and-profiles).

## nozzleplume `case.yaml` (r-4d)

`plumetools.nozzleplume.config`; `model: nozzle_exit_plume`. Required sections
`nozzle`, `exit`, `gas`, `domain`, `mesh`; optional `plate` (absent or `null` =
no plate), `dsmc`, `checks`, `post`, `diagnostics`, `execution_profile`
(`smoke`/`production`), `meta`.

| Key | Default | Constraint |
|---|---|---|
| `nozzle.exit_diameter_m` | required | |
| `exit.temperature_K` | required | |
| `exit.number_density_per_m3` / `mass_density_kg_per_m3` / `knudsen` | — | **exactly one** |
| `exit.speed_m_per_s` / `speed_ratio` | — | **exactly one** |
| `gas.molecular_mass_kg`, `molecular_diameter_m` | required | |
| `gas.species_name`, `omega`, `t_ref_K` | `gas`, 0.5, 273.0 | `omega` in [0.5, 1.0] |
| `gas.internal_degrees_of_freedom`, `mean_free_path_convention` | 0, `vhs` | whole number; `vhs` or `hard_sphere` |
| `plate.standoff_m`, `thickness_m`, `size_y_m`, `size_z_m` | required if a plate | |
| `plate.centre_y_m`, `centre_z_m`, `temperature_K` | 0, 0, 300 | |
| `domain.x_max_m`, `y_half_m`, `z_half_m` | required | |
| `domain.symmetry_planes` | `()` | subset of `y`, `z`, no repeats; a plate must stay centred on a symmetry plane |
| `mesh.background_cell_size_m` | required | the background is cubic |
| `mesh.plate_refinement_level`, `n_cells_between_levels`, `max_cells` | 2, 3, 3 000 000 | |
| `mesh.refinement_regions` | `()` | list: `name` (unique, not a patch name), `level`, `x_max_m`, `shape` (`cone` default, or `box`), `x_min_m` (0); cone: `radius1_m`, `radius2_m`, optional `inner_radius1_m`, `inner_radius2_m`; box: `y_min_m`, `y_max_m`, `z_min_m`, `z_max_m` |
| `mesh.patch_names` | `nozzle`, `upstreamVacuum`, `vacuum`, `plate`, `symmetry` | distinct |
| `dsmc.n_equivalent_particles` | null | written exactly as given; null = derived for `target_particles_per_exit_cell` (20) |
| `dsmc.collision_model` | `VariableHardSphere` | `VariableHardSphere`, `LarsenBorgnakkeVariableHardSphere`, `NoBinaryCollision` (written as `none`) |
| `dsmc.relaxation_collision_number` | 5.0 | |
| `dsmc.wall_interaction_model`, `diffuse_fraction` | `MaxwellianThermal`, 1.0 | `MaxwellianThermal`, `SpecularReflection`, `MixedDiffuseSpecular`; fraction in [0, 1] |
| `dsmc.delta_t_s`, `courant_target` | null, 0.2 | |
| `dsmc.end_time_s`, `average_start_s`, `write_interval_s` | null | derived from transits when null |
| `dsmc.transient_domain_transits`, `sampling_domain_transits` | 4.0, 2.0 | |
| `dsmc.initial_number_density_per_m3`, `initial_temperature_K` | 1e10, 300 | |
| `dsmc.n_subdomains`, `decomposition_method` | 4, `scotch` | `scotch` or `simple` |
| `checks.*` | `max_cell_over_mfp` 1, `max_time_step_over_collision_time` 1, `max_courant` 1, `warn_courant` 0.5, `min_particles_per_cell` 5, `min_nozzle_faces` 12, `max_nozzle_area_error` 0.05, `warn_parcels` 3e7 | |
| `post.centreline_radius_over_D`, `centreline_points` | 0.25, 200 | points ≥ 2 |

## Slice-imagery specs

`plumetools.viz.spec` loads `viz.yaml` files and the package default
`plumetools/viz/slices.yaml`. Sections: `sampling` (`time`: `latest`, `first`,
`all` or a value; `time_min`, `time_max`, `prefer_mean` (true),
`mesh_regions`, `field_type` `CELLS`/`POINTS`, `skip_constant_fields`,
`length_scale_m`), `image` (size, zoom, margins, offsets, `preset`, colour-bar
settings, background), `output` (`directory` `results/viz`, `subdirectory`,
`filename` `{field}_{plane}`, `log_directory`, `manifest`), `video`
(`enabled`, `framerate` 4, `codec` `libx264`, `container` `mp4`, `quality` 18,
`min_frames`, `ffmpeg`), `planes` (list: `name`, `normal`, `camera_up`,
`origin` or `origin_over_L`, per-plane image overrides), `fields` (list:
`name`, `enabled`, `kind`, `planes`, `prefer_mean`, `component` (`''`,
`Magnitude`, `X`, `Y`, `Z`), `log`, `log_decades` (6), `range`, `preset`,
`legend_title`, `color_bar`, `field_type`), `derived` (list: `name`,
`expression`, `kind`, `component`, `units`, `log`). What each setting does and
how it can mislead: [Slice imagery](../viz-slices.md).

## Custom OpenFOAM libraries

The dictionary contracts the generated cases rely on. The generators write
these entries; they are listed so a hand-built case can use the libraries.

**`plumeFieldInflow`** (`libplumeDsmcBoundaryModels.so`):

```
constant/dsmcProperties
    InflowBoundaryModel   plumeFieldInflow;
    plumeFieldInflowCoeffs
    {
        patches ( nozzle );                         // explicit list; a 'wall' patch is refused
        numberDensityFields { Ar  boundaryNumberDensity_Ar; }   // a field per species
    }
system/controlDict
    libs ( "libplumeDsmcBoundaryModels.so" );
```

The named field is a `volScalarField` whose boundary values on the listed
patches carry the number density in m⁻³; `boundaryU` and `boundaryT` carry the
per-face velocity and temperature as for `FreeStream`. A missing library is a
hard failure. Behaviour, failure modes, restarts and parallel runs:
[plumeFieldInflow](../plume-field-inflow.md).

**`hybridPlumeThermo`** (`libhybridPlumeThermo.so`), used by the hybrid
family's `rhoCentralFoam` precursor:

```
constant/thermophysicalProperties
    thermoType
    {
        type            hePsiThermo;
        mixture         pureMixture;
        transport       tabulated;
        thermo          hConst;
        equationOfState perfectGas;
        specie          specie;
        energy          sensibleInternalEnergy;
    }
system/controlDict
    libs ( "libhybridPlumeThermo.so" );
```

It registers OpenFOAM's own `tabulatedTransport` for `psiThermo`, which v2512
compiles only for `rhoThermo`; the μ(T) and κ(T) tables are written by
`plumetools.hybridplume.cfd` from the same gas object as the DSMC dictionaries.
Build each library in the environment that runs the case.

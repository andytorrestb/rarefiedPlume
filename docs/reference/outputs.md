# Outputs

What each family writes, where, with which columns and keys, and which command
consumes it. Paths are relative to the case or study directory named. Read the
[quantity distinctions](#quantities-that-are-easy-to-confuse) first: most
misreadings of these files come from them.

## Quantities that are easy to confuse

| Pair | Here |
|---|---|
| **number density vs parcel count** | `rhoN` (and `rhoNMean`) is the real number density n [m⁻³]. `dsmcRhoN` is the **parcel count in the cell** — dimensionless, not a density: `DSMCCloud::calculateFields` adds one per parcel and never divides by the volume. Measured on Cai Kn100: `rhoN · V / dsmcRhoN` equals `nEquivalentParticles` to 1.5e-9 |
| **number vs mass density** | `rhoM` = m n [kg/m³]; for one species the same picture as `rhoN` in other units |
| **dimensional vs normalised velocity** | the solver writes momentum density (`momentum`, kg m⁻² s⁻¹); velocity is `momentum / rhoM` [m/s]. Cai files use U₁√β₀ (dimensionless, β₀ = 1/(2RT₀)); nozzle-plume files use U/U_exit |
| **scalar vs directional temperature** | `Ttra` and the T in the comparison files are the scalar translational temperature from the averaged moments, 3RT = ⟨v²⟩ − \|u\|². In an expanding plume the parallel and perpendicular temperatures differ (the hybrid collisionless smoke recorded T∥/T⊥ from 0.94 to 3.75 on the axis); a scalar agreement does not establish tensor agreement |
| **instantaneous vs averaged** | `X` is one step's sample (shot noise); `XMean` is the `fieldAverage` running mean from `average_start_s`; `XPrime2Mean` its variance. Post-processing reads the means |
| **gas pressure vs wall pressure** | wall pressure is the wall-normal component of `fDMean` (momentum delivered per area per time), not n k T of the adjacent cell |
| **`boundaryT` type** | a `volScalarField` in the standard dialect, a `volVectorField (T 0 0)` in the MNF dialect; writing the wrong one is a read failure |

## OpenFOAM fields a dsmcFoam run writes

From `plumetools/viz/catalog.py`, which states what each field is:

| Field | Units | Meaning |
|---|---|---|
| `rhoN` | m⁻³ | number density |
| `rhoM` | kg/m³ | mass density |
| `momentum` | kg m⁻² s⁻¹ | momentum density ρu |
| `linearKE` | J/m³ | translational kinetic-energy density ½ρ⟨v²⟩ |
| `internalE`, `iDof` | J/m³, m⁻³ | internal-energy and internal-DOF density; zero for a monatomic gas |
| `dsmcRhoN` | parcels/cell | parcel count |
| `dsmcSigmaTcRMax` | m³/s | the NTC collision-selection ceiling |
| `q`, `fD` | W/m², Pa | heat flux and force density on boundaries; zero unless the case has a wall |
| `boundaryT`, `boundaryU`, `boundaryNumberDensity_<sp>` | K, m/s, m⁻³ | the inflow state per boundary face (written into `0/` by the generators, then AUTO_WRITE in every time) |
| `C`, `V` | m, m³ | cell centres and volumes, from `postProcess -func writeCellCentres`/`writeCellVolumes` (the `Allpost` scripts run them) |

The `Mean` and `Prime2Mean` variants come from the generated `fieldAverage`
function object; `uniform/functionObjects/` in each time directory holds its
accumulators, which is why averaging continues across a restart.

## Shared

| Path | Writer | Contents | Consumers |
|---|---|---|---|
| `<case>/<time>/mfp`, `KnD`, `KnCellMax`, `KnCellV`, `cellOverMfpMax`, `cellOverMfpV`, `KnGradN`, `KnGradT`, `cellSizeTarget`, `cellSizeMax`, `cellSizeV`, `parcelsPerStep`, `sampledParcels`, `mfpResolved`, `diagValid`, `diagValidGradN`, `diagValidGradT`, `diagStatus` | `plumetools.rarefaction` (nozzle plume `Allpost`, Kn extension `Alldiagnose.py`, `python -m plumetools.rarefaction`) | per-cell diagnostics; `invalid_value` (−1) where undefined; `mfpResolved` is 1 where λ/Δ_max ≥ 3, 0 where not, −1 no data | ParaView; `viz-diagnostics.yaml` |
| `results/diagnostics_summary.yaml` | the same | sampling window, valid fractions, Δ/λ statistics, worst violations by location, a row per cell size, permitted octree levels, a lower-bound cell count for a mesh meeting λ/3 | the remeshing decision ([Rarefaction diagnostics](../rarefaction-diagnostics.md)) |
| `results/viz/*.png`, `results/viz/manifest.yaml` | `render_slices.py` | images; the manifest records time and colour range per image | `video.py`, `dissect.py` |
| `<new>/iteration.yaml`, `<new>/provenance/` | `plumetools.iteration` | parent case, note, the parent's `case.yaml` | `diff provenance/parent-case.yaml case.yaml` |

## cai2012

| Path | Contents |
|---|---|
| `manifest.yaml` | per case: `cai_case`, particle level and multiplier, baseline weight, Kn, D, λ₀, T₀, U₀, n₀, weight, minimum cell, `deltaT`, `mesh_id` (equal across a Kn group), `startTime`/`endTime`/`writeControl`/`writeInterval` and the multipliers, the `fieldAverage` set-up |
| `Cases/<c>/results/centerline.csv` | `x_over_D`, `n_over_n0`, `n_over_n0_analytical`, `n_over_n0_analytical_axis`, `U_sqrt_beta0`, `U_sqrt_beta0_analytical`, `U_sqrt_beta0_analytical_axis`, `T_over_T0`, `T_over_T0_analytical`, `T_over_T0_analytical_axis`, `n_cells`. `_analytical` is Cai's solution averaged over the same cells with the same weights as the DSMC sample (removing the tube-averaging bias); `_analytical_axis` is exactly on the axis |
| `Cases/<c>/results/density_plane.csv` | `x_over_D`, `<lateral>_over_D` (`z` for `post.plane: xz`), `n_over_n0` at cell centres of the centre plane — no resampling |
| `Cases/<c>/results/metrics.yaml` | `case` (name, Kn, n₀, cell/λ₀, cell count, time, Cai's reported error) and `centerline` (max, mean, RMS relative errors per quantity, the index of the worst point) |
| `Cases/<c>/results/*.png` | `centerline-density.png`, `centerline-velocity.png`, `centerline-temperature.png`, `density-contour.png` |
| `results/study-table.csv` | `case, Kn, n0_per_m3, cell_over_mfp, n_cells, time, density_max_rel_error_percent, density_mean_rel_error_percent, density_rms_rel_error_percent, cai_reported_max_percent, velocity_max_rel_error_percent, temperature_max_rel_error_percent` |

`cai_reported_max_percent` is Cai's published error, carried as a reference
value, never a tolerance.

## cai2012-health

Per case, besides cai2012's files: `results/occupancy.yaml` (measured parcels
per cell per region — exit cells, inside Cai's lowest contour, the centreline
tube, the whole domain — against the floor and target) and
`results/convergence.yaml` (centreline error at every sampled frame), from
`Allaudit.py`. Study level (`analyse.py`): `results/contact-final.png`,
`contact-common.png`, `convergence.yaml` (PSNR per case, `saturated` where two
frames are identical), `sweep-table.csv`, `sweep-error.png`, `overlap.yaml`
(field-by-field byte or numeric comparison of a short run against its long
run). `results/cost-model.yaml` is tracked measured data.

## cai2012-kn-extension

| Path | Contents |
|---|---|
| `manifest-full-suite.yaml` | per case: everything cai2012's manifest records, plus `feasibility` (`status` `production-candidate` or `unsuitable-for-production`, `resolution.cell_over_mean_free_path`, `time_step.dt_over_exit_collision_time`, `cost.estimated_parcels` and memory, `unresolved` reasons), `production_status`, the suite profile, the inheritance report, git provenance, stray directories |
| `<case>/NON-PRODUCTION` | the unresolved criteria of an unsuitable case |
| `<case>/results/comparison.yaml`, `centerline_extended.csv`, `meridional_plane.csv` | `Allcompare.py`: r/D = 10 angular density and mass flux, off-axis axial and radial velocity, temperature and pressure, with exclusion counts |
| `<case>/results/diagnostics_summary.yaml`, `centreline_diagnostics.csv`, `centreline-diagnostics.png` | `Alldiagnose.py` |
| `<case>/results/occupancy.yaml`, `convergence.yaml` | `Allaudit.py` |
| `<case>/results/stationarity.yaml` | `Allconverge.py`: far-field stationarity across frames |
| `results/full-suite/study-table.csv`, `study-table.md`, `knudsen-comparison.png` | `analyse.py`; every case, including any not run (its state says so) |
| `results/full-suite-preparation/smoke.json` | `smoke_cases.py` |
| `interface-manifest.yaml`, `<case>/results/interface-report.yaml`, `interface-identity.yaml`, `interface-inflow.yaml`, `overlap.yaml`, `interface-evidence.yaml` | the interface track ([tutorial](../tutorials/cai2012-interfaces.md#outputs)) |
| `results/scaling/full-suite/` | `inventory.{json,csv}`, `assessment.json`, `stratum*/…/predictor.json`, `README.md`, figures; `smoke/` beneath it for the synthetic mode. JSON is strict: non-finite values are `null` |

Normalisations and estimator definitions of every compared quantity:
[the extension's §7–8](../cai2012-kn-extension.md#7-what-is-compared-and-how).

## cai2012-virtual-collision

In `--output`: `inventory.json`, `inventory.csv`, `skipped.json`; for
`smoke`/`full` also `report.json`, `cost.csv`, `convergence.csv`,
`model_profiles.csv`, `assessment.csv`, `assessment_skips.csv`,
`reference_assessment.json`, `reference_skips.json`, `verification.json`,
model arrays `Kn<k>_*.npz` with `*_diagnostics.json`. The committed
`results/smoke/` is an example of the full set. `compare_exports.py` writes per
case into the Kn extension's `<case>/results/virtual_collision/`.

## cai2012-hybrid

Under `<case>/results/` — one JSON per step: `setup.json`,
`mesh_check_cfd.json`, `cfd_report_cfd.json`, `breakdown.json`,
`interface.json`, `dsmc_setup.json`, `dsmc_fields.json`, `dsmc_post.json`,
`crossing_window.json`, `reference.json`, `overlap.json`, `hpc_scripts.json`;
`summary.json` and PNG figures from `report`. `<case>/hybrid.resolved.yaml` is
the merged configuration every step reads. Field-by-field meaning:
[hybrid tutorial](../tutorials/cai2012-hybrid.md#outputs).

## markelov1999

| Path | Contents |
|---|---|
| `manifest.yaml` | every case and its inputs, every `[ASSUMPTION]` |
| `Cases/gap06in/<c>/case-summary.json` | `case`; from `runInflow.py` the inflow, flux verification and DSMC checks; from `Allpost` `mesh`, `particles`, `pressure` (`windows`, one entry each), `resolution_audit` (cells, mean/median/max parcels per cell, occupied cells, total, `region_occupancy`, `target`, `meets_target`) |
| `results/study-table.csv` | one row per case |
| `results/windward_pressure.png`, `leeward_pressure.png`, `pressure_ratio.png` | plots |

## nozzleplume (r-4d)

| Path | Columns / keys |
|---|---|
| `results/plate_loads.csv` | `x_m, y_m, z_m, area_m2, side` (+1 upstream face, −1 lee face, 0 edge), `pressure_pa, shear_pa, heat_flux_w_per_m2, fD_x_pa, fD_y_pa, fD_z_pa` |
| `results/plate_summary.yaml` | per side: face count, area, force, peak pressure and its location; totals `force_on_plate_n`, `force_on_modelled_part_n` (a symmetric quarter or half), `heat_load_w`, `jet_thrust_n`, `jet_power_w`, `heat_load_over_jet_power` |
| `results/centreline.csv` | `x_m, x_over_D, n_over_n_exit, U_over_U_exit, T_over_T_exit, n_cells, cell_size_m` |
| `results/plate-pressure.png`, `plate-heat-flux.png`, `centreline.png` | figures |
| diagnostics, `results/viz/` | as under *Shared* |

## 3d-inflow

`0/boundaryU`, `0/boundaryT`, `0/boundaryNumberDensity_Ar` and, in the
standard dialect, nine zeroed measurement fields. No post-processing outputs
(RP-17). The regression goldens (`tests/regression/golden/`) are the
pre-refactor values: `3d_inflow_v0.npz` (`labels`, `rhoN`, `U`, `T` per inflow
face, in insertion order) and the three `boundary*.txt` files.

## Run bookkeeping

| Path | Meaning |
|---|---|
| `log.<application>` | one log per OpenFOAM tool; `log.dsmcFoam` holds `Number of dsmc particles` per output step — check it has levelled off |
| `processor*/` | decomposed state; what a parallel resume reads |
| `<time>/uniform/time` | written by `Time` on every rank at every write; part of the checkpoint completeness test |
| `system/decomposeParDict.Allrun` | the per-run decomposition override written by `-np` |
| `output.log`, `error.log` | PBS job output (Kn extension) |

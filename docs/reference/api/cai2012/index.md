# `plumetools.cai2012`

The Cai & Wang 2012 collisionless circular-plume family: its own `case.yaml`
schema, argon VHS gas, the re-derived collisionless solution, the graded
Cartesian mesh, the study and manifest, DSMC quality checks, post-processing,
the statistical-health sweep and the post-run audit. Its generated output is
pinned byte for byte, so `cai2012ext` builds on it without editing it, and
`nozzleplume` re-implements its exit state rather than importing it.

Depends on: `numpy`, `scipy` (analytical quadrature), `PyYAML`, `matplotlib`
(figures in `post`); `plumetools.foamio`, `plumetools.mesh`. Used by: the
[cai2012](../../../tutorials/cai2012.md), [cai2012-health](../../../tutorials/cai2012-health.md)
and [Kn-extension](../../../tutorials/cai2012-kn-extension.md) workflows.

::: plumetools.cai2012
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - PAPER
        - MODEL
        - CAI_REPORTED_MAX_DENSITY_ERROR_PERCENT

## Modules

| Module | Purpose | Principal interfaces | Workflow step |
|---|---|---|---|
| [`plumetools.cai2012.config`](config.md) | `case.yaml` schema, loading, validation | `load_case_config`, `CaiCaseConfig`, `CaiConfigError` | every step |
| [`plumetools.cai2012.constants`](constants.md) | physical constants | `BOLTZMANN_J_PER_K`, `AVOGADRO_PER_MOL` | — |
| [`plumetools.cai2012.gas`](gas.md) | argon VHS properties; Kn ↔ λ₀ ↔ n₀ | `VhsSpecies`, `ARGON`, `mean_free_path`, `number_density_from_kn`, `knudsen_from_number_density` | derivation |
| [`plumetools.cai2012.geometry`](geometry.md) | exit disk at x = 0 and the vacuum box | `CaiGeometry`, `from_config` | derivation |
| [`plumetools.cai2012.inflow`](inflow.md) | the derived exit state and run settings | `ExitState`, `from_config`, `RunSettings`, `derive_run_settings` | derivation, `runInflow.py` |
| [`plumetools.cai2012.analytical`](analytical.md) | Cai's collisionless solution, re-derived | `centerline_density_ratio`, `centerline_axial_speed_ratio`, `centerline_temperature_ratio`, `centerline_profile`, `density_ratio`, `moments`, `describe` | post-processing |
| [`plumetools.cai2012.mesh`](mesh.md) | graded multi-block mesh plan, uniform-grid cost, nozzle-face prediction | `plan`, `MeshPlan`, `uniform_cost`, `predicted_nozzle_faces`, `estimate_particles` | `Allmesh` |
| [`plumetools.cai2012.checks`](checks.md) | mesh, time-step and boundary checks | `run_checks`, `CheckReport`, `DsmcCheckError` | `Allmesh`, generation |
| [`plumetools.cai2012.dictionaries`](dictionaries.md) | `dsmcProperties`, `controlDict`, … | `write_case_dictionaries`, `render_control_dict`, `simple_split`, `AVERAGED_FIELDS` | `Allmesh` |
| [`plumetools.cai2012.foamfields`](foamfields.md) | the `0/` fields | `write_all`, `write_inflow_fields`, `write_measurement_fields` | `runInflow.py` |
| [`plumetools.cai2012.study`](study.md) | `study.yaml`, the Kn × particle matrix, the manifest | `load_study`, `StudyConfig`, `apply_case_parameters`, `write_manifest`, `FORBIDDEN_OVERRIDES` | `generate_cases.py` |
| [`plumetools.cai2012.verify`](verify.md) | check a generated tree against its study | `verify_tree`, `verify_group`, `verify_schedule`, `VerificationError` | `validate_cases.py` |
| [`plumetools.cai2012.post`](post.md) | the four validation outputs and their scores | `read_case`, `centerline`, `density_plane`, `error_metrics`, `write_results`, `study_row`, `STUDY_COLUMNS` | `Allpost` |
| [`plumetools.cai2012.health`](health.md) | the 3 × 3 statistical-health sweep and its cost model | `load_study`, `clone_cases`, `fit_cost_model`, `cost_report` | health `generate_cases.py`, `AllrunCases` |
| [`plumetools.cai2012.audit`](audit.md) | measured occupancy, determinism, the sweep table | `occupancy_audit`, `region_occupancy`, `compare_field`, `weight_consistency` | `Allaudit.py`, `analyse.py` |

## The derivation chain, in code

Only `exit.knudsen` varies per case. Everything the solver runs at follows:

<!-- docs-example: cai2012-derivation -->
```python
from pathlib import Path
from plumetools.cai2012 import gas, inflow
from plumetools.cai2012.config import load_case_config

cfg = load_case_config(Path("cases/cai2012/baseCase"))      # validated; rejects unknown keys
exit_state = inflow.from_config(cfg)                          # lambda0 = Kn * D, then n0, then U0
assert exit_state.mean_free_path_m == cfg.exit.knudsen * cfg.nozzle.diameter_m

# Kn -> n0 and back, for argon at 300 K with the VHS convention
n0 = gas.number_density_from_kn(0.1, 0.2, gas.ARGON, 300.0)
assert abs(gas.knudsen_from_number_density(n0, 0.2, gas.ARGON, 300.0) - 0.1) < 1e-12
```

Units: lengths m, densities m⁻³, speeds m/s, temperatures K. `ExitState`
carries `mean_free_path_m`, `number_density_per_m3`, `velocity_m_per_s`,
`speed_ratio` and the temperature; `derive_run_settings(cfg, exit_state, geom,
min_cell_size_m=…, exit_cell_volume_m3=…)` adds the particle weight, `deltaT`
and the schedule, and needs the mesh plan's smallest cell and exit-cell volume
(`mesh.plan(cfg, geom, exit_state)`).

## The analytical solution

`analytical` evaluates Cai's collisionless field for a uniform drifting
Maxwellian on the disk r ≤ R₀ at x = 0. The centreline has closed forms
(C0–C4 in the module docstring); off the axis the density is a Gauss–Legendre
quadrature over the disk (C5). Inputs are dimensional (`x`, `y`, `z` in m,
`nozzle_radius_m`), `x` may be an array, and results are the normalised
quantities n/n₀, U₁√β₀ and T/T₀:

$$
\frac{n}{n_0}\Big|_{x\to 0} = \tfrac12\left(1 + \operatorname{erf} S_0\right),
\qquad t_0 = \frac{x}{\sqrt{x^2 + R_0^2}} .
$$

<!-- docs-example: cai2012-analytical -->
```python
import math
from plumetools.cai2012 import analytical

R0, S0 = 0.1, 2.0
near_exit = analytical.centerline_density_ratio(1e-9, nozzle_radius_m=R0, speed_ratio=S0)
assert abs(near_exit - 0.5 * (1 + math.erf(S0))) < 1e-9       # not 1: only the outgoing half is present

profile = analytical.centerline_profile([0.1, 0.5, 1.0], nozzle_radius_m=R0, speed_ratio=S0)
# keys: x_m, x_over_D, n_over_n0, U_sqrt_beta0, T_over_T0

# the off-axis quadrature agrees with the closed form on the axis
on_axis = analytical.density_ratio(0.2, 0.0, 0.0, nozzle_radius_m=R0, speed_ratio=S0)
assert abs(on_axis - analytical.centerline_density_ratio(0.2, R0, S0)) < 1e-9
```

`tests/unit/test_cai2012_analytical.py` asserts the limits listed in the
module docstring and the agreement of (C5) with (C1) on the axis to 1e-10.
What those tests do and do not establish: [Evidence](../../../explanation/evidence.md).

## Post-processing contract

`post.read_case` reads the **time-averaged** fields (`REQUIRED_FIELDS`) with
the cell centres and volumes; `centerline` averages the cells within
`post.centerline_radius_over_D` of the axis, volume-weighted, and averages the
analytical solution over the same cells (`*_analytical` columns) so the tube
bias cancels; `error_metrics` scores max/mean/RMS relative error. Readers raise
`post.PostError`, a subclass of `foamio.errors.FoamDataError`. Column names:
[Outputs](../../outputs.md#cai2012).

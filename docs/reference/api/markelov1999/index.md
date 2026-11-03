# `plumetools.markelov1999`

The AIAA 99-3455 flat-plate-in-a-cylinder-wake family: a source-flow model as
printed in the paper (deliberately separate from the frozen legacy
`plumetools.sourceflow`), the geometry and its paper-to-SI conversions, the
snappyHexMesh setup, inflow-flux verification, particle-resolution estimate and
audit, DSMC checks, dictionaries, `0/` fields, the CaseFoam-based study, and
post-processing of wall pressures.

This package is **frozen**: general code goes in new packages, and defects
found in it are reported rather than fixed (see
[Known gaps](../../../development/known-gaps.md)). Its `case.yaml` is read by
the shared legacy loader, [`plumetools.config`](../legacy/config.md), with
`model: markelov1999_axisymmetric`.

Depends on: `numpy`, `scipy`, `PyYAML`, `matplotlib`; **`casefoam`** for
`study.clone_cases`; `plumetools.config`, `plumetools.foamio`,
`plumetools.mesh`. Workflow: [markelov1999](../../../tutorials/markelov1999.md).

::: plumetools.markelov1999
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - PAPER

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.markelov1999.constants`](constants.md) | constants and unit conversions | `psi_to_pa`, `inch_to_m`, `molecular_mass_kg`, `PSI_TO_PA`, `INCH_TO_M` |
| [`plumetools.markelov1999.sourceflow`](sourceflow.md) | the source-flow model as printed (P1 …) | `limiting_angle`, `limiting_velocity`, `angular`, `normalization_coefficient`, `sonic_density_ratio`, `number_density`, `radial_velocity`, `MODEL_NAME` |
| [`plumetools.markelov1999.geometry`](geometry.md) | the configuration, paper dimensions and assumptions | `MarkelovGeometry`, `from_config`, `PAPER_DIMENSIONS_IN`, `ASSUMED_DIMENSIONS`, `report_paper_conversions` |
| [`plumetools.markelov1999.inflow`](inflow.md) | evaluate the model over the inflow patch | `compute_inflow`, `InflowResult`, `area_weighted_number_density`, `faces_outside_the_cone` |
| [`plumetools.markelov1999.flux`](flux.md) | integrated inflow verification: mass, momentum, energy | `verify`, `FluxVerification`, `analytical_face_fluxes`, `quadrature_face_fluxes` |
| [`plumetools.markelov1999.mesh`](mesh.md) | background box + snappyHexMesh carving of the cavity, cylinder and plate | `write_mesh_setup`, `refinement_levels`, `surface_cell_sizes`, `location_in_mesh`, `describe` |
| [`plumetools.markelov1999.resolution`](resolution.md) | particle-resolution estimate per region, and the post-run audit | `estimate`, `ResolutionEstimate`, `audit`, `AuditResult`, `regions_of_interest` |
| [`plumetools.markelov1999.checks`](checks.md) | DSMC mesh and time-step checks | `run_checks`, `CheckReport`, `DsmcCheckError` |
| [`plumetools.markelov1999.foamdicts`](foamdicts.md) | `dsmcProperties`, `controlDict`, … | `write_case_dictionaries`, `render_dsmc_properties`, `render_control_dict` |
| [`plumetools.markelov1999.foamfields`](foamfields.md) | the `0/` fields | `write_all`, `write_inflow_fields`, `read_patch_field` |
| [`plumetools.markelov1999.study`](study.md) | `study.yaml`, CaseFoam cloning, manifest | `load_study`, `require_casefoam`, `clone_cases`, `apply_case_parameters`, `write_manifest` |
| [`plumetools.markelov1999.verify`](verify.md) | the built mesh against the declared geometry; CLI | `verify`, `MeshVerificationError`, `main` |
| [`plumetools.markelov1999.postprocess`](postprocess.md) | mesh and particle statistics, windowed wall pressures, study table | `wall_pressure`, `pressure_report`, `PressureWindow`, `parse_solver_log`, `summary_row`, `write_study_table`, `plot_study` |

## The source-flow model

<!-- docs-example: markelov-sourceflow -->
```python
import math
from plumetools.markelov1999 import sourceflow

theta_lim = sourceflow.limiting_angle(1.4)                     # rad, for gamma = 1.4 (N2)
v_lim = sourceflow.limiting_velocity(1.4, 300.0, 4.65e-26)     # m/s at T0 = 300 K
assert 0.0 < theta_lim < math.pi and v_lim > 0.0
```

`number_density(theta, radius_m, *, gamma, T0_K, mass_kg, p0_pa,
orifice_radius_m=ORIFICE_RADIUS_M, exponent_offset=0.41, quadrature_points=0)`
returns molecules/m³ at polar angle θ (from the plume axis) and radius r from
the orifice. The six places where this model and the legacy code disagree, and
the provenance of the 0.41, are in
[The AIAA 99-3455 family §2](../../../markelov1999-case.md#2-the-source-flow-model).
The flux verification (`flux.verify`) checks the integrated mass, momentum and
energy through the meshed inflow surface two independent ways.

!!! note "A docstring to read with care"
    `resolution.mean_free_path` uses the hard-sphere expression with the VHS
    reference diameter and says so; its docstring writes the omitted VHS factor
    as (T_ref/T)^(ω−½), which is the factor on the cross-section. On λ it is
    (T/T_ref)^(ω−½), as `cai2012.gas` applies. The value returned is unaffected.

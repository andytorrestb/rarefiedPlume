# `plumetools.cai2012ext`

The Cai 2012 Knudsen extension: matrix policy, idempotent generation,
provenance, pre-run feasibility, the off-axis and angular comparisons,
stationarity, the nozzle-plume diagnostics wired to a Cai case, the PBS
wrappers — and two further tracks, the **interface** (removing the dense core
and injecting downstream) and **field scaling** (empirical coordinate maps of
Cai's density). It imports `cai2012`, `markelov1999`, `nozzleplume` and
`rarefaction` and edits none of them.

Depends on: `numpy`, `scipy`, `PyYAML`, `matplotlib` (figures);
`plumetools.cai2012`, `plumetools.nozzleplume`, `plumetools.rarefaction`,
`plumetools.markelov1999.sourceflow` (interface closure), `plumetools.foamio`.
Workflows: [direct cases](../../../tutorials/cai2012-kn-extension.md),
[interface track](../../../tutorials/cai2012-interfaces.md),
[field scaling](../../../tutorials/cai2012-scaling.md),
[virtual collisions](../../../tutorials/cai2012-virtual-collision.md) (as a data source).

::: plumetools.cai2012ext
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - EXCLUDED_KNUDSEN
        - EXTENSION_KNUDSEN
        - FULL_SUITE_KNUDSEN
        - PARTICLE_MULTIPLIERS
        - REFERENCE_KNUDSEN
        - FAMILY_DIRECTORY
        - GENERATOR

## Direct-exit cases

| Module | Purpose | Principal interfaces | Used by |
|---|---|---|---|
| [`plumetools.cai2012ext.policy`](policy.md) | suite profiles (legacy four-Kn, `full_decades`), template inheritance from the health study, allowed overrides, schedulability | `check_study`, `check_template`, `suite_name`, `suite_knudsen`, `inheritance_report`, `schedulable`, `DELIBERATE_OVERRIDES`, `ExtensionError` | generation, runners, `hpc guard` |
| [`plumetools.cai2012ext.direct_case`](direct_case.md) | one direct-exit case, derived once | `DirectCase`, `derive`, `derive_config` | generation, feasibility, `hpc label`, post |
| [`plumetools.cai2012ext.feasibility`](feasibility.md) | production criteria and cost, before meshing | `assess_case`, `assess`, `Criterion`, `load_cost_model`, `record_measurement`, `STATUS_CANDIDATE`, `STATUS_UNSUITABLE` | generation, `hpc label`, `hpc record` |
| [`plumetools.cai2012ext.generation`](generation.md) | stage, compare and classify every case; write the manifest | `classify`, `decide`, `apply`, `has_solver_output`, `build_entry`, `write_manifest`, `STATUS_NEW`, `STATUS_IDENTICAL`, `STATUS_STALE`, `STATUS_CONFLICT` | `generate_cases.py` |
| [`plumetools.cai2012ext.provenance`](provenance.md) | what every generated case was built from | `git_provenance`, `sources`, `independent_knudsen` | manifest |
| [`plumetools.cai2012ext.verification`](verification.md) | check a generated tree against the extension's promises | `verify_extension` | `validate_cases.py`, `smoke_cases.py` |
| [`plumetools.cai2012ext.hpc`](hpc.md) | what one PBS job does around `./Allrun` | `guard`, `label`, `reached_end`, `record`, `MARKER` | `runCase.sh` |
| [`plumetools.cai2012ext.diagnostics`](diagnostics.md) | the nozzle-plume rarefaction diagnostics, on a Cai case | `run`, `load_settings`, `centreline_profile` | `Alldiagnose.py` |
| [`plumetools.cai2012ext.comparison`](comparison.md) | off-axis and angular comparison with Cai's field | `angular_profile`, `angular_metrics`, `centerline_extended`, `meridional_plane`, `plane_metrics`, `masked_metrics` | `Allcompare.py` |
| [`plumetools.cai2012ext.convergence`](convergence.md) | far-field stationarity across the sampling window | `assess`, `regional_densities`, `linear_trend`, `DEFAULT_TOLERANCE` | `Allconverge.py` |
| [`plumetools.cai2012ext.figures`](figures.md) | figures on fixed normalised scales | `plot_angular`, `plot_plane`, `plot_stationarity`, `status_banner` | post |
| [`plumetools.cai2012ext.summary`](summary.md) | the study table and the Knudsen comparison | `study_row`, `write_table`, `write_markdown`, `plot_knudsen`, `COLUMNS` | `analyse.py` |
| [`plumetools.cai2012ext.vizspec`](vizspec.md) | per-Kn imagery specs from the manifest | `group_specs`, `write_group_specs` | `AllpostCases` |

## Interface track

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.cai2012ext.interface`](interface.md) | the pure model: equivalent source, Bird's parameter, sphere selection, inlet faces | `equivalent_source`, `bird_parameter`, `select_sphere`, `face_inflow`, `flux_balance`, `exit_targets`, `SphereSelection` |
| [`plumetools.cai2012ext.interface_config`](interface_config.md) | `interface.yaml`: read and validate; the one place the mode is decided | `InterfaceConfig`, `parse_interface`, `read_settings`, `InterfaceConfigError` |
| [`plumetools.cai2012ext.interface_case`](interface_case.md) | derive, mesh and inject one interface-mode case; the CLI | `derive`, `Derived`, `InterfacePlan`, `guard`, `label`, `write_mesh_dictionaries`, `write_inflow`, `verify_mesh`, `check_decomposed`, `main` |
| [`plumetools.cai2012ext.interface_assessment`](interface_assessment.md) | the pre-run verdict; writes nothing | `assess`, `STATUS_INFEASIBLE`, `STATUS_PENDING`, `STATUS_VALIDATION`, `STATUS_SMOKE` |
| [`plumetools.cai2012ext.interface_evidence`](interface_evidence.md) | what a production coupling verdict may rest on | `coupling_validation`, `compatibility`, `evidence_provenance`, `REQUIRED` |
| [`plumetools.cai2012ext.interface_identity`](interface_identity.md) | restart compatibility: what a mesh and its parcels were built for | `identity`, `write_identity`, `identity_problem`, `IDENTITY_FILE` |
| [`plumetools.cai2012ext.interface_mesh`](interface_mesh.md) | the truncated mesh: the Cai box with the core carved out | `plan`, `InterfaceMeshPlan`, `write_mesh_setup`, `render_snappy_hex_mesh_dict` |
| [`plumetools.cai2012ext.interface_report`](interface_report.md) | report file, manifest row, console lines, evidence collection | `report`, `manifest_entry`, `describe`, `write_report`, `collect` |
| [`plumetools.cai2012ext.interface_study`](interface_study.md) | `interfaces.yaml`, its cases and manifest | `load`, `InterfaceStudy`, `Variant`, `stage`, `classify`, `write_manifest` |
| [`plumetools.cai2012ext.overlap`](overlap.md) | a truncated run against the full-domain benchmark | `run`, `parse_overlap`, `OverlapConfig`, `compare`, `verdict`, `main` |

## Field scaling

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.cai2012ext.scaling`](scaling.md) | coordinate-only deformations of the unchanged Cai field | `CoordinateMap`, `CaiPredictor`, `kn_feature`, `FAMILIES` |
| [`plumetools.cai2012ext.scaling_fit`](scaling_fit.md) | separate whole-field and contour fits | `prepare_case`, `fit`, `evaluate`, `FitResult`, `sensitivity_report` |
| [`plumetools.cai2012ext.scaling_inventory`](scaling_inventory.md) | read-only inventory and opt-in registration | `inventory`, `read_registration`, `reference_suitability`, `prepare_missing_specs` |
| [`plumetools.cai2012ext.scaling_report`](scaling_report.md) | numerical exports and paired figures | `render_report`, `run_sensitivities` |
| [`plumetools.cai2012ext.scaling_workflow`](scaling_workflow.md) | study-scoped orchestration; never launches DSMC | `main`, `run_assessment`, `load_config`, `synthetic_cases` |

## Deriving and assessing a case

`direct_case.derive(case_dir)` returns a `DirectCase` with `cfg`, `geom`,
`exit_state`, `plan`, `run`, `nozzle`, `particles` and `report`; it still
unpacks in the old 8-tuple order for code written against
`generation.derive`. `feasibility.assess_case` reads it against the production
criteria and returns a verdict dictionary — the same one the manifest stores:

<!-- docs-example: cai2012ext-feasibility -->
```python
from pathlib import Path
from plumetools.cai2012ext import direct_case, feasibility

case = direct_case.derive(Path("cases/cai2012-kn-extension/baseCase"))   # the template, Kn = 10
verdict = feasibility.assess_case(case)
assert verdict["status"] in (feasibility.STATUS_CANDIDATE, feasibility.STATUS_UNSUITABLE)
# keys: status, production_ready, criteria, unresolved, resolution, time_step, particles, cost, note
```

`verdict["unresolved"]` lists the criteria a case misses; it is what `hpc
label` writes into a `NON-PRODUCTION` marker. Nothing here meshes or runs.

## Generation semantics

`generation.classify` stages each case and compares it with the disk;
`generation.decide(states, refresh=…, regenerate=…, overwrite=…)` sets the
action (`create`, `keep`, `keep-stale`, `refresh`, `replace`, `refuse`);
`generation.apply` carries it out and raises if anything is refused.
`has_solver_output` is wider than `cai2012.study.case_is_complete`: a time
directory in the case **or any processor directory**, or a `log.dsmcFoam`,
counts. The table of statuses and flags is in the
[tutorial](../../../tutorials/cai2012-kn-extension.md#generation-exactly).

# Legacy source flow (3d-inflow)

The top-level modules that implement the **legacy** plume source-flow model and
its case: the frozen equations E1–E8, face geometry, the gas table, the
`case.yaml` loader (also used by `markelov1999`), the orchestration that
evaluates the model over the inflow patch, and the mesh-geometry verifier. The
regression golden pins their output for `cases/3d-inflow` to the pre-refactor
script, defects included.

Depends on: `numpy`, `scipy`, `PyYAML`; `plumetools.mesh`, `plumetools.foamio`.
Workflow: [3d-inflow](../../../tutorials/3d-inflow.md). The model:
[The plume source-flow model](../../../source-flow-model.md).

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.config`](config.md) | `case.yaml` loading and validation (source-flow and Markelov schemas) | `load_case_config`, `CaseConfig`, `validate_against_case`, `ConfigError`, `ConfigWarning`, `MESH_TYPES`, `INFLOW_MODELS`, `DIALECTS` |
| [`plumetools.sourceflow`](sourceflow.md) | equations E1–E8, pure functions — the frozen legacy model | `limiting_angle`, `limiting_velocity`, `angular_separable`, `angular_theta_only`, `normalization_coefficient`, `number_density`, `velocity` |
| [`plumetools.geometry`](geometry.md) | centroids, normals, areas, spherical coordinates | `centroids`, `normals`, `areas`, `spherical` |
| [`plumetools.inflow`](inflow.md) | evaluate the model over a case's inflow patch | `compute_inflow`, `InflowResult`, `read_patch_geometry`, `area_weighted_number_density` |
| [`plumetools.species`](species.md) | the gas property table | `get`, `Species` |
| [`plumetools.constants`](constants.md) | constants and unit conversions | `BOLTZMANN_J_PER_K`, `AVOGADRO_PER_MOL`, `AMU_TO_KG`, `PSI_TO_PA` |
| [`plumetools.verify_mesh`](verify_mesh.md) | the built mesh against `case.yaml`; CLI | `verify`, `MeshVerificationError`, `main` |

## Evaluating the inflow

`inflow.compute_inflow(case_dir, cfg)` reads the inflow patch (the
`constant/polyMesh/sets/<patch>` face set if present, else the patch's face
range — the same order), evaluates the model at each face centroid, and returns
an `InflowResult` with `labels`, `rhoN` (m⁻³), `U` (m/s, shape `(n, 3)`) and
`T` (K). `foamio.fields.write_inflow_fields` writes it in the dialect
`output.dialect` names. Under the standard dialect,
`area_weighted_number_density` gives the single number density stock
`FreeStream` can take; the per-face density is then not used by the solver.

!!! warning "Conventions are frozen, not chosen"
    θ is measured from +z, not from the +x plume axis (SM-04); the model gas is
    N₂ while the solver simulates argon (SM-03); the pressure is the literal
    475 psi (SM-01). Each is a `LEGACY` value in `case.yaml`. The equations and
    the full defect list: [The plume source-flow model](../../../source-flow-model.md).

The family pages for every other package are under [API reference](../index.md).

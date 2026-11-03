# `plumetools.foamio`

OpenFOAM files, shared by every family: reading a field's values with their
metadata, numeric time directories, restart checkpoints across ranks, the
patch-field entry each patch type takes, and the dictionary rendering the
meshing families share. Also the legacy 3d-inflow `0/` writer and mesh
dictionaries, frozen by the regression golden.

Depends on: `numpy` (field values, block mesh); nothing else outside the
standard library. Used by every family.

::: plumetools.foamio
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - write_mesh_setup
        - mesh_pipeline

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.foamio.errors`](errors.md) | the one error the readers raise (a `ValueError`) | `FoamDataError` |
| [`plumetools.foamio.times`](times.md) | which subdirectories are time directories | `time_directories`, `time_names`, `time_value`, `latest_time` |
| [`plumetools.foamio.volfield`](volfield.md) | a field's `internalField` or one patch's `value`, with metadata | `read_internal_field`, `read_patch_field`, `read_internal`, `read_patch`, `FieldValues`, `parse_value` |
| [`plumetools.foamio.patchfields`](patchfields.md) | the patchField type each patch type takes | `patch_field_type`, `nonuniform_value_entry`, `CONSTRAINT_PATCH_TYPES` |
| [`plumetools.foamio.checkpoints`](checkpoints.md) | restart checkpoints, complete on every rank or not, and what a resume must do; CLI | `scan`, `Checkpoints`, `plan_restart`, `RestartPlan`, `reached_end`, `main` |
| [`plumetools.foamio.primitives`](primitives.md) | searchable sphere/cylinder/cone/box and dictionary helpers | `SearchableSphere`, `SearchableCylinder`, `SearchableCone`, `SearchableBox`, `render_snappy_dict`, `render_box_block_mesh_dict`, `render_mesh_quality_dict` |
| [`plumetools.foamio.blockmesh`](blockmesh.md) | `blockMeshDict` for the 3d-inflow O-grid | `render_block_mesh_dict`, `write_block_mesh_dict`, `build_vertices` |
| [`plumetools.foamio.snappy`](snappy.md) | snappyHexMesh setup for the box-with-hemispherical-cavity | `write_snappy_setup`, `background_divisions`, `surface_cell_size`, `location_in_mesh` |
| [`plumetools.foamio.fields`](fields.md) | the legacy `0/` inflow and measurement fields | `write_inflow_fields`, `write_measurement_fields`, `MEASUREMENT_FIELDS` |

The package itself exports `write_mesh_setup` and `mesh_pipeline`, which pick
the generator and the OpenFOAM command list from a legacy `case.yaml`'s
`mesh.type`.

## Rules the readers keep

- A list must hold the count it declares and, through `FieldValues.over(n, what)`,
  one value per cell or face; only a `uniform` entry is expanded. A binary
  file, a non-number, or a list from another mesh raises `FoamDataError` saying
  which. (Before 2026-10-27 a list declaring four values and holding one was
  silently broadcast.)
- Readers stream and close; nothing is cached between calls, so a remeshed case
  is never read through a stale copy.
- Time directory names are kept as written; comparison is numeric, identity is
  by name.

Arrays: `read_internal_field` returns `(n,)` for a scalar field and `(n, k)`
for a vector or tensor; `read_patch_field(path, patch)` the same for one patch.

## Checkpoints

A time is a **complete checkpoint** on a rank when `T/uniform/time` exists,
every field file in `T/` and in `T/lagrangian/<cloud>/` ends with OpenFOAM's
end divider, `T/` holds every field file any other rank holds at `T` (and the
same for the cloud), and a rank whose previous time carried a cloud carries one
at `T`. A decomposed checkpoint is complete on every `processor*`. Only the
newest times are inspected; nothing is deleted, moved or edited.

`plan_restart(case_dir, n_subdomains)` returns a `RestartPlan`: `resume` (the
time to continue from, or `None`), `case_time`, `proc_time`, `proc_count`,
`reached_end`, `decompose` (`none` — reuse processor directories; `all` — fresh
start; `latest` — decompose only the resume time), `reconstruct_first` (the
processors hold newer complete state that `decomposePar -force` would delete),
and `refusal` (why the case must not be continued, or `None`).

<!-- docs-example: foamio-checkpoints -->
```python
import tempfile
from pathlib import Path
from plumetools.foamio import checkpoints

with tempfile.TemporaryDirectory() as tmp:
    case = Path(tmp)
    (case / "system").mkdir()
    (case / "system" / "controlDict").write_text("endTime 1;\ndeltaT 0.1;\n")
    plan = checkpoints.plan_restart(case, 4)          # nothing written yet: a fresh start
    assert plan.resume is None and plan.decompose == "all" and plan.refusal is None
```

How `Allrun` uses the plan, and how to recover from a refusal:
[Restart and recovery](../../../guides/restart-and-recovery.md).

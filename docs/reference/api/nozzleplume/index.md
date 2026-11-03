# `plumetools.nozzleplume`

Nozzle-exit plumes into vacuum, with or without a plate. It is the Cai problem
with the paper taken out: a uniform exit state on a circular exit, given
**directly** (density, speed, temperature — not derived from a Knudsen number),
any single VHS gas, an optional flat plate, and a snappyHexMesh mesh. It
re-implements its gas and exit state rather than importing `cai2012`; a unit
test holds the two to the same exit state for Cai's inputs.

Depends on: `numpy`, `PyYAML`, `matplotlib` (figures); `plumetools.foamio`,
`plumetools.mesh`, `plumetools.rarefaction`. Workflow:
[r-4d](../../../tutorials/r-4d.md).

::: plumetools.nozzleplume
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - MODEL

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.nozzleplume.config`](config.md) | `case.yaml` schema, loading, validation | `load_case_config`, `CaseConfig`, `ExitConfig`, `PlateConfig`, `RegionConfig`, `ConfigError` |
| [`plumetools.nozzleplume.constants`](constants.md) | constants | `BOLTZMANN_J_PER_K`, `STANDARD_GRAVITY_M_PER_S2` |
| [`plumetools.nozzleplume.gas`](gas.md) | one VHS species and the kinetic relations | `VhsSpecies`, `mean_free_path`, `collision_time`, `maxwellian_number_flux`, `maxwellian_momentum_flux`, `maxwellian_energy_flux` |
| [`plumetools.nozzleplume.exitstate`](exitstate.md) | the exit state and everything it implies (Kn, Mach, mass flow, thrust, Isp) | `ExitState`, `from_config` |
| [`plumetools.nozzleplume.geometry`](geometry.md) | exit disk, vacuum box, optional plate | `PlumeGeometry`, `Plate`, `from_config` |
| [`plumetools.nozzleplume.mesh`](mesh.md) | cubic background, refinement cones and boxes, plate carving, cell estimate | `plan`, `MeshPlan`, `regions_from_config`, `estimate_cells`, `exit_level`, `write_mesh_setup` |
| [`plumetools.nozzleplume.run`](run.md) | time step, particle weight, run length | `derive`, `RunSettings` |
| [`plumetools.nozzleplume.checks`](checks.md) | DSMC mesh, time-step and boundary checks | `run_checks`, `CheckReport`, `expected_patch_types` |
| [`plumetools.nozzleplume.dictionaries`](dictionaries.md) | `dsmcProperties`, `controlDict`, … | `write_case_dictionaries`, `render_control_dict`, `simple_split` |
| [`plumetools.nozzleplume.foamfields`](foamfields.md) | the `0/` fields | `write_all`, `write_inflow_fields` |
| [`plumetools.nozzleplume.verify`](verify.md) | the built mesh against `case.yaml` | `verify`, `MeshVerificationError` |
| [`plumetools.nozzleplume.post`](post.md) | plate loads and the centreline | `plate_loads`, `PlateLoads`, `centreline`, `write_plate_results`, `write_centreline_results`, `plot_results` |
| [`plumetools.nozzleplume.diagnostics`](diagnostics.md) | rarefaction and mesh-resolution diagnostics for a case | `run`, `reference_lengths`, `centreline_profile` |
| [`plumetools.nozzleplume.workflow`](workflow.md) | the steps `Allmesh`, `Allrun`, `Allpost` call, and the CLI | `load`, `Case`, `describe`, `setup`, `verify`, `fields`, `post`, `diagnostics`, `main` |
| [`plumetools.nozzleplume.meshread`](meshread.md) | compatibility names; geometry now lives in `plumetools.mesh` | (re-exports only) |
| [`plumetools.nozzleplume.__main__`](__main__.md) | `python -m plumetools.nozzleplume` | (calls `workflow.main`) |

## Loading a case

`workflow.load(case_dir)` reads `case.yaml` and derives everything that follows
from it into a typed `Case` with `case_dir`, `cfg`, `geom`, `exit_state`,
`plan`, `run`, `nozzle`, `checks` and `label`. Every step consumes that object.

<!-- docs-example: nozzleplume-load -->
```python
from pathlib import Path
from plumetools.nozzleplume import workflow

case = workflow.load(Path("cases/r-4d/baseCase"))
assert case.cfg.plate is not None                    # baseCaseFreePlumeHPC has plate: null
mfp = case.exit_state.mean_free_path_m               # m, from the exit density and the gas
```

`workflow.describe(case_dir)` prints the same derivation and writes nothing —
it is `./Allmesh --describe`. Mesh planning rules (cones written as
`searchableCone`, the exit at one refinement level, the staircase nozzle) are
explained in [r-4d, explained](../../../r-4d-case.md#the-mesh).

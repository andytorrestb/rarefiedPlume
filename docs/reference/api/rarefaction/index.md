# `plumetools.rarefaction`

After a run: the local mean free path, Knudsen numbers and mesh resolution as
volume fields written back into the sampled time directory, plus a summary for
planning the next mesh. It works on any reconstructed `dsmcFoam` case; the
nozzle-plume family and the Kn extension call it with settings from their own
configuration, and `python -m plumetools.rarefaction` is the door for every
other case.

Depends on: `numpy`, `PyYAML`; `plumetools.foamio`, `plumetools.mesh`,
`plumetools.nozzleplume.gas` (the VHS species). Workflows:
[r-4d](../../../tutorials/r-4d.md), [Kn extension](../../../tutorials/cai2012-kn-extension.md).
The fields, their thresholds and the run → inspect → remesh workflow:
[Rarefaction diagnostics](../../../rarefaction-diagnostics.md).

::: plumetools.rarefaction
    options:
      show_root_heading: false
      members: false

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.rarefaction.settings`](settings.md) | thresholds, lengths, field names | `DiagnosticSettings`, `SettingsError`, `DEFAULT_FIELD_NAMES`, `PARTICLE_COUNT_BASES` |
| [`plumetools.rarefaction.core`](core.md) | the numerics and the diagnostic registry | `mean_free_path`, `translational_temperature`, `least_squares_gradient`, `cell_status`, `register`, `evaluate`, `knudsen_on`, `CellGraph` |
| [`plumetools.rarefaction.foam`](foam.md) | read a sampled time; write fields back into it, streamed | `load`, `CaseData`, `read_window`, `SamplingWindow`, `write_fields`, `weight_check` |
| [`plumetools.rarefaction.run`](run.md) | one case, one sampled time: read, compute, write, summarise | `diagnose`, `Diagnosis` |
| [`plumetools.rarefaction.summary`](summary.md) | condense the fields into what a remeshing decision needs | `build`, `report`, `write_yaml`, `mesh_estimate`, `worst_cells`, `required_level` |
| [`plumetools.rarefaction.__main__`](__main__.md) | `python -m plumetools.rarefaction CASE --settings FILE` | `load_settings_file`, `main` |

## Mean free path, element by element

<!-- docs-example: rarefaction-mfp -->
```python
import math
from plumetools.nozzleplume.gas import VhsSpecies
from plumetools.rarefaction.core import mean_free_path

gas = VhsSpecies(name="gas", mass_kg=3.40055e-26, diameter_m=3.5e-10, omega=0.5, t_ref_K=273.0)
lam = mean_free_path([1.306e22, 0.0], [498.04, 300.0], gas)     # m; a hard sphere, so T drops out
assert abs(lam[0] - 1.0 / (math.sqrt(2.0) * math.pi * 3.5e-10**2 * 1.306e22)) < 1e-15
assert math.isnan(lam[1])                                        # n = 0 is "no data", not infinity
```

`mean_free_path(number_density, temperature, species, convention="vhs")`
applies (T/T_ref)^(ω−½) to the hard-sphere value for the `vhs` convention and
nothing for `hard_sphere`; with ω = ½ the two agree. It returns NaN where the
input is undefined, because here a missing value is a cell with no data that a
statistic must not count. The temperature is the translational temperature of
the averaged moments (`translational_temperature`), and a temperature at the
round-off of its moments is treated as no temperature (`TEMPERATURE_RESOLUTION`).

## Adding a diagnostic

A diagnostic is one function registered in `core` with `register`; a new input
is read in `foam` with the shared readers. The procedure and the conventions it
must keep are in [Rarefaction diagnostics § Adding a diagnostic](../../../rarefaction-diagnostics.md#adding-a-diagnostic).

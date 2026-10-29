# `plumetools.hybridplume`

A steady one-way hybrid CFD → DSMC free plume in full 3-D, for Cai & Wang's
circular argon source: a `rhoCentralFoam` precursor, continuum-breakdown
indicators on its cells, a fixed coupling surface chosen and assessed against
explicit criteria, the CFD state mapped onto the DSMC inflow faces, and DSMC
measurements, budgets and an independent collisionless reference. Each step is
a function in `workflow` and a CLI command; OpenFOAM tools are run by
`cases/cai2012-hybrid/Allrun`, never from Python.

Depends on: `numpy`, `scipy`, `PyYAML`, `matplotlib` (report figures);
`plumetools.foamio`, `plumetools.mesh`. Workflow:
[cai2012-hybrid](../../../tutorials/cai2012-hybrid.md). Method and verification
record: [`plumetools.hybridplume`, explained](../../../hybridplume.md).

::: plumetools.hybridplume
    options:
      show_root_heading: false
      members: false

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.hybridplume.config`](config.md) | `hybrid.yaml`: load, overlay a profile, validate | `load`, `CaseConfig`, `validate`, `deep_merge`, `build_source`, `departures_from_paper`, `ConfigError` |
| [`plumetools.hybridplume.physics`](physics.md) | the single source of every derived physical quantity | `Gas`, `SourceState`, `BIRD_ARGON`, `MEAN_FREE_PATH_CONVENTIONS` |
| [`plumetools.hybridplume.mesh`](mesh.md) | the shared 3-D Cartesian mesh and the staircase source | `plan`, `MeshPlan`, `predicted_source`, `write_mesh_setup` |
| [`plumetools.hybridplume.meshgeom`](meshgeom.md) | whole-mesh geometry: faces, cells, patches, adjacency | `read_mesh`, `MeshGeometry`, `cell_geometry` |
| [`plumetools.hybridplume.cfd`](cfd.md) | the `rhoCentralFoam` precursor case, with tabulated transport | `write_case`, `transport_table`, `render_thermophysical_properties`, `background_state` |
| [`plumetools.hybridplume.breakdown`](breakdown.md) | Kn_ρ, Kn_T, Kn_U, Kn_GLL and B_exp on the CFD cells | `read_state`, `indicators`, `continuum_mask`, `threshold_sensitivity`, `INDICATORS` |
| [`plumetools.hybridplume.interface`](interface.md) | select the coupling surface and assess it | `extract`, `assess`, `manual_region`, `overlap_layers`, `Surface`, `STATUS_SUITABLE`, `STATUS_UNSUITABLE` |
| [`plumetools.hybridplume.mapping`](mapping.md) | CFD interface state → DSMC inflow faces, with the flux error | `match_interface`, `inflow_from_fields`, `flux_change`, `consistency` |
| [`plumetools.hybridplume.halfrange`](halfrange.md) | exact half-range fluxes of a drifting Maxwellian | `normal_moments`, `incoming_flux`, `net_flux`, `bird_eq_4_22_number_flux` |
| [`plumetools.hybridplume.fluxes`](fluxes.md) | Navier–Stokes–Fourier flux densities through faces | `continuum_flux`, `viscous_stress` |
| [`plumetools.hybridplume.dsmc`](dsmc.md) | the standard-`dsmcFoam` case: run plan, dictionaries, `0/` | `plan_run`, `RunPlan`, `write_case`, `write_keep_set`, `write_fields` |
| [`plumetools.hybridplume.particles`](particles.md) | moments and boundary crossings from the written parcels | `read_cloud`, `MomentAccumulator`, `classify_exits`, `crossing_window` |
| [`plumetools.hybridplume.budgets`](budgets.md) | control-volume budgets: precursor, DSMC region, interface | `cfd_budgets`, `dsmc_budget`, `reconstructed_patch_fluxes` |
| [`plumetools.hybridplume.collisionless`](collisionless.md) | an independent collisionless reference for the source | `reference`, `raw_moments`, `macroscopic`, `quadrature` |
| [`plumetools.hybridplume.post`](post.md) | DSMC moments, profiles, azimuthal statistics, comparisons | `dsmc_fieldaverage`, `line_profile`, `azimuthal_statistics`, `compare` |
| [`plumetools.hybridplume.resources`](resources.md) | resource and resolution assessment before launch | `assess_configuration`, `cfd_estimate`, `dsmc_estimate`, `measured_rates`, `rank_imbalance` |
| [`plumetools.hybridplume.provenance`](provenance.md) | repository SHAs, source hashes, solver binaries | `collect`, `git_state`, `source_hashes`, `toolchain` |
| [`plumetools.hybridplume.hpc`](hpc.md) | resumable batch scripts sized from the estimates | `write_scripts`, `walltime_hours`, `PBS_DEFAULTS` |
| [`plumetools.hybridplume.report`](report.md) | every step's JSON into one summary, and the figures | `summarise`, `plots` |
| [`plumetools.hybridplume.workflow`](workflow.md) | the explicit steps | `setup`, `check_mesh`, `cfd_report`, `breakdown`, `interface`, `dsmc_setup`, `dsmc_fields`, `dsmc_post`, `crossings`, `reference`, `overlap`, `hpc`, `report`, `WorkflowError` |
| [`plumetools.hybridplume.__main__`](__main__.md) | `python -m plumetools.hybridplume <step> <case>` | `main` |

## Configuration and the half-range moments

<!-- docs-example: hybrid-config -->
```python
import math
from pathlib import Path
from plumetools.hybridplume import config, halfrange

cfg = config.load(Path("cases/cai2012-hybrid/hybrid.yaml"), "smoke-kn0p01-coarse")
assert cfg.get("source.knudsen") == 0.01 and cfg.get("mesh.core_cell_over_D") == 0.1

# J1, J2, J3: half-range moments of the normal velocity over v.n > 0, per unit
# density, for a normal drift (m/s) and a thermal speed a = sqrt(2RT) (m/s)
a = 300.0
J1, J2, J3 = halfrange.normal_moments([0.0, 500.0], a)
assert abs(J1[0] - a / (2.0 * math.sqrt(math.pi))) < 1e-9     # zero drift: the one-sided flux
```

The half-range closed forms are in the module's docstring; with zero drift J1
is the one-sided Maxwellian number flux per unit density, a/(2√π). These are what make the injected flux
at the interface exact rather than sampled, and what
`tests/unit/test_hybridplume_physics.py` and `test_hybridplume_interface.py`
check against quadrature.

## The interface verdict

`interface.assess` evaluates a candidate surface against `interface.criteria`
in the configuration — returning mass fraction, area fraction with a subsonic
normal Mach number, contact with the outer boundary, number of overlap layers
— and records `STATUS_SUITABLE` or `STATUS_UNSUITABLE` with the measured values
(`STATUS_FULL_DSMC` when `interface.mode` is `none`). `workflow.dsmc_setup`
refuses an unsuitable surface unless `interface` ran with
`allow_unsuitable=True`, in which case every downstream report is labelled a
machinery smoke test.

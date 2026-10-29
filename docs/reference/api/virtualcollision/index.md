# `plumetools.virtualcollision`

Reduced virtual pair collisions for a single-species uniform disk plume: a
uniform disk flux source with its geometric (non-Maxwellian) ballistic
distribution, a verified single-species elastic VHS kernel matching OpenFOAM
v2512, mutually thinned representative pair events in a finite plume, raw
moment algebra, and the read-only assessment against Cai and DSMC.

Depends on: `numpy`, `scipy`, `PyYAML`, `matplotlib`. Workflow:
[virtual-collision assessment](../../../tutorials/cai2012-virtual-collision.md).
Mathematical contract: [Reduced virtual collisions](../../../cai2012-virtual-collision.md).

::: plumetools.virtualcollision
    options:
      show_root_heading: false
      members: false

## Modules

| Module | Purpose | Principal interfaces |
|---|---|---|
| [`plumetools.virtualcollision.source`](source.md) | uniform disk flux source and ballistic f₀ | `DiskSource` |
| [`plumetools.virtualcollision.kernel`](kernel.md) | the single-species elastic VHS kernel of OpenFOAM v2512 | `VHSKernel`, `ZeroKernel`, `elastic_pair`, `UnsupportedCollisionLaw` |
| [`plumetools.virtualcollision.model`](model.md) | reduced, mutually thinned representative pair events | `VirtualCollisionModel`, `NumericalOptions`, `ray_exit`, `MODEL_VERSION` |
| [`plumetools.virtualcollision.moments`](moments.md) | raw moment algebra; no clipping or covariance repair | `derive`, `aggregate` |
| [`plumetools.virtualcollision.demonstration`](demonstration.md) | bounded predictive-only convergence, conservation and cost run | `run_model_demonstration` |
| [`plumetools.virtualcollision.assessment`](assessment.md) | read-only inventory and DSMC-referenced metrics | `inventory`, `run_reference_assessment`, `dsmc_reference_metrics`, `time_block_standard_error`, `contour_location_discrepancy` |
| [`plumetools.virtualcollision.export_geometry`](export_geometry.md) | recover comparison geometry from verified Cartesian mesh inputs | `load_export_context`, `point_volumes`, `meridional_geometry`, `UnsupportedExportGeometry` |
| [`plumetools.virtualcollision.exported`](exported.md) | three-way comparisons from compact Cai-family exports | `run_case`, `run_exports`, `PointPrediction`, `BallisticCache` |
| [`plumetools.virtualcollision.export_figures`](export_figures.md) | three-method figures from preserved CSVs and predictions | `plot_case`, `plot_density_contours`, `DENSITY_CONTOUR_LEVELS` |
| [`plumetools.virtualcollision.export_summary`](export_summary.md) | index and four-Kn views of saved runs | `summarize`, `read_table` |

## What is verified and what is not

`kernel` reproduces the OpenFOAM v2512 VHS collision law for one species, with
its Boltzmann constant (`OPENFOAM_V2512_BOLTZMANN_J_PER_K`); any other law
raises `UnsupportedCollisionLaw` rather than being approximated. The unit tests
(`tests/unit/test_virtualcollision_*.py`) check the kernel, the model's
conservation and sampling properties, the export-geometry recovery and the
figures. The comparison of model and DSMC is **exploratory**: the cell-centre
quadrature and the continuum ring estimator have a reported spatial
discrepancy, so agreement or disagreement is not yet a statement about the
model's accuracy.

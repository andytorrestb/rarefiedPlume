# Virtual-collision numerical assessment against DSMC

This separate read-only workflow compares the existing Cai evaluator, DSMC,
and the **Virtual-collision model**. It never generates cases, starts DSMC,
submits a scheduler job, or changes the existing family policy gates.
The derivation, verified VHS kernel provenance, closure, and validity limits are
in [the scientific-method document](../../docs/cai2012-virtual-collision.md).
Actual timings, convergence sensitivities, tests and local-data limitations are
in [the bounded run report](run-report.md).

From the repository root, with the existing Python environment and plotting
dependencies installed:

```bash
python cases/cai2012-virtual-collision/assess.py --help
python cases/cai2012-virtual-collision/assess.py inventory
python cases/cai2012-virtual-collision/assess.py smoke --output cases/cai2012-virtual-collision/results/smoke
python cases/cai2012-virtual-collision/assess.py full --output cases/cai2012-virtual-collision/results-full
python cases/cai2012-virtual-collision/compare_exports.py --samples 2048 --grid-refinement 1
# Select multiple cases; no DSMC is launched:
python cases/cai2012-virtual-collision/compare_exports.py --kn 10 1 --samples 2048
# Refresh the overview from saved predictions without evaluating the model:
python cases/cai2012-virtual-collision/summarize_exports.py
```

`smoke` is bounded by 64 flux samples, four independent randomized batches,
eight hazard steps, and a 4×8×12 collision-partner quadrature; convergence runs
change one numerical control at a time. All thread counts are capped at eight
(one by default). `full` uses the seven requested Kn values and the larger
numerical settings in `assessment.yaml`; it is intended for HPC and still
performs postprocessing only. No physical parameter is fitted to reference data.

For a manual PBS run, follow the existing environment convention:

```bash
cd cases/cai2012-virtual-collision
qsub runAssess.sh
# Another site's environment:
qsub -v VC_HPC_ENV=/path/to/your/hpc_env.sh runAssess.sh
```

The template defaults to the existing extension family's `hpc_env.sh`, which
contains that site's paths. The new template has no hardcoded user path. The
batch process needs Python and the package, with no LLM service.

`inventory.json` and `inventory.csv` discover actual case folders, including
ignored `Cases*` trees and interface cases. They retain configured and inverse
checked Kn, the nozzle-diameter and mean-free-path conventions, actual collision
dictionaries where available, exit/source/geometry parameters, available times
and exports, resolution/occupancy/stationarity evidence, runtime provenance,
quality labels, and explicit exclusion reasons. Kn 0.01 is the additional
Cai-family comparison omitted from the original six-point list. Interface
sources are excluded from this nozzle-source correction. Kn 0.001 and 0.0001
direct-exit cases remain non-production regardless of apparent stationarity.

`report.json`, `cost.csv`, `convergence.csv`, model raw arrays, and the model-only
plot describe the bounded numerical demonstration. Preprocessing and evaluation
costs are separate. A demonstration without a compatible reference establishes
neither accuracy nor a speed advantage over DSMC production. Unavailable DSMC
resource accounting and confidence intervals remain unavailable.

When full DSMC averaged fields and physical `C`/`V` geometry exist, the workflow
also writes per-case raw arrays, DSMC-reference metrics, same-level three-way
density contours, signed error maps, central-bin profiles, angular distributions,
contour-location distances, and an error-versus-Kn summary. Every comparison uses
the same intended cylindrical bins and common masks. Density errors use physical
bin-volume weights and the DSMC global norm; velocities and temperatures use
stated physical scales. Raw first and second moments are aggregated before
forming central moments; scalar DSMC temperature does not establish tensor
agreement. The cell-centre quadrature and continuum ring estimator have a
reported spatial discrepancy, so these comparisons remain exploratory until
spatial convergence has been established. Raw arrays retain invalid values;
display masks and metric exclusions never alter predictions.

The four local result sets contain compact averaged exports rather than original
`C/V` and full sampled fields. A later audit established that their tensor-product
Cartesian geometry **can be recovered** from `case.yaml` and the archived
`system/blockMeshDict`. `compare_exports.py` checks vertices, block cell counts,
grading, exported point coordinates, tube counts and angular footprints before
using this reconstruction. It writes three-method plots, augmented CSVs, model
raw moments/tensors, numerical uncertainty and diagnostics to each case's
`results/virtual_collision/`. Original DSMC/Cai exports are preserved. The study
index is `cases/cai2012-kn-extension/results/virtual_collision/README.md`.
The added `density_contour_comparison.png` (and SVG) overlays all three methods
at the legacy levels `n/n0 = 0.1, 0.01, 0.001`: DSMC red solid, Cai black dashed,
and the virtual-collision model orange dash-dot. Missing contour crossings and
common display-mask counts are recorded in its JSON file. The study index also
includes `density_contours_all_kn_three_way.png`. Refresh these views from saved
predictions with `summarize_exports.py`; no model evaluation is needed.

The collision correction uses a finite cylindrical-ring grid, interpolating
signed raw moments before adding the exact Cai point field. This smoothing is
an additional numerical approximation; batch arrays propagate sampling
covariance through interpolation and profile aggregation. Profiles retain the
legacy molecule-weighted local-cell temperature estimator, with a separately
named raw-aggregate temperature. Original DSMC tensor components and independent
instantaneous time samples remain unavailable. The graded Kn 1/10 stationarity
is inconclusive, and Kn 0.001/0.0001 references remain non-production. These are
exploratory comparisons, not established accuracy claims. New uniform Kn 1/10
folders without result CSVs are explicitly skipped; archived graded results are
processed separately.

The [results review and refinement plan](refinement-plan.md) records the mixed
accuracy results and numerical limitations. For the next implementation session,
use the [targeted sampling implementation prompt](sampling-implementation-prompt.md).
It specifies event-measure consistency, direction averaging, location
stratification and uncertainty fixes while preserving the existing physical
closure and saved predictions.

Focused checks:

```bash
python -m pytest tests/unit/test_virtualcollision_assessment.py -q
python -m pytest tests/unit/test_virtualcollision_kernel.py tests/unit/test_virtualcollision_model.py -q
```

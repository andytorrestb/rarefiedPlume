# Virtual-collision assessment: Cai, DSMC and a reduced collision model

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | Compare three predictions of the Cai free plume — Cai's collisionless evaluator, DSMC, and a **reduced virtual-collision model** that adds mutually thinned representative pair events to the collisionless field |
    | **Model** | `plumetools.virtualcollision`: uniform disk flux source, a verified single-species OpenFOAM v2512 VHS kernel, a finite-plume pair-event model, raw moment algebra |
    | **Evidence** | Unit-tested kernel, model, assessment and exports; a committed **bounded smoke** result (`results/smoke/`). Comparisons with DSMC are **exploratory**: spatial convergence of the comparison is not established |
    | **Runs on** | Python only (numpy, scipy, PyYAML; matplotlib for figures). Never generates cases, runs DSMC, or changes another family's policy gates |
    | **Start in** | the repository root (the scripts locate it themselves) |

The derivation, the kernel's provenance, the closure and its validity limits
are in [Reduced virtual collisions](../cai2012-virtual-collision.md). The
case directory also holds the working records of this effort:
[`run-report.md`](../../cases/cai2012-virtual-collision/run-report.md) (bounded
run, timings, sensitivities),
[`refinement-plan.md`](../../cases/cai2012-virtual-collision/refinement-plan.md)
(mixed accuracy results and numerical limits),
[`sampling-refinement-design.md`](../../cases/cai2012-virtual-collision/sampling-refinement-design.md)
and [`sampling-implementation-prompt.md`](../../cases/cai2012-virtual-collision/sampling-implementation-prompt.md)
(the next planned changes).

## Inputs

| File | Role |
|---|---|
| `assessment.yaml` | `families` to inventory (`cai2012`, `cai2012-kn-extension`), the Kn list (`100 … 0.0001`), `source_case` (the extension's `baseCase`, which supplies D, S₀, T₀ and the gas), `smoke` and `full` numerical settings (flux samples, hazard steps, partner quadrature, batches, seed, domain, chunk), `assessment_bins` |
| DSMC references | discovered, never generated: case folders of the two families, including ignored `Cases*` trees and interface cases |

No physical parameter is fitted to reference data. Interface sources are
excluded from this nozzle-source correction, and the Kn = 0.001 and 0.0001
direct-exit cases stay non-production whatever their apparent stationarity.

## Modes of `assess.py`

```bash
python cases/cai2012-virtual-collision/assess.py {inventory|smoke|full} \
    [--root REPO] [--config assessment.yaml] [--output DIR] [--no-plots]
```

| Mode | Does | Writes into `--output` |
|---|---|---|
| `inventory` | discovers every case folder and records configured and back-checked Kn, conventions, collision dictionaries, times and exports, resolution/occupancy/stationarity evidence, quality labels and exclusion reasons | `inventory.json`, `inventory.csv`, `skipped.json` |
| `smoke` | the above, then a bounded model demonstration (Kn = 10; 64 flux samples, 4 batches, 8 hazard steps, a 4×8×12 partner quadrature; one control changed at a time in the convergence runs) and the reference assessment where a compatible DSMC field exists | plus `report.json`, `cost.csv`, `convergence.csv`, model arrays (`*.npz`), `*_diagnostics.json`, the assessment CSV/JSON files, plots |
| `full` | the same with the seven Kn and the larger `full` settings — intended for a batch node | the same set |

!!! warning "Choose `--output` deliberately"
    The default is `cases/cai2012-virtual-collision/results/`, and
    `results/smoke/` holds the **committed** smoke result. Running
    `smoke --output cases/cai2012-virtual-collision/results/smoke` (the case
    README's example) rewrites those tracked files. For a check, write to a
    scratch directory instead, e.g. `--output /tmp/vc-smoke`.

Thread counts are capped at eight and default to one (`OMP_NUM_THREADS` and
friends are clamped before numpy loads).

## Cheap checks

```bash
python cases/cai2012-virtual-collision/assess.py inventory --output /tmp/vc-inventory
python cases/cai2012-virtual-collision/assess.py smoke --output /tmp/vc-smoke --no-plots
python -m pytest tests/unit/test_virtualcollision_kernel.py tests/unit/test_virtualcollision_model.py -q -p no:cacheprovider
```

The inventory is meaningful only where DSMC results exist; in a fresh clone it
records every case as missing and excluded, which is the correct answer. The
smoke run demonstrates the numerics. **A demonstration without a compatible
reference establishes neither accuracy nor a speed advantage over DSMC.**

## Comparing against saved DSMC exports

Four local result sets of the Kn extension were saved as compact averaged
exports rather than full fields with `C` and `V`. Their tensor-product
Cartesian geometry can be recovered from `case.yaml` and the archived
`system/blockMeshDict`; `compare_exports.py` checks vertices, block counts,
grading, exported coordinates, tube counts and angular footprints before using
the reconstruction.

```bash
python cases/cai2012-virtual-collision/compare_exports.py --samples 2048 --grid-refinement 1
python cases/cai2012-virtual-collision/compare_exports.py --kn 10 1 --samples 2048
python cases/cai2012-virtual-collision/summarize_exports.py      # redraw from saved predictions, no model evaluation
```

| Script | Options | Writes |
|---|---|---|
| `compare_exports.py` | `--study` (default `cases/cai2012-kn-extension`), `--samples` (2048), `--hazard-steps` (24), `--grid-refinement {0,1,2}` (1), `--kn …`, `--no-plots` | per case: `<case>/results/virtual_collision/` (three-method plots, augmented CSVs, model moments/tensors, uncertainty, diagnostics, `density_contour_comparison.{png,svg}`); study index `cases/cai2012-kn-extension/results/virtual_collision/README.md` and `density_contours_all_kn_three_way.png` |
| `summarize_exports.py` | — | refreshes the overview from saved predictions |
| `audit_hazard.py` | `--samples` (128), `--seed` (2026), `--steps` (24 48 96), `--kn` (10 1), `--output` | a hazard-sampling audit, by default under `cases/cai2012-kn-extension/results/virtual_collision/review/` |
| `review_results.py` | — | reviews saved comparisons without rerunning transport: `metrics_review`, `regional_review` and `validity_review` (`.json`/`.csv`) under `cases/cai2012-kn-extension/results/virtual_collision/review/` |

These write **into the Kn-extension family's result directories** (gitignored).
Original DSMC and Cai exports are preserved. In a fresh clone the four local
result sets do not exist and the cases are skipped; new uniform Kn 1/10 folders
without result CSVs are explicitly skipped.

## Batch run

```bash
cd cases/cai2012-virtual-collision
qsub runAssess.sh                                       # assess.py full --output results-hpc
qsub -v VC_HPC_ENV=/path/to/your/hpc_env.sh runAssess.sh
```

The job requests one node, 8 cores, 24 GB and 2 h on `shortq` (the owner's
settings) and by default sources the extension family's `hpc_env.sh` — which
holds the owner's site paths; pass `VC_HPC_ENV` for any other site. It is
post-processing only. `results-full/` and `results-hpc/` are gitignored.

## Outputs and how to read them

`report.json`, `cost.csv` and `convergence.csv` describe the bounded numerical
demonstration; preprocessing and evaluation costs are reported separately.
Where full DSMC averaged fields with physical `C`/`V` exist, the assessment
also writes per-case raw arrays, DSMC-referenced metrics, same-level three-way
density contours, signed error maps, central-bin profiles, angular
distributions, contour-location distances and an error-versus-Kn summary.

Conventions that matter when reading them: every comparison uses the same
cylindrical bins and common masks; density errors use physical bin-volume
weights and the DSMC global norm; velocities and temperatures use stated
physical scales; raw first and second moments are aggregated **before** central
moments are formed; a scalar DSMC temperature does not establish tensor
agreement; raw arrays keep invalid values, and display masks never alter
predictions. The cell-centre quadrature and the continuum ring estimator have a
reported spatial discrepancy, so these comparisons remain exploratory until
spatial convergence is shown. Unavailable DSMC resource accounting and
confidence intervals stay unavailable rather than estimated.

## Resources

`smoke` is bounded by the settings above; the committed `results/smoke/` is its
output. `full` is configured for a batch node (`runAssess.sh` above). Measured
timings for the bounded run are in `run-report.md`; none are recorded here.

## Figure placeholder

!!! figure "Figure placeholder — three-way density contours"
    **Subject:** `density_contour_comparison.png` — DSMC (red solid), Cai (black dashed) and the virtual-collision model (orange dash-dot) at n/n₀ = 0.1, 0.01, 0.001.
    **Source run:** `compare_exports.py --kn 10 1` on the saved Kn-extension exports (local archive, not in the repository).
    **Plotting command:** `python cases/cai2012-virtual-collision/compare_exports.py --kn 10 1 --samples 2048` (verified in source; writes the PNG per case).
    **Destination:** `docs/assets/figures/virtual-collision-contours.png`
    **Suggested caption:** "Collisionless, DSMC and reduced virtual-collision density contours (exploratory; spatial convergence not established)."

## Related

- [Reduced virtual collisions](../cai2012-virtual-collision.md) — the mathematical contract and the assessment protocol
- [Kn extension](cai2012-kn-extension.md) — where the DSMC references come from
- API: [`virtualcollision.model`](../reference/api/virtualcollision/model.md), [`kernel`](../reference/api/virtualcollision/kernel.md), [`assessment`](../reference/api/virtualcollision/assessment.md), [`exported`](../reference/api/virtualcollision/exported.md)

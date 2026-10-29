# Kn extension: empirical whole-field scaling of Cai's density

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | Fit coordinate-only maps W_Kn that deform Cai's unchanged collisionless density field, N̂(X, P; Kn) = N_C(W_Kn(X, P)), to sampled DSMC density across Kn |
    | **Model** | `cai2012ext.scaling` (map families), `scaling_fit` (separate whole-field and contour fits), `scaling_inventory`, `scaling_report`, `scaling_workflow` |
    | **Evidence** | Unit-tested inventory, fitting and workflow, and a bounded **synthetic** smoke mode. **There is currently no whole-field scientific fit** to the seven-decade suite: its DSMC fields do not exist yet |
    | **Runs on** | Python only (numpy, scipy, PyYAML; matplotlib for figures). Never launches DSMC |
    | **Start in** | `cases/cai2012-kn-extension/` (or pass the script path from the repository root) |

This is an **empirical description**, not a physical model: the maps compose
scalar fields and do not conserve particle number, momentum or energy, and a
good density fit says nothing about force, heating, pressure or temperature.
Map families, objectives, diagnostics and limits are in
[Empirical whole-field scaling](../cai2012-scaling.md).

## Inputs

| File | Role |
|---|---|
| `scaling.yaml` | `registration` (optional), `output` (`results/scaling/full-suite`), `resources` (`cores` ≤ 8, per-input and total byte caps, chunk size, smoke limits), `fit` (tracks F and C, subsets `exploratory`/`suitable`, map families, bounds, regularisation, bins, `max_nfev`, hold-out Kn), `evaluation` (grids, temperature transfer, sensitivity sets) |
| `scaling-registration.yaml` | an explicit, **analysis-only** opt-in for extra result directories (currently two local graded archives under `Cases_run1_graded/`, which a fresh clone does not have — they are then reported missing). Not needed for the active suite |
| `study.yaml`, `manifest-full-suite.yaml` | the active suite; every generated case is registered automatically |
| case fields | reconstructed OpenFOAM fields with actual `C` and `V` and density, or a portable NPZ (schema in the [explanation](../cai2012-scaling.md#inventory-registration-and-source-compatibility)) |

Adequacy for the `suitable` subset requires the existing feasibility criteria,
measured mesh and mean-free-path diagnostics, occupancy above the minimum, a
consistent particle weight, `stable-window-demonstrated` stationarity and
matched coverage. Unknown evidence does not pass. The Kn ≤ 0.001 direct-exit
cases can enter exploratory fits but never the suitable one while their inflow
surface limitation stands.

## Commands

`Allscale.py` pins BLAS/OpenMP to one thread before importing numpy, then calls
`plumetools.cai2012ext.scaling_workflow.main` with this directory as the study
root. With no mode it runs `inventory`.

| Mode | Reads | Writes (under `output`, default `results/scaling/full-suite/`) | Cost |
|---|---|---|---|
| `inventory` | configuration, manifests, which fields exist, diagnostics YAML | `inventory.json`, `inventory.csv` | seconds |
| `prepare` | the inventory | analysis-only copies of case specifications under `specifications/`; refuses to change existing copies | seconds |
| `smoke` | nothing from DSMC: **synthetic** Cai fields under a known axial map, two Kn, ≤ 48 points per case, ≤ 8 optimiser evaluations | `smoke/` — inventory, `assessment.json`, figures, all titled synthetic | bounded |
| `fit --full` | every compatible, registered field that exists | `assessment.json`, per-stratum predictors, figures | unbounded by design; `--full` is required |
| `report` | a stored `assessment.json` (no refit) | the summary `README.md` and figures | seconds |
| `evaluate` | `--predictor predictor.json --points X,P.csv --kn K` | `predictions.csv`, `evaluation-timing.json` | seconds |

Common options: `--config PATH`, `--output PATH`, `--no-plots`; `fit` and
`smoke` also take `--sensitivities` (refit the configured numerical
sensitivity sets).

### Cheap checks

```bash
cd cases/cai2012-kn-extension
python Allscale.py inventory            # what exists, what is missing, and why
python Allscale.py smoke                # synthetic, bounded; output in results/scaling/full-suite/smoke
python Allscale.py report --output results/scaling/full-suite/smoke
```

The smoke mode demonstrates that the numerics recover a known synthetic map.
It is not a fit to any DSMC result and its outputs say so.

### Full assessment

Only once reconstructed fields with `C`, `V` and density exist for the cases
you want to fit:

```bash
python Allscale.py fit --full
python Allscale.py report
```

Missing fields are skipped with reasons; an empty `suitable` subset is reported
as unidentifiable, not filled in. The run is serial, single-threaded, in chunks
of 64 Cai points, with a 512 MiB per-input and 1 GiB aggregate cap; exceeding
the aggregate aborts rather than silently choosing a subset. Raise the caps in
`scaling.yaml` only after checking the machine. On PBS, `qsub runScaling.sh`
(one node, 8 cores, `shortq`, 1 h — the owner's settings) runs
`inventory`, `fit --full` and `report`.

### Evaluating a fitted predictor

```bash
python Allscale.py evaluate \
    --predictor results/scaling/full-suite/smoke/stratum0/exploratory/F/axial/predictor.json \
    --points points.csv --kn 0.1 --output predictions
```

`points.csv` is headed `X,P` with X > 0 and P ≥ 0 (X = x/D, P = ρ/D); the
off-axis disk quadrature is undefined on the source plane. The predictor
inherits the source speed ratio from its fit. The path above exists only after
`smoke`; a real predictor exists only after `fit --full`.

## Outputs

| Path | Contents |
|---|---|
| `inventory.json`, `inventory.csv` | per case: source, Kn, whether fields exist, adequacy evidence and why a case is excluded |
| `assessment.json` | per stratum (source-problem hash), subset, track and map family: fitted coefficients, objectives, diagnostics, timings |
| `stratum*/…/predictor.json` | a serialised `CaiPredictor` |
| `README.md` | the written summary `report` produces, including the limits of the method |
| figures | paired Cartesian and source-centred scientific figures (`scaling_report`) |

JSON is strict: unavailable or non-finite diagnostics are written as `null`,
never `NaN`. Densities are n/n_e (exit number density) and temperatures T/T_e.

**Judging the output.** A better descriptive fit alone does not establish
physical correctness, predictive validity or a self-consistent kinetic
solution. Read the `suitable` track for anything you would cite, and check
which cases it actually contains.

## Resources

`scaling.yaml` declares `cores: 8` (an upper bound — the optimiser is serial)
and the byte caps above. No measured timing of a full fit exists, because no
full fit has been run.

## Figure placeholder

!!! figure "Figure placeholder — fitted map against sampled density"
    **Subject:** sampled n/n_e and the fitted axial map's prediction on the meridional plane, one Kn.
    **Source run:** `Allscale.py fit --full` on real reconstructed fields of at least two candidate cases.
    **Plotting command:** `python Allscale.py report` (figures from `scaling_report`; not run on real data).
    **Destination:** `docs/assets/figures/kn-scaling-fit.png`
    **Suggested caption:** "Empirical coordinate map of Cai's density against DSMC at Kn = … (subset, track and family stated)."

## Related

- [Empirical whole-field scaling, explained](../cai2012-scaling.md)
- [Kn extension, direct cases](cai2012-kn-extension.md) — where the fields come from
- API: [`scaling`](../reference/api/cai2012ext/scaling.md), [`scaling_fit`](../reference/api/cai2012ext/scaling_fit.md), [`scaling_inventory`](../reference/api/cai2012ext/scaling_inventory.md), [`scaling_workflow`](../reference/api/cai2012ext/scaling_workflow.md)

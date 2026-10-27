# Empirical whole-field coordinate scaling of Cai density

This experiment evaluates `N_hat(X,P;Kn) = N_C(W_Kn(X,P))`, where
`X=x/D_e`, `P=rho/D_e`, `xi=hypot(X,P)` and `theta=atan2(P,X)`.
`W` maps a physical comparison point to a Cai sampling point. The original
`plumetools.cai2012.analytical` disk evaluator is unchanged. There is no fitted
density amplitude, normalization, determinant multiplier, or virtual-collision
model. The mappings compose scalar fields; they do not conserve particle number,
momentum or energy. Density agreement establishes no force, heating, pressure,
or tensor-temperature agreement.

Implementation audit revision: `6e013bcf6099eb84828050c2f0b498ff53147176`.
The initial working tree contained unrelated changes and archived results; these
were preserved. No applicable `AGENTS.md` was found. No production DSMC,
scheduler submission, large result copy, or mesh generation was performed.

The follow-up full-suite setup selects `meta.knudsen_suite: full_decades` in
the extension study and generates Kn_e = 100, 10, 1, .1, .01, .001, .0001 into
`Cases_full_suite`, with `manifest-full-suite.yaml`. Scaling reads those active
study paths automatically. The older `Cases` tree and `manifest.yaml` remain
separate; their presence cannot override the selected full-suite inputs.

## Commands and bounded execution

Run from the repository root using an environment with the project dependencies,
pytest, and optional matplotlib (`pip install -e ".[test,plots]"`). On this
Windows workspace the available interpreter is `python`; the repository `.venv`
contains Linux binaries and is not the Windows interpreter.

```powershell
python cases/cai2012-kn-extension/Allscale.py inventory
python cases/cai2012-kn-extension/Allscale.py prepare
python cases/cai2012-kn-extension/Allscale.py smoke
python cases/cai2012-kn-extension/Allscale.py report --output cases/cai2012-kn-extension/results/scaling/full-suite/smoke
```

No arguments means inventory only. `smoke` uses **synthetic** Cai fields under a
known axial map, two Kn values, at most 48 points per case and 8 optimizer
evaluations per fit. It writes beneath `results/scaling/full-suite/smoke`; titles and JSON
identify the data as synthetic. Whole-field fits require the explicit command:

```powershell
python cases/cai2012-kn-extension/Allscale.py fit --full
python cases/cai2012-kn-extension/Allscale.py report
```

These commands inspect every compatible, registered available Kn. Missing input
fields are skipped with reasons, and an empty suitable subset is reported as
unidentifiable. The implementation uses a serial optimizer, one BLAS/OpenMP
thread, chunks of 64 Cai points, no interpolated Cai cache, a 512 MiB per-input
cap, and a 1 GiB aggregate loaded-input cap. Exceeding aggregate memory aborts
the assessment; it must not silently select a convenient subset. Configuration
rejects core declarations above eight. A complete assessment is never the
default laptop action. Input files are checked against uncompressed NPZ size
before loading. Increase explicit limits only after checking HPC resources.

The PBS template inherits `hpc_env.sh`, the study's `shortq` convention, and one
node with eight cores. Submit manually on the cluster; this was **not executed**:

```bash
cd cases/cai2012-kn-extension
qsub runScaling.sh
# Optional predeclared sensitivity refits, also analysis only:
python Allscale.py fit --full --sensitivities
```

`--config PATH`, `--output PATH`, and `--no-plots` work for the CLI modes.
`report` reads stored fits and does not refit. For a portable fitted predictor,
create a CSV headed `X,P` and evaluate without DSMC input:

```powershell
python cases/cai2012-kn-extension/Allscale.py evaluate --predictor cases/cai2012-kn-extension/results/scaling/full-suite/smoke/stratum0/exploratory/F/axial/predictor.json --points points.csv --kn 0.1 --output predictions
```

The predictor inherits the source speed ratio from its fit. Points must lie
downstream with `X>0`, `P>=0`; the off-axis disk quadrature is undefined on the
source plane. The 3D mapping API preserves azimuth and supports explicit solver
axis, nozzle center, and diameter adapters.

## Inventory, registration and source compatibility

Inventory reads study configuration, local generated manifests, actual fields,
and diagnostic YAML. It never uses README run status. `scaling-registration.yaml`
is an explicit **analysis-and-specification-only** opt-in named by `scaling.yaml`.
It registers the local graded Kn10/Kn1 archives in addition to the active suite.
Removing `registration` selects the active study alone. All seven generated
full-suite cases are registered from `study.yaml` and `manifest-full-suite.yaml`
without adding a row to the analysis registry. Reconstructed OpenFOAM fields
are discovered automatically after the future HPC runs.

Generation, meshing and run guards recognize the explicit full-decade profile
carried in the study, manifest and generated case metadata. Absent that profile,
the historical four-case policy and paper-Kn exclusions remain intact. Numerical
feasibility checks are unchanged. The analysis registry is not a production
manifest and grants no separate execution permission.

`Allscale.py prepare` writes analysis-only copies under
`results/scaling/full-suite/specifications`, including feasibility and unresolved
constraints. All new full-suite scaling outputs use `results/scaling/full-suite`
to preserve the initial `results/scaling` assessment. It refuses changed existing copies. It does not generate runnable
case trees; use `generate_cases.py` for the active seven-case suite. Historical
analysis copies may no longer match the new suite; use a fresh `--output` path
instead of overwriting them.

To supply a new study-local result, add a registration row with `name`, `path`,
`Kn`, and optional `field_npz`. Paths must remain inside this study. A portable
NPZ must contain numeric arrays with no pickle objects:

| Key | Shape | Meaning |
|---|---|---|
| `points` | `(n,2)` | `X,P`, dimensionless, actual physical cell locations |
| `density` | `(n,)` | measured `n/n_e`, raw values retained |
| `volumes` | `(n,)` | actual cell volumes in m³ |
| `temperature` | `(n,)`, optional | scalar translational `T/T_e` |
| `occupied` | `(n,)`, optional | meaningful density-estimator mask |
| `domain_valid` | `(n,)`, optional | fluid-domain mask excluding solid/outside points |
| `parcel_count` | `(n,)`, optional | measured occupancy; zero excludes absent estimators |
| `velocity_sqrt_beta0` | `(n,3)`, optional | source-aligned velocity for independent diagnostics |
| `cylindrical_velocity_sqrt_beta0` | `(n,2)`, optional | axial and cylindrical radial components |

The loader also supports original sampled OpenFOAM fields with actual `C` and
`V`; reconstructed density and geometry are required. Plane/profile CSVs lack
these estimators and cannot substitute for full-domain data. Missing data are
never synthesized as zeros. A zero estimator may be retained only with explicit
meaningful occupancy; zero parcel count indicates no estimator. Low positive
sampled densities remain in the exploratory fit and floor-active counts are
reported. Occupancy alone is not a statistical interval.

Direct-exit cases are grouped by a hash of diameter, speed ratio, temperature,
gas law/parameters, geometry and mean-free-path length convention. Different
source problems are fitted separately. Downstream interface sources are visible
and excluded, including unreadable interface configurations reported honestly.
External Cai/health configurations are provenance only, with `included: false`;
their result fields are never read or silently imported.

Adequacy requires existing feasibility criteria, measured mesh/mfp diagnostics,
plume occupancy above the configured minimum, consistent particle weight,
`stable-window-demonstrated` stationarity, and matched full-field coverage.
Unknown evidence does not pass. Every model/track is run separately on all
computable compatible cases (exploratory) and the adequate subset (suitable).
The acknowledged wrong direct-exit inflow surface at Kn_e <= .001 is explicitly
retained for this setup and deferred to future work. Those cases remain available
to exploratory scaling when computable field arrays arrive, but this source
limitation prevents them passing the suitable-reference track even if other
diagnostics pass. No downstream-interface replacement is implemented or mixed
into the direct-exit Kn-only fit.

## Map families and numerical settings

The default empirical Kn feature is `h=1/(1+Kn)`. The small replaceable feature
configuration also permits `inverse_kn_squared` (`1/(1+Kn)^2`); this is a basis
choice, not a collision probability or similarity law. Both saturate at low Kn,
so they cannot distinguish decades arbitrarily deep in the continuum limit.

* Baseline: identity.
* Axial: `log(a_x)=c_x1*h` by default, `W=(X/a_x,P)`.
* Axial/transverse: independent positive `a_x,a_rho`,
  `W=(X/a_x,P/a_rho)`. Scales may contract or expand. Setting `degree: 2` adds
  `h^2` coefficients as an explicit sensitivity experiment.
* Radial/angular: preserve theta and integrate `eta=integral(exp(g),0,xi)`.
  `g=h*sum(c_jl*B_j*A_l)`, with `B=[1,t/(1+t)]` and
  `A=[1,cos(2theta),cos(4theta)]`. `angular_order: 0` or `1` tests truncations.
  Coefficients are ordered by radial basis first, angular basis second.

`a_x<1` compresses displayed axial contours: a physical X samples a larger Cai X.
Affine candidates can change sampling angle; the theta-preserving radial family
does not contain general affine anisotropic scaling. The radial integral uses
deterministic Gauss-Legendre quadrature after `y=log(1+t)`, default order 24.
Tests check convergence across radii `1e-4` to `1e4`. `quadrature_convergence`
compares q, 2q, 4q; disk-evaluator convergence is reported separately. Identity
short-circuits mapping exactly and evaluates Cai directly.

Default global log-scale/log-derivative bounds are ±2; a separate coefficient
box is ±1. Both are serialized. Global extrema are checked without clipping the
map. Invalid candidates receive a fixed-size penalty. Bound-active fits and
optimizer termination are reported; an exhausted iteration cap is not convergence.
The explicit coefficient-size penalty is `1e-4*sum(c^2)`; physical diagnostic
losses exclude it. The configured sensitivity suite includes zero and `1e-3`
regularization, bounds 1/3, map orders 12/48, degree2 and angular truncations.
These are numerical choices, not physics.

For `xi>0`, `J=(eta/xi)^2*exp(g)` is the axisymmetric 3D volume Jacobian;
the affine Jacobian is `1/(a_x*a_rho^2)`. Positive radial derivative prevents
radial folding. Angular dependence in the radial basis at the origin can leave
no unique differentiable origin Jacobian; the API returns NaN and serializes
unavailable diagnostics as null. The map remains continuous with eta(0)=0.
Source consistency outputs explicitly check movement of the physical lip
`(X,P)=(0,.5)` and nearby points. Preserving the exit plane does not preserve the
whole emitting disk boundary. Near-plane evaluator-cap counts flag unresolved
quadrature sensitivity; such computable points remain in the fixed field mask.

## Separate objectives and full-domain diagnostics

F fits log density on matched physical cell locations. Fixed bands come from
baseline Cai density and boundaries `.1,.01,.001`, including both tails. Each
case has equal importance, each occupied angular-band stratum has equal total
weight, and actual cell volumes set relative weights within a stratum. Empty
strata are omitted with counts and coverage. The floor is `1e-10`, fixed before
fitting; diagnostics also check floor multiples and configured `1e-12,1e-8`.
No sample may change strata or leave the comparison mask as coefficients change.

C fits squared `log(xi_hat/xi_DSMC)` at the same three levels. Each case, occupied
level and available angular bin has equal importance. Rays are fixed angular/
log-radius volume estimators, not digitized images or claims about exact samples
on the central ray. Missing shells interrupt a profile. Missing, truncated,
plateau/ambiguous and multiple crossings are flagged; no branch is forced.
Interpolation is linear in log radius/log density within a single observed
crossing bracket. A model losing an established reference crossing is penalized.
The full report evaluates **both F and C for both fitted tracks**.

No historical lip exclusion or density-floor mask is applied. Disjoint excluded
categories account for invalid file values, solid/outside domain, undefined
analytical boundary, absent estimators and undefined baseline. Regions include
near exit, near lip, centerline/off axis, radial bands and angular bands, all
using common physical masks. Volume-weighted normalized L2, bias, absolute
errors and log loss accompany fitting objectives. Cartesian and source-centered
maps have equal space, scales and residual treatments. Raw field CSVs preserve
reference/predicted density, actual volume, coordinates, mapped points and J.

`holdout_kn: []` disables unseen-Kn validation. When enabled, every case at the
selected Kn is separated before DSMC-dependent preparation; training and held
metrics are stored separately. Interpolation/extrapolation is labeled. Numerical
sensitivity settings are predeclared and are not selected using held-out data.

Scalar temperature is evaluated under the same density-derived map with no
coefficient refit. There is no separate temperature fit or velocity correction
in this implementation. Where normalized velocity is available it is compared
independently against unchanged Cai velocity; angular particle flux and
corrected-density/unchanged-velocity **hybrid** flux are labeled diagnostics,
not conservative predictions or a joint kinetic distribution.

Statistical uncertainty is unavailable for averaged fields alone. Optional
`block_resamples` defaults to four bounded draws, skipped immediately when
qualified blocks are unavailable. Set it to zero to disable. Portable metadata must
declare independent block fields `density_blocks`, `independent_blocks: true`,
the `temporal_correlation_method`, `block_duration_s`, and
`integrated_autocorrelation_time_s`. Duration must be at least five integrated
correlation times. The code does not infer independence from parcel occupancy
or treat correlated cumulative means as independent snapshots. Admissible
blocks are jointly resampled across training cases with a fixed seed and bounded
refits; parameter spread is reported separately from numerical sensitivities.

## Local assessment and scientific limits

The initial implementation's local inventory found partial exports at Kn10/Kn1 in the
explicitly registered graded archives and at Kn.001/Kn.0001 in the generated
direct-exit study. None contains the matched whole-field arrays and actual cell
volumes needed by this objective. Thus **all seven target decades currently lack
eligible whole-field results**, and neither exploratory nor suitable scientific
coefficients were fitted. There is no assessment yet of whether this scaling
describes DSMC disagreement. The follow-up full-suite setup generates fresh
inputs in `Cases_full_suite`; missing results remain missing until HPC sampling
is completed. The new suite's manifest is the automatic scaling entry point.

The graded archives have inconclusive stationarity and roughly 70–71% of plume
cells below the configured occupancy floor. The low-Kn configurations retain
their nonproduction warnings: cell/lambda approximately 26.17 and 269.51, and
dt/tau approximately 1.42 and 14.60 at Kn.001 and Kn.0001. A stationary average
does not repair unresolved spatial or collision-time scales. The prepared Kn.01
specification is also budget limited (cell/lambda approximately 2.69). No
requirements were relaxed to make these cases appear feasible.

The bounded synthetic example confirms numerical execution and known-map
recovery only. Its 16 fit records cover all four families, two objectives and
two suitability subsets; the synthetic suitable/exploratory inputs deliberately
coincide. A more flexible family improving descriptive error cannot establish
physical meaning, novelty, causal near-source collisions or unseen-Kn accuracy.

Outputs are separate from legacy results: `inventory.json/csv`,
`assessment.json`, `coefficients.csv`, per-fit `fit.json` and `predictor.json`,
paired Cartesian/polar density/residual maps, contour overlays and displacements,
profiles, source/Jacobian diagnostics, subset comparisons and separate numerical
timings. Production simulation timing is never folded into fitting timing.

## Validation record

Focused numerical tests cover identities/high-K, scaling direction, bounds,
monotonicity, finite-difference 3D Jacobians, independent-point known-map recovery,
crossing flags, equal fixed weights, full-domain accounting, whole-Kn separation,
quality/source filtering, temperature transfer, serialization, missing results,
time-block correlation gates and independent physics diagnostics. Run:

```powershell
python -m pytest tests/unit/test_cai2012ext_scaling.py tests/unit/test_cai2012ext_scaling_fit.py tests/unit/test_cai2012ext_scaling_inventory.py tests/unit/test_cai2012ext_scaling_workflow.py -q -p no:cacheprovider --basetemp cases/cai2012-kn-extension/results/scaling-test-tmp
```

The initial scaling implementation's legacy regression run returned 66 passes,
one failure and 17 setup errors. At that initial audit, `policy.py`, `generation.py`,
`study.yaml` and the failing legacy tests have hashes identical to HEAD.
`test_cai2012ext_policy` expects the parcel-target override to be rejected,
although HEAD's allowlist permits it; `test_cai2012ext_generation` expects an
invariant target while the shipped Kn10/Kn1 overrides are 2000 versus 40 in other
rows. They are documented rather than silently changing existing production
policy as part of this experiment. A broader Cai run (excluding that generation
file and the known policy assertion) completed with 807 passes, four skips,
one deselection and 26 failures: 24 shell tests could not start the Windows
system Bash/WSL (`E_ACCESSDENIED`), and two untouched interface-policy tests
disagree with their existing admissibility verdicts. The new PBS template
passes `bash -n` when invoked directly with Git Bash. This is syntax validation,
not cluster execution. The final focused scaling run passed **71 tests in
8.60 seconds**. The final bounded synthetic fit completed all 16 model/track/
subset records successfully (7.78 seconds assessment, 5.69 seconds inventory and
loading on this machine; report rendering is timed separately). Its axial
coefficient estimates were -0.199943 (F) and -0.199743 (C) for a synthetic truth
of -0.2, including the declared regularization. These are numerical checks,
not measured DSMC improvements. Actual scientific fitting wrote an explicit
zero-eligible-input assessment, and seven configuration pairs were prepared.

Added source files are `scaling.py`, `scaling_fit.py`, `scaling_inventory.py`,
`scaling_workflow.py`, and `scaling_report.py` in `plumetools/cai2012ext`.
Study additions are `Allscale.py`, `scaling.yaml`, `scaling-registration.yaml`,
`runScaling.sh`, this document, and four focused unit-test modules. The study
README links here. Existing baseCase comparison scripts and production policy
sources were not edited during that initial scaling implementation. The later
seven-case setup extends the explicit generation/profile policy and its
verification; its smoke checks are recorded in the extension study README.

The completed seven-case follow-up passed 237 focused policy/generation/HPC/
diagnostics/scaling tests with three Unix-locking skips. Its actual preparation
produced 70 dictionaries, passed 84 independent checks, and verified that all
seven active case paths appear in the scaling inventory. No meshes or DSMC
results were produced. The earlier policy/parcel-target and Windows Bash test
mismatches described above were addressed in this explicitly authorized suite
extension; the separate historical interface-validation tests were outside
this focused run.

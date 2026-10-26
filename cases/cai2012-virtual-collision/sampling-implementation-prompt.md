# Implementation prompt: refine virtual-collision sampling before changing physics

You are working in `https://github.com/andytorrestb/rarefiedPlume`. Implement the
targeted changes below, run focused tests and a bounded numerical pilot, and
report actual results. Do not stop after proposing a plan. This is a numerical
assessment against DSMC, not experimental validation.

## Objective and scope

Resolve numerical sampling and event-measure errors in the existing reduced
virtual-collision model before deciding whether its physical closure improves
free-plume predictions over Cai. Preserve the existing predictions as the
original-sampler baseline. Implement a separately identified revised sampler
with the **same physical closure**, verified VHS kernel, disk source, gas,
geometry, finite domain and comparison conventions.

The user expects existing unresolved DSMC results to retain useful trends, with
approximately 10% resolution sensitivity. This is a provisional sensitivity
assumption, not a measured confidence interval or a target for parameter fitting.
The collision correction can be much smaller than 10% of the field: numerical
error must be measured relative to that correction too. Improvement is a
hypothesis; reproducible degradation or lack of improvement is a valid result.

Do not change DSMC, fit collision strength or distance, introduce local-Maxwellian
partners, implement a time-stepped collision solver, or add further physical
collisions in this pass. Consistent depletion and second-collision closures in
`sampling-refinement-design.md` are subsequent, separately versioned work.

## Audit and existing implementation

Read applicable `AGENTS.md`, inspect `git status` and record the current revision
before editing. Preserve unrelated work, local datasets and previous predictions.
The historical scientific audit used
`ad4cdc636b094b3bf1e0b601e7967b8b927370cc`; the virtual-collision implementation
was subsequently added as focused commits. Verify the current files rather than
resetting to that old revision. Generated case/result trees may be ignored and
available only locally. Do not regenerate missing cases or results automatically.

Read:

- `docs/cai2012-virtual-collision.md`: derivation and collision-law provenance.
- `cases/cai2012-virtual-collision/{README,refinement-plan,sampling-refinement-design}.md`.
- `plumetools/virtualcollision/{source,kernel,moments,model}.py`.
- `plumetools/virtualcollision/{export_geometry,exported,export_figures,export_summary}.py`.
- `cases/cai2012-virtual-collision/{compare_exports,review_results,audit_hazard}.py`.
- `tests/unit/test_virtualcollision_*.py`.

Important existing contracts:

- `DiskSource.sample_flux` samples the disk and outward flux `vx*f_e`, using
  scrambled Sobol source coordinates. `local_quadrature` represents the actual
  ray-restricted ballistic distribution `f0`, not a Maxwellian reconstructed
  from moments. Preserve these distinctions.
- `VHSKernel` matches the installed OpenFOAM v2512 gamma-normalized collision
  kernel and isotropic scattering. Keep its independently checked reference
  value and parameters. The exact Hou 2025 reference has not been identified;
  do not claim this implements that published model.
- Current `MODEL_VERSION` is `mutual-first-event-1`. The symmetric continuum
  pair-event measure is

      dR = 0.5 f0(x,v) f0(x,w) K(|v-w|)
                 exp[-H0(x,v)-H0(x,w)] dx dv dw,

  where `K(g)=sigma(g)*g` and `H0` is ballistic exposure from the emitting disk
  to the event. Both incoming populations lose a trajectory; both products
  stream freely after the event. This is a mutually thinned first-event
  approximation, not a complete first-order Boltzmann integral or DSMC solver.
- Disabled collisions and `ZeroKernel` reproduce the existing Cai evaluator
  exactly at the original requested coordinates. Preserve that regression.
- Raw `M0`, `M1`, `M2` precede central moments. Preserve raw negative/invalid
  predictions and report flags. No clipping, renormalization or fallback to Cai.
- Enabled ring transport uses local cylindrical `(x,r,azimuthal)` moments.
  Recovered Cartesian cells, finite tube selection and angular footprint
  fractions are verified by `export_geometry.py`. Cylindrical radial velocity
  and spherical outward velocity are different observables.

## 1. Correct the discretized event target and proposal

`model.py::_hazard` currently uses quadratic time panels and midpoint collision
frequencies. `preprocess` samples a conditional event location from this panel
hazard, selects partners with the event-point frequency, then recomputes partner
survival on panels rescaled to a different endpoint. The source survival prefix
and partner endpoint survival are therefore not the same numerical functional.
Pairwise conservation does not detect this event-distribution error.

Write a short derivation **before implementing the correction**. Define a shared
endpoint-hazard evaluator `H0(x,v)` for either incoming particle, with a stated
quadrature/error control and source visibility conditions. Distinguish the
proposal frequency `nu_p` and prefix `H_p` from the target event-point frequency
`nu0` and common endpoint exposure `H0`.

For a conditional proposal `q_t=nu_p*exp(-H_p)/P`,
`P=1-exp(-H_p,end)`, and partner proposal `f0(w)*K/nu0(x,v)`, the pair importance
weight for the stated target is

    a = 0.5 P (nu0/nu_p) exp[H_p-H0(x,v)-H0(x,w)].

Derive this from the actual implemented discretized measure, including source
flux and residence Jacobians; do not copy the expression without checking those
assumptions. A frequency-ratio correction alone omits the survival inconsistency.
A consistent alternative quadrature construction is acceptable if derived and
tested against the same target. Use stable log weights. Diagnose proposal support,
zero rates, large weights, effective sample size and nonfinite evaluations. Do
not cap weights or normalize away an event-budget failure.

Test a nonconstant analytically integrable hazard, source/partner exchange
symmetry and convergence at 24/48/96 panels. Keep the original numerical sampler
available for explicitly labelled A/B comparisons. Original proposal probabilities
must not be presented as corrected physical ensemble survival probabilities.

## 2. Reduce outgoing-direction variance

For each sampled pair, begin with three orthogonal relative-velocity axes in a
uniformly random rotated frame, each with weight `1/3`. Retain both equal-mass
elastic products along every axis. This integrates the isotropic local pair
second-moment change exactly:

    E[Delta(vv^T + ww^T)] = 0.5 [(g^2/3) I - g_vector g_vector^T].

Verify this against an independent analytic benchmark. Spatial residence is not
integrated exactly by three axes; assess 3/6/12 distinct axes, for example by
averaging independently rotated orthogonal frames. Opposite signs merely swap
the two identical outgoing products and do not add information.

Angular branch weights sum to the original pair-event weight. Subtract the two
incoming downstream trajectories once with that total weight, not once per
direction with full weight. Preserve the physical isotropic scattering law.
Measure variance times cost; do not assume more directions are always cheaper
than more source samples.

## 3. Stratify possible first-event locations

Support two and four equal-probability strata in the conditional first-event
CDF, using a randomized quantile in each. Include stratum mass in each event
weight. These are alternative quadrature branches for **one possible physical
collision**, not two/four successive collisions. Sample the conditional partner
from `f0(w)*K/nu0` at each location and apply its location-specific survival and
full target/proposal correction.

Make source samples, location branches, direction branches and physical pair
event weights explicit. Replace the current implicit assumptions of four rays
per event and two remaining-exposure entries per event. Update every affected
diagnostic: event-rate/emission ratio, ensemble loss/survival, exposure summaries,
remaining-exposure weighting and boundary fluxes. Branch refinement must not
multiply the represented emission or collision count. Maintain paired number,
momentum and energy source/sink balance and count return/backscatter crossings.

Use deterministic random streams keyed by source/batch/branch purpose so that
changing chunk size or direction count does not accidentally change source and
event-location draws. Keep reproducibility and useful common random inputs for
controlled A/B differences; use independent scrambles to assess uncertainty.

## 4. Repair uncertainty and diagnose support

Current 2,048-source runs use eight batches and a 1,760-ring interpolation grid:
only 256 trajectories per batch. Many rings have zero or very sparse signed
corrections. Derived `U` and `T` from near-zero individual batch densities can
have ratio poles. Zero sampled correction with zero estimated SE is not evidence
of zero numerical error.

Retain raw batch moments and their covariance. Implement documented uncertainty
for derived fields using independent full-estimate scrambles or a justified
leave-one-batch-out/raw-batch resampling method. Recompute nonlinear profile
estimators from the resampled raw arrays. Report density conditioning, resample
failures and unresolved support; suppress misleading Gaussian bands with an
explicit uncertainty-unavailable flag, while retaining the underlying raw data.
Do not clip samples to force finite ratios or use occupancy as DSMC uncertainty.

Preserve `exported.py`'s matched Cartesian tube and angular-footprint estimators.
Legacy molecule-weighted averages of local-cell temperature and temperature
derived after raw-moment aggregation are distinct; retain their names. An ideal
smooth cylinder/sphere may be an additional diagnostic, never a silent replacement.
Do not refine the interpolation grid before establishing adequate correction
support. Use the same saved common mask for baseline and revised errors, including
invalid-result counts; do not let each prediction select where it is assessed.
Keep legacy metric semantics and add explicitly named DSMC-reference metrics.

## 5. Configuration, output preservation and cost

Expose seed, batches, source count, hazard controls, partner quadrature orders,
location strata, angular directions, original/revised sampler and an output
destination through the existing driver/config conventions. Validate numerical
settings. Bump model/cache schema versions and include every new numerical
control plus source, geometry, gas and kernel in cache keys and saved metadata.

Write pilot outputs to a new versioned destination rather than overwriting
`results/virtual_collision/` or the saved original raw arrays. Record source
revision, case paths, full options, hashes of reference inputs where practical,
and cache state. Update CLI help and the existing manual PBS template for the
new controls if necessary; do not submit a job or add an LLM dependency.

Measure preprocessing and field/profile evaluation wall time separately,
core-seconds, peak memory where supported, source/event/branch counts and bins
evaluated. Compare original and revised samplers on identical requested outputs.
Report accuracy, uncertainty and variance-times-cost together. DSMC postprocessing
time is not production simulation cost; unavailable production cost stays
unavailable.

## 6. Bounded pilot and reference handling

Start with focused tests and a tiny model-only smoke. Then run a bounded pilot
on Kn=10 and Kn=1, initially 2,048 sources, comparing original, corrected event
measure, direction averaging, location stratification and their combination.
Use a staged/coarse-observable matrix to avoid multiplying full-resolution
exports by every setting. Compare 8,192 sources only for the informative variants;
32,768 requires a measured cost justification and an explicit resource bound.
Use common inputs for differences and at least two independent scrambles for
the selected comparison, clearly reporting the limited uncertainty evidence.
Keep local parallelism at most eight cores (one thread by default), memory below
32 GB with chunking, and document a wall-time/run-count bound before execution.
Report incomplete or failed runs; do not manufacture a convergence conclusion.

Historically available compatible direct-exit exports are:

- `cases/cai2012-kn-extension/Cases_run1_graded/Kn10_np1x/results/`
- `cases/cai2012-kn-extension/Cases_run1_graded/Kn1_np1x/results/`
- `cases/cai2012-kn-extension/Cases/Kn0p001_np1x/results/`
- `cases/cai2012-kn-extension/Cases/Kn0p0001_np1x/results/`

Recheck availability and actual geometry/Kn metadata. New uniform `Cases/Kn10_np1x`
and `Cases/Kn1_np1x` lacked the three compact result CSVs at the last audit and
are not controlled refinement replicas of the graded cases. Do not mix them.
Dense direct-exit Kn=0.001/0.0001 remain non-production breakdown tests, not
headline accuracy references. Interface-source cases are a different boundary
problem and remain excluded. Preserve the extension family gates excluding
Kn=100, 0.1 and 0.01; discover those in the Cai family if data exist. Kn=0.01 is
additional to the original six-point inventory. Skip missing cases/quantities
with machine-readable reasons. Do not launch DSMC or generate cases.

Saved original plots are indexed at
`cases/cai2012-kn-extension/results/virtual_collision/README.md`; review metrics
and the hazard audit are in its `review/` subdirectory. Existing aggregate
density/velocity/flux errors mostly worsen, despite selected angular-temperature
improvements. Kn=10 plane density normalized error was about 1.443% Cai versus
1.535% original virtual; Kn=1 about 1.483% versus 4.435%. Do not set these as
acceptance targets. Some apparent improvements are below numerical uncertainty;
Kn=1 has moment-invalid angular bins. Recompute scores without hiding them.

## Tests and completion criteria

Add meaningful tests for:

1. Nonconstant-hazard target/proposal integration, incoming exchange symmetry
   and numerical refinement; conservation alone is insufficient.
2. Analytic isotropic second-moment recovery and equal-mass elastic conservation.
3. Location/direction branching normalization and closed-control-volume number,
   momentum and energy balance without double-counting events or emission.
4. Exact disabled Cai and zero-kernel regressions, weak-collision scaling,
   reproducible chunking and unchanged source draws under angular refinement.
5. Raw-moment uncertainty with known distributions, ill-conditioned ratio and
   zero-support flags, and retained invalid arrays without fallback.
6. Cache sensitivity to new controls, common masks, matched profile estimators,
   missing-data skips and unsupported-law flags.

Run the existing `test_virtualcollision_*.py` suite and affected Cai/extension
tests. Historical verification had 44 virtual-collision tests passing. Three
legacy failures reproduced independently at the old pristine source revision:
the policy `target_particles_per_cell` override assertion and interface manifest/
guard `CONFIGURATION_INCOMPLETE` assertions. Recheck and report current outcomes;
do not weaken policy gates to conceal them.

Deliver working code and tests, derivation/diagnostic documentation, a bounded
original-versus-revised pilot with raw arrays and reproducible tables/plots,
timing/convergence/conservation/support reports, exact laptop and manual HPC
commands, and focused commits. Keep density contours, centerline n/U/T, angular
n/U/normalized flux and supported tensor diagnostics. Distinguish model numerical
uncertainty, unavailable DSMC uncertainty and the provisional 10% reference
sensitivity. Finish with actual improvements/degradations, failed/skipped cases,
and whether numerical evidence is sufficient to proceed to a separately derived
physical-closure variant. Do not declare usefulness solely because plots look
smoother or because one error decreased.

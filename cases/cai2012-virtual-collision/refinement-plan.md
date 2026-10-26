# Interpretation and refinement plan, 2026-10-26

**Current judgement:** the conservative event-transport implementation is working,
but predictive usefulness over Cai has not been established. Most aggregate
density, velocity and flux errors increase. Angular temperature is a promising
specific lead, especially at Kn=10, but numerical noise and reference limitations
prevent treating it as a validated improvement. The next stage should isolate
numerical error, reference error and closure error before changing physical
assumptions or adding a fitted collision-strength parameter.

This review uses the saved 2,048-source-sample runs and their original common
masks. It adds postprocessing and a bounded source-only hazard audit; it changes
no physical model, case configuration, DSMC run or saved prediction. Source
revision remains `ad4cdc636b094b3bf1e0b601e7967b8b927370cc` plus the delivered
working-tree implementation.

## What the present data say

Density and flux entries below are volume/solid-angle weighted normalized L2
errors against DSMC, expressed in percent. Temperature is weighted RMS error
divided by exit T0, also expressed in percent. These are different error scales,
not quantities to average into a single score. All remain exploratory.

| Observable | Kn=10: Cai to virtual | Kn=1: Cai to virtual |
|---|---:|---:|
| Meridional density | 1.443% to 1.535% | 1.483% to 4.435% |
| Centerline density | 1.123% to 1.278% | 1.155% to 4.565% |
| Angular density | 0.419% to 0.879% | 2.156% to 5.230% |
| Angular particle flux | 0.543% to 0.802% | 1.344% to 4.586% |
| Angular local-cell temperature | 1.436% to 1.360% | 2.234% to 1.570% |

The apparent Kn=1 angular-temperature gain includes four moment-invalid common
bins. It cannot presently count as evidence of a physically useful prediction.
No invalid bins were removed to improve these aggregate scores.

There are local leads worth testing. At Kn=10, angular temperature error falls
46% in 0-15 degrees and 44% in 15-30 degrees; in absolute units these are about
1.32 K to 0.71 K and 0.83 K to 0.46 K. Density and particle-flux errors fall about
26% and 24% over 30-60 degrees, but the respective correction-to-reported-model-SE
RMS ratios are only 0.86 and 1.10. Their sign and size need confirmation with new
independent samples and numerical refinement. These fixed diagnostic angle/X/D
regions were applied to every case and observable; they are exploratory findings,
not a newly selected headline validation mask.

For Kn=10 angular flux, 65% of the scored solid-angle measure has a smaller
pointwise absolute error, yet the aggregate error increases 48%. A broad region
of small improvements can coexist with larger weighted errors elsewhere. This
explains why a contour or selected profile segment can look better while the
overall score gets worse. Retain both diagnostics instead of selecting one.

The complete [aggregate review](../cai2012-kn-extension/results/virtual_collision/review/metrics_review.csv),
[regional review](../cai2012-kn-extension/results/virtual_collision/review/regional_review.csv),
and [error-ratio figure](../cai2012-kn-extension/results/virtual_collision/review/error_ratios.png)
preserve degradations, masks, reference quality, model SE and invalid counts.

## Why the present comparison is inconclusive

**Sparse signed estimates.** There are 2,048 source trajectories across eight
batches and 1,760 spatial rings: only 256 trajectories per batch. At Kn=10,
372 rings have zero density correction in all batches and 516 have nonzero
correction in at most two batches. These are support diagnostics, not literal
event counts. A zero sampled correction and zero reported SE can mean an
unresolved rare contribution, not an exactly collisionless result. A finer grid
at the same source count can worsen this issue.

On the actual common centerline mask at Kn=10, density-correction RMS is 0.00203
n0 and reported model SE RMS is 0.00213 n0. The temperature correction/SE ratio
is 0.84. At Kn=1, occasional batch densities approach zero and the derived
velocity/temperature ratios become ill-conditioned. Huge bands are a warning
about the estimator; they are not reliable Gaussian confidence intervals.

**Finite-panel event inconsistency.** The location CDF uses a midpoint hazard
frequency, while conditional partners are normalized at the sampled event.
At finite order the sampled pair measure therefore carries an additional
nu_mid/nu_event factor. This is a numerical approximation that disappears under
refinement. Pairwise conservation alone does not test it.

A six-run source-only audit (128 trajectories per run, one thread, 39.5 s total)
found the following event-weighted fraction with frequency ratio differing
from unity by more than 10%:

| Hazard panels | Kn=10 | Kn=1 |
|---:|---:|---:|
| 24 | 29.3% | 22.6% |
| 48 | 7.5% | 8.5% |
| 96 | 2.2% | 0% |

Aggregate exposure and event-rate changes were only about 0.6-0.7%. This does
not explain all the error degradation. It identifies a necessary consistency
check, not a proven cure. The
[audit output](../cai2012-kn-extension/results/virtual_collision/review/hazard_audit.json)
retains settings and distributions.

**Reference error matters at the signal size being sought.** About 96.7% of
Kn=10's aggregate centerline baseline density squared error comes from its first
point. That point has a 4.82% density deficit; the second layer is nearly nominal.
The installed v2512 source inserts particles before movement, gives them a full
timestep and samples fields after movement. A first-cell/end-step residence
effect of approximately U dt/(2 dx)=4.85% is a strong explanation to test.
This is not yet a verified cause for the historical executable. Do not tune the
collision correction to compensate for this possible source/timestep artifact.
Removing the first point also does not make the current correction succeed:
the fixed 0.5-2, 2-5 and 5-10 D centerline regions still degrade.

Stationarity is inconclusive in Kn=10's 4-6 and 7-10 D centerline regions and
Kn=1's 7-10 D region; judged shell sectors are stable. There are 37 saved
instantaneous regional density samples, so some limited blocking analysis is
already possible. Maximum block-SE estimates range approximately 0.2-1.8% for
centerline regions and 0.04-0.22% for shell regions. Only four blocks remain at
the largest tested size and no convincing far-centerline plateau is established.
These are regional diagnostics, not confidence intervals for profile bins or
for the complete 9,401-step average. Per-bin U/T/tensor block histories remain
missing.

**The one-event regime is narrow.** The event-weighted mean of
1-exp(-remaining frozen exposure) is 2.67% at Kn=10, 24.31% at Kn=1, and about
99.94% in both dense cases. This is an omitted-collision-opportunity proxy in
the frozen background, not a measured probability of actual further collisions.
It makes Kn=10 the cleanest starting point and Kn=1 a marginal test.

Kn=0.001/0.0001 also have inadequate DSMC references: exit cell/mean-free-path
ratios about 26.2/270 and dt/collision-time ratios about 1.42/14.6. They belong
in breakdown figures. More averaging or model samples will not fix this.
Moreover the mutually thinned closure retains approximately 54% ballistic
source-number weight in these optically thick runs. In the simple limit with
equal incoming hazards, total modeled source loss tends to
integral nu exp(-2 nu t) dt = 1/2. Its dense-limit saturation is structural.

## Recommended next stages and decision gates

1. **Freeze physical parameters and make the estimator auditable.**
   Expose seed, batch count and the new quadrature orders in a reproducible
   configuration/CLI; current export CLI does not expose independent seed runs.
   Record event/support counts, effective sample information, zero-hit flags and
   density conditioning before dividing raw moments. Keep failed raw results.
   Audit the hazard proposal/partner normalization and its refinement; do not
   insert an unverified corrective weight. Remove the small azimuthal Cai-cache
   discrepancy from the plane comparison. Apply the separate angular flux mask
   explicitly, carry the review's added plane-velocity scores into the normal
   report, and keep these named metric additions separate from legacy results. Compare flux from the
   interpolated fields against the event-crossing control-volume benchmark;
   conserved signed rays do not certify an interpolant's flux.

2. **Resolve Kn=10 numerically before widening the physical study.**
   Use nested source counts 2,048, 8,192 and 32,768 at fixed batch count, first
   on fixed broad observables. Test hazard panels 24/48/96, partner quadrature,
   ring resolution and domain separately, using paired random inputs for
   refinement differences and independent scrambles for repeatability. Increase
   support before using finer rings. Compare changes in the correction and error
   difference, not only the full field that is dominated by Cai. This is a
   numerical ladder, not a physical fit or an assumption that 32,768 is enough.

   Average several **distinct** outgoing scattering directions per expensive
   event, and test 2-4 stratified location quantiles. These remain quadrature of
   a single-event distribution. The paired loss must retain the same total
   weight. Opposite directions merely swap identical products and add no
   information. A new angular-moment benchmark can use
   E[delta(vv^T+ww^T)] = [(g^2/3)I - gg^T]/2 for fixed incoming velocities.
   Three orthogonal relative directions with equal weights match that local
   second moment; spatial transport still needs angular convergence.

   The current implementation already estimates a correction about analytical
   Cai, so adding a Cai control variate is not a new fix. Conditional averaging
   and paired comparisons are supported by standard
   [variance-reduction principles](https://artowen.su.domains/mc/Ch-var-basic.pdf);
   their variance-times-cost benefit here must be measured. Near-equilibrium
   [deviational Boltzmann methods](https://web.mit.edu/ngh/www/LVDSMC.pdf)
   provide useful ideas, not a replacement for this non-Maxwellian plume closure.

   **Gate:** raw moments and derived observables must be realizable and stable
   under independent refinements at an uncertainty smaller than the effect being
   assessed. If this cannot be reached economically, the current estimator is
   not useful at the requested resolution even if its continuum closure might be.

3. **Improve the reference with controlled diagnostic runs.**
   First test the first-layer deficit by halving dt at unchanged source mesh and
   particle weight, or by a separate uniform-injection-time boundary benchmark.
   A confirmed timestep trend would justify a source/sampling correction, not a
   fitted plume-collision parameter. Then use matched-geometry collision-on/off
   auxiliary DSMC cases to compare increments: virtual-minus-Cai against
   DSMC-on-minus-DSMC-off. Preserve the official collision-on/Cai comparison and
   all policy gates. Save compact non-overlapping per-bin sums of number,
   momentum and second-moment trace, optionally all tensor components, plus
   emission/downstream/side/return crossings. Choose block lengths using
   autocorrelation and stationary-window checks.

   The new uniform Kn=10/1 setups are not controlled refinement evidence: they
   change from 10-mm to 20-mm exit cells, 316 to 80 nozzle faces, area ratio
   1.005859 to 1.018592, target 40 to 2,000 exit parcels, and double dt while
   refining outer cells. Vary mesh, dt and particle weight separately.

4. **Confirm a useful weak-collision window on new evidence.**
   Predeclare density, axial velocity and particle flux as safeguards, and test
   the angular-temperature benefit at 0-30 degrees as a specific hypothesis.
   Use new seeds and adequate references. If needed, add controlled Kn≈3 and 5
   in a separate study without weakening extension gates; their expected frozen
   incoming mean exposures are roughly 0.15 and 0.09 for this fixed source/gas.
   Exposure diagnostics, not Kn alone, define the tested regime.

   **Gate:** an improvement must persist after model numerical uncertainty,
   DSMC uncertainty and discretization sensitivity are accounted for. Report
   absolute error reduction and runtime, not just relative improvement from an
   already tiny baseline error. A temperature-only benefit is a legitimate
   restricted use case if its uncertainty/cost and other-field limitations are
   explicit. Consistent degradation after convergence is also a valid conclusion.

5. **Decide on closure changes only after these gates.**
   Keep the one-event model for a measured weak-collision niche if evidence
   supports it. If repeated collisions are material, document the boundary of
   this approximation before deriving a separate limited multi-event closure.
   Do not silently turn it into DSMC, a full collision-integral solver, or an
   empirically rescaled Cai correction. Keep dense cases as breakdown checks.

## Immediate next deliverable

A controlled Kn=10 convergence/variance report and a short source/timestep
verification study would be more decisive than extending the current dense-case
plots. Kn=1 should accompany it as a marginal-regime stress test. No production
DSMC or scheduler jobs were run for this review.

Reproduce the saved-data review and bounded numerical audit from the repository
root in the configured environment:

```bash
python cases/cai2012-virtual-collision/review_results.py
python cases/cai2012-virtual-collision/audit_hazard.py
```

The review tables retain all common-mask points, including invalid predictions.
The optional gap-greater-than-two-reported-model-SE diagnostic is not a confidence
interval and omits DSMC/systematic uncertainty. It must not become a new
prediction-dependent scoring mask.

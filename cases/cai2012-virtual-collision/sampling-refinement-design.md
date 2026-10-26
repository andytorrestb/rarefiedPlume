# Proposed sampling refinements and subsequent closure tests

These are proposed changes, not implemented model changes. The existing
predictions remain the reproducible baseline. Current DSMC cases can guide
broad trends under the user's provisional approximately 10% resolution
sensitivity assumption. That assumption is not a measured confidence interval;
the numerical sampler must resolve the collision increment, which may be much
smaller than 10% of the full field.

## First implementation pass: same physical closure

1. **Use a consistent event measure.** Define one endpoint-hazard evaluator
   H0(x,v) for either incoming particle. The present source CDF uses midpoint
   panel frequency while partners use the event-point frequency; the source
   survival prefix and partner survival also use differently rescaled panels.
   A frequency-ratio correction by itself is incomplete.

   For conditional time proposal q_t=nu_p exp(-H_p)/P, partner proposal
   f0(w)K/nu0(x,v), and intended symmetric pair measure
   (1/2)f0(v)f0(w)K exp[-H0(x,v)-H0(x,w)], the importance weight is

       a = (1/2) P (nu0/nu_p) exp[H_p-H0(x,v)-H0(x,w)].

   Derive the actual discretized target/proposal before choosing between this
   reweighting and a consistent quadrature construction. Evaluate weights in
   log space, verify proposal support and finite moments, and report large
   weights. Do not cap or renormalize a failed estimate. Require a nonconstant
   analytically integrable hazard benchmark, incoming-particle exchange checks,
   and 24/48/96-panel convergence. Conservation remains necessary but does not
   verify the event distribution.

2. **Integrate several outgoing directions per sampled pair.** Start with three
   orthogonal relative-velocity axes in a randomly rotated frame, each with
   weight 1/3; retain both products. This exactly integrates the isotropic local
   second-moment change for fixed incoming velocities. It does not exactly
   integrate spatial residence, which needs angular refinement. Test 3/6/12
   distinct axes and measure variance times cost. Opposite directions merely
   swap identical outgoing products and add no information. Keep the verified
   isotropic VHS physics unchanged.

   The outgoing branch weights must sum to the original pair-event weight.
   Subtract each original downstream trajectory with its total event weight;
   do not multiply the loss or event count by the direction count.

3. **Stratify the possible first-event location.** Test two and four equal-mass
   conditional-hazard strata, using a randomized quantile in each. Include the
   stratum probability in every gain/loss weight. These are alternative
   quadrature branches for one possible collision, not multiple collisions on
   one trajectory. Partners remain weighted by f0 sigma(g)g at each location.

4. **Improve support and ratio uncertainty.** Compare source counts 2,048 and
   8,192 first, with 32,768 only if needed. Hold geometry, physical parameters
   and batch count fixed. Add explicit seed and quadrature CLI/configuration
   controls. Compare variance reductions at matched runtime rather than assuming
   that multiplying all quadrature orders is efficient.

   Report correction support, effective weights and raw-moment uncertainty.
   Use independent full-estimate scrambles or leave-one-batch-out/raw-batch
   resampling for derived quantities, with density-conditioning checks. Ratios
   from eight sparse 256-source batches can have poles and do not give reliable
   normal error bands. Zero observed corrections with zero SE require an
   unresolved-support flag.

   Preserve the reconstructed Cartesian-cell and angular-footprint estimators.
   A smooth ideal tube/shell estimator is a separately named numerical
   diagnostic, not a replacement for the legacy exported estimator. Similarly,
   raw-aggregate temperature differs from the legacy average of local-cell T.
   Refine the ring grid only after its correction is adequately supported.

Benchmark the original sampler, direction averaging, location stratification,
and their combination on Kn=10 and 1. Keep the physical kernel, gas, source,
domain and output masks fixed. Use common random inputs for refinement
differences and independent scrambles for reproducibility. Retain all failures.
Judge both the direction/size of the mean correction and the numerical error;
cleaner error bands alone are not a closer physical prediction.

## Next physical variant: consistent surviving ballistic population

After sampling is stable, a separate no-fit closure can let the particles that
have never collided deplete consistently:

    fb = f0 S,
    dS/dt = -S nu_b,
    nu_b(x,v) = integral fb(x,w) K(|v-w|) dw.

Only fb-fb pairs collide. Each pair contributes both elastic products, which
then stream freely. The pair measure is (1/2)fb(v)fb(w)K. Its two-particle loss
matches the attenuation equation, so event/source balance follows without an
inexhaustible partner reservoir or fitted strength coefficient.

After solving this stationary survival problem, the conditional event density
is nu_b S/P_b, the partner density is fb K/nu_b, and the pair weight is P_b/2.
Do not also multiply by the old partner-survival factor: it is already in fb.
Stationary fixed-point iteration needs convergence/residual checks and possibly
damping or brackets; naive iteration is antitone and can oscillate.

This addresses the current dense-limit retention of approximately half the
ballistic population. In a homogeneous benchmark, S=1/(1+H0) tends to zero.
It still omits collisions of already scattered particles. Its weak-collision
difference from the current closure begins at third order, so it is unlikely to
repair Kn=10 sampling errors. Stronger corrections may improve or worsen DSMC
agreement; compare it as a separate model version.

For Kn=1, a later limited second-collision variant may be more consequential,
but must derive losses for both participating populations and preserve paired
source/sink accounting. More than one physical collision is a new closure,
distinct from taking more quadrature branches of a single event. Dense cases
remain breakdown tests for a few-event model.

## Comparison with the unresolved DSMC data

Use the existing cases to assess whether density redistribution, angular flux,
axial acceleration and cooling move in reproducible directions. Retain exact
scores and add a separately named approximately 10% reference-sensitivity
scenario if desired; do not treat its width as independent random error bars.
Keep VHS coefficients and exit parameters fixed. The objective is a physically
derived, converged correction; closer agreement is measured after the change.

The first concrete deliverable should be an original-versus-revised sampler
comparison at Kn=10/1, including raw moments, uncertainty/support flags,
realizability, closed-volume flux checks and matched-output cost. Physical
closure variants should follow once that comparison is numerically resolved.

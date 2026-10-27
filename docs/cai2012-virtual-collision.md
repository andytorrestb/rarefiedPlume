# Reduced virtual collisions: mathematical contract

Source revision audited: `ad4cdc636b094b3bf1e0b601e7967b8b927370cc`.
This derivation precedes implementation. This is a numerical assessment against
DSMC, not experimental validation or an implementation of Hou's published model.

The source is the existing uniform disk in x=0, emitting into x>0 with the
drifting Maxwellian f_e and flux measure v_x f_e dA dv. The ballistic f0 is
f_e restricted to backward rays intersecting the disk; no local Maxwellian
moment closure is used. The finite collision/transport domain is a half-cylinder
0<x<L, r<B. Crossings of its side, downstream cap, and return plane are exits.

For a ballistic ray define H(x,v)=integral_from_exit^x nu0 dt, with
nu0=integral f0(x,w) K(|v-w|) dw and K=sigma(g)g [m^3/s]. The proposed reduced
pair event measure is

    dR = (1/2) f0(x,v) f0(x,w) K(g) exp[-H(x,v)-H(x,w)] dx dv dw.

We sample source trajectories from the *flux* distribution, integrate their
frozen hazard, and select one location from the conditional first-event CDF
(1-exp[-H(t)])/(1-exp[-H_end]). Its weight is 1-exp[-H_end]. Conditional partners
are drawn with weights f0 K/nu0; an additional exp[-H_partner] weight and the
factor 1/2 yield the symmetric measure above. This is a mutually thinned,
at-most-one-event approximation to each original trajectory, not an effective
multiple-collision model. It omits collisions involving already scattered
particles. Frozen survival thinning also omits some permissible first events.

At every representative event we transport four signed rays: +v', +w', -v,
-w, all starting at the same event. Pair scattering conserves v+w and
|v|^2+|w|^2. The 1/2 counts unordered pairs once. Both original particles are
removed, so the partner is not an inexhaustible reservoir. In the continuum
closure, loss on a source ray is bounded by integral nu0 exp(-H) dt <=1.
Finite sampled signed estimates need not remain positive: raw values are kept
and nonpositive density or nonphysical tensors are flagged, never repaired.

Number, Cartesian momentum and energy corrections integrated over **all** exits
of this closed cylinder cancel event by event. Forward-only flux need not equal
emission. These identities do not establish local positivity, numerical
convergence, or validity when remaining exposure on scattered rays is large.

Steady density moments use signed residence time divided by bin volume; flux
uses crossings. Cylindrical ring bins use volume pi(r_hi^2-r_lo^2)(x_hi-x_lo).
Raw M0, M1, M2 are accumulated before deriving U=M1/M0,
T=(m/k_B)(M2/M0-U outer U), and scalar T=trace(T)/3. The moment basis is
(x,cylindrical radial,azimuthal); conservation uses Cartesian vectors. Rotating
each trajectory contribution before bin aggregation is necessary.

Numerical controls (samples, hazard panels, directional/speed quadrature,
finite domain and bin sizes) are not fitted physical parameters. Independent
scrambled source/event samples supply batch standard errors; these quantify
sampling only. Domain, hazard and partner quadrature errors require separate
refinement. DSMC uncertainty requires time information and is never inferred
from occupancy. No accuracy claim is permitted before those checks and reference
quality assessment. Remaining exposure is measured inside the finite domain;
outside-domain truncation requires an enlarged-domain comparison.

## Source, Jacobians and implementation

The existing Cai source has diameter D=0.2 m, S_e=U_e/sqrt(2kT_e/m)=2 and
T_e=300 K. Its n_e is the full exit Maxwellian density before restriction to
v_x>0. For arbitrary supported source parameters the distribution is

    f0(x,v) = n_e (beta/pi)^(3/2) exp[-beta |v-U_e e_x|^2] chi_exit,
    beta = m/(2 k_B T_e).

The local partner quadrature integrates actual visible-disk directions and speed:
d^3v=c^2 dc dmu dphi, mu=cos(theta). Disk intersection lengths set the mu bounds.
Equivalently, a disk-area parameterization has dOmega=x dA/|x-q|^3; direct solid
angle integration avoids the narrow area-integrand peak near the exit. No local
Maxwellian is reconstructed from moments. Speed integration ends at
c_max=sqrt(2kT_e/m)(S_e+speed_tail_thermal); the default tail is eight thermal
speeds beyond drift. Refining this cutoff is a numerical truncation check.
The omitted ballistic density is bounded absolutely by n_e times
erfc(b)+2 b exp(-b^2)/sqrt(pi), b=speed_tail_thermal (1.46e-27 at b=8).
This centered-Maxwellian tail bound does not bound collision-frequency or
corrected-field quadrature errors.

Emission samples have uniform disk area, transverse Maxwellian velocities, and
axial probability density proportional to v_x exp[-(v_x/a-S_e)^2] for v_x>0.
This v_x factor is essential. The steady change of variables from emission
location and flight time obeys f0 dx dv = v_x f_e dA dt dv; it supplies the source
flux factor used in event rates and residence-time moments. Axial inverse-CDF
sampling uses a 32769-point interpolation and a finite tail of at least nine
thermal speeds beyond drift. Independent scrambled Sobol batches reduce source
variance; representative location, partner and scattering samples use seeded
pseudorandom draws. The finite interpolation and velocity quadrature introduce
separate deterministic errors that batch standard errors cannot measure.

The hazard is integrated in midpoint panels with edges t_j=t_end(j/N)^2,
concentrating panels near the source. A representative event is sampled
from the resulting piecewise-constant hazard CDF, then its partner is sampled
with the local quadrature's normalized f0 sigma(g)g weights. This converges to
the stated event measure as hazard panels and partner quadrature are refined;
at finite order it is an additional approximation. After the event both products
travel ballistically to their first boundary exit. Remaining frozen-background
exposure on those paths measures the interactions being omitted. It is not
used to apply another collision or tune the correction.

At continuum resolution this mutually thinned closure is nonnegative: along
each original trajectory the combined symmetric sink is no larger than its
first-event density nu0 exp(-H); the two gain populations are nonnegative.
It does not claim to sample a collision history of a mutually interacting
particle system. Finite representative-event samples give signed estimators,
which can fail positivity or moment realizability even when the continuum
closure is positive. Such a failure remains visible in saved raw moments.
The reported exp(-tau) and 1-exp(-tau) are the frozen first-event proposal's
uncollided/event weights. Mutual partner thinning reduces actual events, so
these are not the final closure's ballistic-survivor fractions.
`ensemble_event_rate_over_emission` reports sum(event weights)/N; the associated
`ensemble_uncollided_source_number_fraction` is 1-2 sum(event weights)/N,
including both particles removed by every event. Samples are divisible by
the number of independent batches so boundary and moment estimators share
the same weights.

Point evaluation with collisions explicitly disabled calls the existing Cai
evaluator and reproduces its density, velocity and scalar temperature arrays
exactly. Enabled evaluation requires stated cylindrical bin widths; it is a
ring-volume estimator rather than a pointwise collision density. Comparisons
must use the same bins for the collisionless prediction and DSMC reference.
Bin volume quadrature is uniform in x and r^2, with Jacobian
dV=2 pi r dr dx. Tensor contributions are rotated to the local (x,r,phi) basis
before aggregation. The disabled point API returns the Cartesian basis, matching
the existing Cai API. Basis metadata distinguishes them. Pointwise tensor
results extend the existing disk integral with v_i v_j; they do not assert that
DSMC stores a comparable tensor.

## Collision law and reference provenance

The supported collision law is OpenFOAM v2512 `VariableHardSphere`, for a single
monatomic species. Its installed implementation was audited at
`src/lagrangian/DSMC/submodels/BinaryCollisionModel/VariableHardSphere/VariableHardSphere.C`.
The audited source SHA256 is
`d43f948867a526f0452a2aabe5ffdb610a2db318ae427c600b83dc1429246731`.
For identical species with reduced mass mu_r=m/2 the kernel is

    K(g) = pi d_ref^2 / Gamma(5/2-omega)
           * [2 k_solver T_ref/(mu_r g^2)]^(omega-1/2) * g.

Thus K has units m^3/s and nu0 has units 1/s. For repository argon
(m=6.63e-26 kg, d_ref=4.17e-10 m, omega=0.74, T_ref=273 K),
K(500 m/s)=2.8978791142e-16 m^3/s with the installed solver constant
k_solver=1.38065e-23 J/K, verified in its `etc/controlDict`. Cai source and
temperature moments retain the repository SI value k_B=1.380649e-23 J/K.
The constants differ by 7.24e-7 relative; neither legacy formulas nor solver
collision parameters are silently changed. With exact SI k_B in the kernel,
the independently checked reference is 2.897878610477046e-16 m^3/s.
The gamma normalization comes from the solver implementation; it is not inferred
from equilibrium mean free path. The solver's isotropic scattering was also
checked. In center-of-mass coordinates V=(v+w)/2 and g=v-w, the new relative
vector has magnitude |g| and uniform solid angle; products V+-g'/2 conserve
pair momentum and kinetic energy. Unsupported law names are rejected explicitly.
The zero kernel exists for numerical regression of the vanishing-collision
limit; it is not a fitted collision-strength coefficient.

The requested exact Hou 2025 reference was **not found** in the repository audit.
The untracked root `2927.pdf` is *Moon Rover Thermal and Power Analysis* by Creel
(LEAG 2023), not the requested Hou paper. No expression from that document or an
unverified Hou reference is imported. The model is Hou-inspired only in the
broad representative-event idea: the source, conditional partner law, symmetric
self-collision thinning and paired source/sink accounting are derived here.
There is no ambient/coflow reservoir and no substitution of plume density into
an ambient formula. A comparison with Hou's exact formulation remains pending
the exact bibliographic reference.

## Numerical assessment protocol

The scientific hypothesis is that this reduced correction improves selected
free-plume observables within a measured regime of validity, at a measured cost.
No improvement threshold is an acceptance test. Report baseline and corrected
DSMC-reference errors, their difference and signed bias; retain degradations.
Compute relative improvement only when baseline error is distinguishable from
zero. For density use volume-weighted normalized L2 error

    E_n^2 = sum V_i (n_model-n_DSMC)^2 / sum V_i n_DSMC^2,

on one common geometric/reference mask for both predictions. A failed corrected
bin is reported rather than silently removed to improve the corrected score.
Use global or fixed physical scales for components crossing zero. The existing
`masked_metrics` denominator and legacy output semantics remain unchanged.
Angular particle flux retains the repository's one-way exit-flux normalization;
cylindrical radial velocity and spherical outward velocity are distinct.
Common contour levels and any missing contour crossings must be reported.

The requested original six Kn values are {100,10,1,0.1,0.001,0.0001}; 0.01 is an
additional available Cai-family comparison. Discovery aggregates families and
does not alter extension generation/run exclusions. Inventory recomputes Kn
using each case's stated mean-free-path convention and nozzle diameter;
Kn alone does not establish similarity or validity. Source speed ratio, gas law,
geometry and source distribution must also agree. Interface-source cases change
the boundary-value problem and are excluded from the primary nozzle-source
assessment. Low-Kn direct-exit non-production labels remain in force; actual
resolution, occupancy and stationarity diagnostics control quality labels.
Exploratory/unresolved references do not support headline accuracy claims.

Convergence requires independently increasing sample count, hazard resolution,
partner quadrature, velocity cutoff, bin quadrature and domain extent on
representative observables. Assess signed moment errors as well as density:
a scalar-temperature match cannot establish tensor agreement. Model batch
uncertainty and DSMC uncertainty are separate. Averaged DSMC fields alone do not
provide a confidence interval; occupancy is not a substitute for time-block
estimates accounting for autocorrelation. Missing quantities and costs remain
unavailable. Numerical validity flags test density, finite values, tensor
symmetry and eigenvalues with an explicit relative roundoff tolerance of 1e-10;
they do not certify approximation accuracy.

Preprocessing and evaluation costs are measured separately, including requested
outputs, numerical resolution, cache state, wall time and process CPU time.
Peak memory is reported only where measurable. A fair Cai/model timing uses
identical bins and outputs. DSMC postprocessing time is not production solver
cost, and no production cost comparison is inferred without credible logs.
Local demonstrations are bounded and use at most eight cores; full assessment
and PBS commands read existing cases and never generate, run or submit DSMC.
There are no fitted physical parameters, training split or LLM dependency on HPC.

## Reproduction and bounded verification

From the repository root in the configured Python environment:

```bash
python -m pytest tests/unit/test_virtualcollision_model.py tests/unit/test_virtualcollision_kernel.py tests/unit/test_virtualcollision_assessment.py -q
python cases/cai2012-virtual-collision/assess.py inventory
python cases/cai2012-virtual-collision/assess.py smoke --output cases/cai2012-virtual-collision/results
python cases/cai2012-virtual-collision/assess.py full --output cases/cai2012-virtual-collision/results-full
```

The model tests independently check exact disabled Cai arrays, zero-kernel and
dilute limits, an eight-corner anisotropic population's known covariance,
between-cell velocity dispersion, invalid raw-moment preservation, the S=0
effusive source's known moments, phase-space quadrature convergence, two-particle
event bookkeeping, three boundary types and a known ray residence time. The
off-axis S=0 quadrature check found scalar-temperature errors about 1.22e-5,
1.6e-6 and 2.2e-7 of T_e at radial orders 10, 20 and 40 (other orders fixed
high); the regression requires the resolved result to improve and agree within
1e-5. This local check does not establish convergence of full plume observables.

`assessment.yaml` records the numerical controls and source-family search.
Machine-readable inventory, skip reasons, costs, convergence and raw moments
are written separately from display masks. Full-field comparisons additionally
write estimator metadata, DSMC-reference metrics, contour-location differences,
angular tables and plots when the needed inputs are available. Compact exported
profiles without cell centers, volumes and unaggregated raw moments cannot
support matched-estimator three-way accuracy scores; these skips are explicit.
Model tensors remain available even when DSMC tensor validation is unavailable.

For a manually initiated PBS assessment, review the template's environment and
resource directives, then run:

```bash
cd cases/cai2012-virtual-collision
qsub runAssess.sh
```

No scheduler job is submitted by the assessment driver. No missing cases are
generated and no production DSMC is launched. The bounded-run report and
inventory should be read together: availability of a saved profile does not
establish that its reference is resolved or that the model is effective.

## Initial local outcome

This section records the initial smoke assessment. The later export workflow
verified Cartesian geometry reconstruction and produced four-case comparisons;
see the [current usefulness review and refinement plan](../cases/cai2012-virtual-collision/refinement-plan.md).
Missing original C/V alone no longer blocks those comparisons. The review also
identifies limited saved regional density time series, distinct from the still
missing per-profile-bin velocity/temperature/tensor uncertainty histories.

The implementation was developed against source revision
`ad4cdc636b094b3bf1e0b601e7967b8b927370cc`; the assessment report separately
records its execution revision and working-tree state. Existing unrelated
workspace changes and case outputs were preserved.

Compatible full-field DSMC inputs for a matched raw-moment comparison were not
available locally. Saved extension profiles and diagnostics were found, rather
than assuming the old README's run-availability statement remained current.
Those compact profiles cannot recover cell geometry or unaggregated moments;
their absence is a specific assessment-input limitation, not evidence for or
against the model. Low-Kn non-production references and downstream interface
sources remain separate from primary direct-exit evidence.

The bounded Kn=10 model demonstration produced realizable moments in all three
requested bins. That is an internal physical check. Its refinement study did
not establish numerical convergence, and no improvement or degradation against
DSMC has been assessed. The measured timing, conservation residuals, exposure,
remaining exposure and refinement differences are retained in the bounded-run
report; they must accompany any discussion of cost or validity. The thesis
hypothesis therefore remains open.

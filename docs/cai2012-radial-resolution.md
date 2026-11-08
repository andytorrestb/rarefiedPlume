# Cai exit radial DSMC resolution study

This study determines local cell and parcel requirements in a narrow expanding
streamtube. The physical reference is the rocket exit `Kn_e = lambda_e / D_e`:
argon at `100, 10, 1, 0.1, 0.01, 0.001, 0.0001`. It inherits nozzle, exit and gas
inputs from the [Kn-extension reference](../cases/cai2012-kn-extension/baseCase/case.yaml),
with snapshots embedded in each child. Legacy cases and their numerical policies
are preserved. Follow the [tutorial](tutorials/cai2012-radial-resolution.md).

## Physics and fixed inlet

The existing `cai2012ext.interface.equivalent_source` supplies the spherical
isentropic expansion and Bird parameter. No new expansion physics is introduced.
The axis anchor fixes the default inlet at `r_in / D_e = 0.7172871`. Its mapped
temperature is 299.845 K, while the physical Cai exit is 300 K; the small
difference comes from the existing gross outward-flux mapping. Inlet drift is
uniform and axial. This radial surrogate does not reproduce the finite exit disk
or angular variation of the full plume.

The default outer boundary is `25 D_e`. At general Bird threshold 0.05, predicted
breakdown is at approximately `2.788 D_e` for Kn 0.001 and `14.992 D_e` for
Kn 0.0001. Kn 0.01 and larger are already noncontinuum at the anchor. Results
distinguish `CROSSING_IN_DOMAIN`, `ALREADY_NONCONTINUUM_AT_INLET`,
`NOT_REACHED_IN_DOMAIN`, and `INVALID_STATE`; the outlet is never substituted for
a missing breakdown crossing. `radial_approx` is an explicitly named alternative.
A sphere at the axis radius is a surrogate surface, not a validated plume surface.

The square frustum has `y,z = +/- r tan(alpha)`, specular sides and vacuum outlet.
The narrow angle is selected once from the initial predicted target and capped
at 0.5 degrees, unless explicitly overridden. Inlet, outlet and numeric angle
remain fixed across numerical revisions. An inlet override must retain an
available upstream continuum interval. Cartesian parcel weighting is used.
Loading a case restores its recorded source coefficients, inlet and physical
sampling times. Small numerical differences between Windows and WSL do not
change that case's identity; an edited case record fails its checksum.

## Predict, measure, correct

The initial mesh targets `min(1.0 lambda_predicted(r), 0.05 D_e)`, including
slanted and transverse edges. These are provisional model estimates. The final
goal is **largest actual edge <= 0.8 measured local lambda**, evaluated using
the same existing VHS convention on measured positive density and temperature.
The initial prediction is not a measurement of the initially empty domain.
There is no automatic 0.8, 0.4, 0.2 mesh sweep.

Mean occupancy must be at least 20 parcels/cell. Initial global weight is
`W = min(n_estimate_i V_i) / P`. A mesh correction uses parent measured density
on smaller child volumes to recompute W at the same requested P. A population
correction preserves the mesh and reduces W using measured occupancy. Empty or
invalid cells require sampling work; no inlet-only replacement is fabricated.
The streamtube's global weight and total parcel count are not full-plume costs.

Timestep checks are independent: defaults are `deltaT/tau <= 0.2` and a conservative
transport ratio `deltaT (|U| + 3 sqrt(R T)) / minimum_cell_length <= 0.2`, with
`tau = lambda / c_bar`. The runner fixes deltaT, rounded down to a whole number
of steps per physical write window. Explicit timestep-only revisions preserve
mesh and weight. Resource gaps are reported; the runner requires an explicit
`--allow-budget-exceedance` to run beyond declared budgets.

## Sampling and evidence

Transit time integrates `dr/U_model(r)` over the fixed domain. Defaults are two
transits of startup, six of sampling, six equal blocks and two disjoint native
write windows per block. Startup is rounded up to a write boundary and the
requested/achieved durations are recorded. Physical durations remain unchanged
under numerical revisions. The native field model keeps `resetAtOutput on`
throughout the run. It never switches later fields to cumulative averages.

This port writes `rhoN_Ar`, `Ttra_Ar`, `U_Ar`, and `dsmcNMean_Ar`; `UMean_Ar` is an
accepted velocity alias. Native `mfp_Ar`/`mct_Ar` are cross-checks. Empty-cell
sentinels cannot count as physical measurements. Window records contain times,
steps, field hashes and mesh/physical identity. Cases have reproducible distinct
hash-derived seeds. Execution is serial; MPI streams are not exposed.

Post-processing combines density/momentum/energy moments, then derives velocity
and temperature. It retains native radial layers, six block profiles and block
standard errors. Drift and strong lag-one correlation are flagged; these
diagnostics do not establish statistical independence. Invalid coverage or
incomplete/unsettled blocks leaves acceptance unresolved. The older rarefaction
diagnostic's strict `>20` mask is labeled separately from this study's inclusive
`>=20` target.

Mesh comparisons hold parcel target fixed; parcel comparisons hold mesh fixed.
Conservative volume-overlap projection of moments onto the coarser intervals
preserves resolution and thermal energy. Optional L2 tolerances are unset by
default. When declared, uncertainty/coverage can leave them inconclusive.
Missing extra comparison cases do not fail an otherwise satisfied sizing case.
Departure from the isentropic reference beyond breakdown is not a mesh failure.

CSV, JSON, Markdown and optional plots report local MFP, actual edges, mean
occupancy, timestep ratios, breakdown, sampling and remaining corrections. The
portable plan is bound to the parent fields; inherited child estimates require a
new run. See [outputs](reference/outputs.md).

## Verification limits

The bounded native mechanics fixture uses Kn 100, outer radius 0.82D, angle
0.5 degrees and geometric cap 0.04D: 3 cells, 656 steps, about 92 predicted
parcels, 12 post-startup windows. It ran with the owner's OpenFOAM v2512
hyStrath port and wrote clouds, direct fields and complete checkpoints.
Post-processing identified one cell just below occupancy 20 without requesting
mesh refinement. This checks interoperability, not physical validation.
The seven full scientific cases have not been run.
The [verification record](../cases/cai2012-radial-resolution/verification.json)
records executable/mesh hashes and observed mechanics results. Local generated
artifacts retain the actual meshes, fields and logs.

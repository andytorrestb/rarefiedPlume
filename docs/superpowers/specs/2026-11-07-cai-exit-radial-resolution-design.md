# Cai exit Knudsen number: radial DSMC resolution study

Date: 2026-11-07 (America/Chicago). Status: approved for implementation planning,
including the owner's requirement that the mesh refinement utility be portable
to other cases. The implementation plan is reviewed before execution.
No new solver or physical expansion law is needed.

## Intent

Provide a simple, repeatable study for choosing DSMC cell sizes and parcel
populations as a function of rocket exit Knudsen number. A user configures a
physical case, runs an initial prediction-based mesh, measures local mean free
path and occupancy, and creates a revised case only where the measured numerical
requirements are missed. Optional comparisons show changes in the computed flow.

Use the expanding streamtube in `cases/1D` as the geometric template and the
organization of `cases/cai2012-kn-extension` as the workflow template. Preserve
those existing cases. Run the new cases with the owner's OpenFOAM v2512 port of
hyStrath's `dsmcFoam+` and `dsmcInitialise+`.

This is a radial surrogate of the Cai plume. Its results establish resolution
requirements for that surrogate and supply local mesh/occupancy guidance for a
later Cai comparison. A total parcel count or a global particle weight for the
streamtube is not a resource estimate for the full three-dimensional Cai domain.

## Agreed decisions

- The physical reference is the rocket exit plane: `Kn_e = lambda_e / D_e`.
- Use argon at `Kn_e = [100, 10, 1, 0.1, 0.01, 0.001, 0.0001]`, with the
  physical exit inputs used by `cases/cai2012-kn-extension`.
- Use the existing equivalent spherical-source mapping to connect the exit to
  a narrow expanding streamtube.
- Choose a fixed inlet for each physical Kn case. Keep its location and
  prescribed state fixed during all numerical refinements of that case.
- Where the mapped flow has an upstream continuum interval, start DSMC inside
  it and calculate the downstream breakdown location.
- Refine mesh and parcel population independently; retain previous runs.
- Check local cell-edge size against mean free path and use 20 parcels per cell
  as the baseline occupancy requirement.
- Treat `Delta_edge_max <= 0.8 lambda_measured` as the final mesh requirement.
  Predicted MFP sizes the initial case; measured MFP determines whether further
  refinement is needed. There is no automatic `0.8 -> 0.4 -> 0.2` sweep.
- Verify and reuse the existing Stewart/isentropic/breakdown implementation.
- Make mesh assessment and refinement reusable by other cases. Keep Cai exit
  physics and radial geometry outside the shared refinement utility.

## Approaches considered

1. **Fixed-inlet radial study (selected).** Derive a physical case once and hold
   its inlet fixed while correcting measured mesh and occupancy gaps. This
   retains direct comparisons and includes the cost of resolving the dense region.
2. **Breakdown-placed inlet.** Move low-Kn injection downstream while retaining
   continuum overlap. This can reduce cost but adds inlet-location sensitivity
   to the study. It was not selected.

## Scope and boundaries

Deliver a shared utility `plumetools.refinement`, a new family at
`cases/cai2012-radial-resolution`, its supporting package
`plumetools.radialstudy`, and guides for both. The family has one editable
`study.yaml`, a clean hyStrath case template, case generation/refinement,
per-case mesh/run/post commands, and study analysis.

In scope: pre-run sizing; static graded meshes; global parcel weights; the
seven-Kn matrix; measured occupancy/MFP/timestep diagnostics; radial profiles;
independent refinement comparisons; missing-evidence reporting.

Out of scope: a new C++ solver, changes to the WSL hyStrath checkout, CFD
coupling, moving the inlet during refinement, AMR, variable radial parcel
weighting, full Cai plume validation, scheduler submission, a GUI, and a cost
fit for an unrun three-dimensional plume. MPI execution may use the same
declared case settings; the first working path and mechanics check are serial.

Failure conditions include silently relaxing a cell criterion to fit a budget,
changing physical inputs between numerical variants, replacing existing runs,
counting accumulated parcel samples as occupancy, or presenting agreement with
an isentropic curve outside its validity as a numerical-convergence requirement.

## Portable refinement utility

The reusable utility owns measured-MFP assessment, target cell-edge limits,
occupancy checks and weight estimates, refinement-plan provenance, and target
export. It must work without an exit Kn, nozzle diameter, spherical source,
radial coordinate, Cai configuration, or this repository's case directory layout.
Its pure Python API accepts cell geometry, normalized measured arrays, gas
properties/MFP convention, and explicit target/sampling settings.

Provide `python -m plumetools.refinement` commands for assessment, plan export,
and explicit mesh materialization. OpenFOAM readers support reconstructed static
ASCII cases using standard `dsmcFoam` averaged moments, hyStrath's direct fields,
and configurable field names. Gas properties, species labels, sampling-window
metadata, and output directories are configurable. Reading another case never
requires editing a Cai or radial-study configuration.

Reuse the existing ASCII mesh/field parsers and VHS physics. Unsupported binary,
decomposed-only, dynamic, multispecies-mixture, or otherwise incompatible input
gets a specific actionable explanation and leaves the case unchanged. Other
solvers can supply normalized arrays through the Python API. Extending readers
does not change the refinement decision rules.

Separate target selection from mesh construction through a small backend
interface. The first implementation provides:

1. A radial-study backend that regenerates its graded streamtube efficiently
   from measured target edges while retaining its fixed physical domain.
2. A standalone OpenFOAM v2512 backend for supported three-dimensional hex
   meshes, using existing `refineHexMesh` cell-set refinement and its conformity
   extension. Only failing cells and cells required for conformity are split.
   This backend's local factor-of-two splits can produce edges below the 0.8
   target; that is a mesh-tool constraint, not an additional target or sweep.

The generic backend writes a new mesh artifact, source/target records, and
inherited planning fields in a fresh output directory. It does not copy or map
the parent's DSMC parcel cloud, fabricate new measured fields, initialize an
unknown solver, or run a simulation. A case family installs that mesh in a new
case through its own initializer. Never modify the source mesh or results.
Mapped target/density fields used during refinement remain labeled as estimates
inherited from the parent; they cannot prove the child flow meets the criteria.

If a topology cannot be refined by the supplied backend, assessment and target
export remain useful and clearly report that a different mesh backend is
needed. Export target edges and source cell identities/geometry so another
mesher can consume them. Verify changed-mesh detection before applying any plan.

Demonstrate portability on an unrelated Cartesian OpenFOAM case using custom
field names, outside a study tree and without any Cai metadata, including a
real local mesh refinement. Keep scope to this concrete API and the two initial
backends rather than a general plugin-discovery framework.

## Existing implementation and verified baseline

| Need | Existing source and use |
|---|---|
| Exit gas, density from Kn, speed from speed ratio | `plumetools.cai2012.gas` and `inflow.from_config` |
| Exit-to-source mapping and finite-temperature radial state | `plumetools.cai2012ext.interface.equivalent_source`, `EquivalentSource.state` |
| Stewart spherical relations and hard-sphere regression | `interface.stewart_spherical`, `stewart_hard_sphere_density`; reuse, no copied equations |
| Bird parameter, general and radial approximation | `interface.bird_parameter` / `EquivalentSource.bird` |
| Surface search | Inspect/reuse `interface.select_sphere` where its interface-specific contract fits; adapt its result semantics for a radial breakdown diagnostic |
| Mean free path and common diagnostics | `plumetools.rarefaction.core`, using a hyStrath input adapter |
| OpenFOAM field/mesh reading and dictionary formatting | `plumetools.foamio`, `plumetools.mesh`, `foamio.primitives` |
| Non-destructive generation | Follow `cai2012ext.generation` staging, identity, and conflict behavior; its Cai-specific policy and configuration assumptions are not the new family's schema |
| Difference metric and uncertainty-aware verdict | `cai2012ext.overlap.normalised_l2` and `verdict`, where applicable |

The existing interface test module was run during reconnaissance: **66 passed**.
That verifies the existing numerical checks; it is not a successful simulation
or validation of this new study.

The WSL Ubuntu checkout `/home/andy/hyStrath-port` exists. After sourcing
`/usr/lib/openfoam/openfoam2512/etc/bashrc`, both solver help commands start and
report OpenFOAM 2512. Sources for the required inflow, deletion, and measurement
models exist. No case was run as part of this design work.

## Physical inputs and reference-state mapping

The default physical reference is
`../cai2012-kn-extension/baseCase/case.yaml`, relative to the new family. Read
its `nozzle`, `exit`, and `gas` inputs through existing Cai code; do not inherit
its three-dimensional mesh, particle targets, run scripts, or historical Kn
exclusion policy. Snapshot the resolved physical inputs and reference-file hash
into each generated case so it remains interpretable after the reference changes.
A user can point to a compatible reference copy to change physical settings.

Defaults from that reference:

| Input | Value |
|---|---|
| Exit diameter | `D_e = 0.2 m` |
| Exit temperature | `T_e = 300 K` |
| Exit speed ratio | `S_e = U_e / sqrt(2 R T_e) = 2` |
| Species | Ar; molecular mass `6.63e-26 kg` |
| VHS diameter, exponent, reference temperature | `4.17e-10 m`, `omega = 0.74`, `273 K` |
| Internal degrees of freedom | 0; hence `gamma = 5/3` |
| MFP convention | Existing temperature-dependent VHS convention |

Only density changes across the default Kn sweep: derive `lambda_e = Kn_e D_e`
and `n_e` using the existing conversion. Do not reinterpret these as sonic-radius
Knudsen numbers or use Stewart's hard-sphere inversion for the production VHS
inputs.

Use the **axis** of the existing equivalent source (`theta = 0`) as the radial
reference. Its anchor `r_a` fixes the default inlet. The current mapping gives
`r_a / D_e = 0.7172871`, `U_in = U_e = 706.9529 m/s`, and
`n_in = n_e`. Its mapped temperature is `T_in = 299.845 K`, slightly different
from the physical exit's 300 K because the mapping conserves the gross outward
exit flux targets. Record both states and this reason; do not silently replace
the mapped temperature or retune the mapping.

The new study uses the radial axis state to prescribe a uniform inlet
Maxwellian, directed along the streamtube axis. It does not reproduce the finite
disk or inject the full Markelov angular profile. The narrow-cone approximation
and the existing exit-to-source assumptions remain part of every case's record.

## Domain and inlet

Builder calls proposed for approval:

- Default nominal radial interval: `r_in = r_a`, `r_out = 25 D_e`. The larger
  outer radius includes the predicted breakdown of the lowest-Kn case. Matching
  Cai's exit inputs does not require retaining Cai's `10 D_e` outer boundary.
- Allow an explicit `r_in / D_e >= r_a / D_e` and an explicit outer radius.
  An override creates a different physical case identity. Its inlet must remain
  upstream of an available breakdown crossing if an upstream continuum interval
  exists at the anchor.
- Geometry is a square expanding streamtube, following the legacy template:
  axial coordinate `x = r`; side planes `y,z = +/- r tan(alpha)`; planar inlet
  and vacuum outlet. Nominal `r` and actual distance from the source are labeled
  separately where they differ. The four sides have specular symmetry behavior,
  using the patch/model combination supported by this hyStrath port.
- Select a narrow angle once per Kn case from the initial predicted mesh target:
  `tan(alpha) <= 0.25 min_r(s_0(r)/r)`, capped at `alpha = 0.5 degrees`, where
  `s_0` is defined below. This keeps baseline transverse edges smaller than the
  radial target. Allow an explicit angle override and report its actual cell
  requirements. Freeze the resolved angle before generating any refinements;
  a revised mesh never narrows an existing case's physical domain.
- Use Cartesian parcel weighting, as in the legacy streamtube. The port's
  `dsmcSpherical` model changes parcel weighting by cloning/deletion; source
  inspection did not establish it as a replacement for this streamtube's
  geometry or spherical particle-motion remapping. Do not enable it implicitly.

All refinements of a Kn case share the same domain, inlet state, gas, source
mapping, startup interval, sampling interval, and physical sampling cadence.

## Predicted MFP, measured MFP, and the mesh acceptance target

Keep two distinct local MFP profiles, evaluated through the same existing VHS
function:

- `lambda_predicted(r)`: calculated before the run from the prescribed inlet
  state and existing radial expansion model. This is the initial mesh estimate.
  It is not the MFP of the initially empty numerical domain.
- `lambda_measured(i)`: calculated after startup from adequately sampled local
  number density and translational temperature in each actual cell. This is the
  evidence used for mesh acceptance and a possible revised mesh. It is not the
  mean free path held constant at the inlet.

The final spatial requirement is

`Delta_edge_max(i) <= 0.8 lambda_measured(i)`.

Once that requirement is met throughout the assessed domain with adequate data,
the default workflow requests no additional mesh refinement. Incomplete or
under-sampled cells make the assessment unresolved; they cannot be silently
omitted to obtain a passing result. Occupancy and timestep requirements remain
separate checks.

Proposed initial mesh `r0` is a static graded radial mesh with

`s_0(r) = min(1.0 lambda_predicted(r), 0.05 D_e)`.

The initial factor 1.0 is a designer proposal and is configurable. It provides
a provisional starting mesh near the predicted mean free path; the final
acceptance ratio remains 0.8 against measured MFP. The geometric cap retains
useful radial resolution where the predicted MFP exceeds the whole domain.

Construct intervals so the **largest actual edge**, including slanted edges,
satisfies the requested predicted target over the interval. Check planned
vertices and volumes; axial spacing or `V^(1/3)` alone is insufficient. With an
explicit wider angle, allocate enough transverse cells for the same criterion,
or report that a declared resource limit prevents the requested mesh.

After measurement, generate a new static mesh only if the final ratio is
exceeded. Derive the required edge limits directly from `0.8 lambda_measured`
in the failing region. Retain existing radial boundaries where practical and
insert as few extra intervals as needed to satisfy those limits. Choose integer
transverse subdivisions from the same actual-edge requirement. The conformal
streamtube topology can require a transverse count increase outside a failing
radial layer; show that cost explicitly. Keep the inlet, cone, and outlet fixed.
Do not require global halving or introduce predefined finer target ratios.

Use a conservative radial representation of the valid measured cell bounds
(including the smallest MFP among transverse cells in an interval); interpolation
must not hide a tighter measured bound. Inadequate measured coverage leaves the
refinement plan incomplete and identifies the needed sampling/occupancy work.
This remains a one-dimensional static mesh generator, not arbitrary polyhedral
AMR. These geometry rules belong to the radial backend; they are not
assumptions in the portable assessment or target-plan API. The builder may
group intervals into graded blocks without relaxing the
requested edge limits. Report unavoidable rounding or extra resolution.

Report predicted and measured MFP side by side, current edge/MFP ratios, failing
cell counts/volume, and the proposed new cell counts and spacing. A resource
budget never causes automatic coarsening. A revised mesh is an estimate based
on its parent's measured flow; run and reassess it because the measured state
can change with the new mesh.

The earlier `lambda/3` to `lambda` band is context for the resolution criterion,
not a request for a mandatory refinement sweep. Smaller cells caused by the
geometric cap or conformal grading remain admissible and are reported. Users
may explicitly request stricter targets for additional studies; default
completion requires the chosen 0.8 measured-MFP target.

## Parcel sizing and independent variation

The default acceptance target is 20 mean parcels per cell. Use a single global
`nEquivalentParticles = W` for each case, without radial weighting. For the
initial planned mesh and target `P`, derive

`W = min_i[n_model(r_i) V_i] / P`.

This targets the lowest predicted occupancy over the domain, rather than only
the inlet cell. Report predicted minimum, distribution, total parcels, and the
cell controlling the weight. These are estimates. Measured mean occupancy is
the post-run evidence.

- A mesh refinement holds `P` fixed and recomputes `W` for the smaller volumes,
  using the parent's valid measured density as the updated population estimate.
- A parcel correction holds the mesh fixed and reduces `W` as needed to reach
  `P` in the least-populated adequately measured cell. For an unchanged mesh,
  estimate `W_new <= W_old min_i(P_measured(i) / P)`; remeasure after running.
- Higher parcel targets such as 40 or 80 remain explicit optional studies.
  Doubling a target halves its planned weight and approximately doubles the
  population. No automatic 20/40/80 matrix is required for completion.
- Record both the requested target and the resulting weight. Reusing a coarse
  mesh's weight unchanged after spatial refinement is not the default.
- Report the formula `P_i = n_i V_i / W` alongside recommendations, so a user
  can apply local results to different cell volumes. Streamtube weights and
  total populations depend on its represented volume.

Use measured `dsmcNMean_Ar` for mean parcels per cell. `dsmcN_Ar` can supplement
it with instantaneous counts. Do not use number density, the solver log's
domain total, or the number of accumulated samples as a per-cell population.
Account for the port's `0.001` empty-cell placeholder using positive measured
density/temperature and the available count fields.

The study's occupancy pass condition is `P_mean >= 20`. The existing generic
diagnostic validity mask uses a strict `> min_sampled_particles`; preserve that
library behavior and label its mask separately rather than silently treating
the two threshold definitions as identical.

Sparse or empty cells need an explicit sampling/population response and remain
unresolved; an undefined density or zero count cannot establish a new global
weight. Preserve the source field/time/window used for each revised estimate.

## Time step, sampling, and run behavior

Derive a fixed timestep for each numerical case from both a collision-time
limit and a transport limit. Proposed safety factors are 0.2 for each. Use the
existing gas/state functions for collision time and characteristic speed, and
the actual minimum cell dimensions for transport. Expose explicit timestep
overrides, recording achieved ratios. No automatic adaptive timestep, mesh
adaptation, load balancing, or parcel weighting is enabled in this first study.

Define a common transit time by integrating `dr / U_model(r)` over the fixed
domain. Proposed defaults: two transits of startup followed by six transits of
sampling, divided into six equal physical-time blocks. These durations and
block boundaries stay fixed across numerical refinements; round only solver
steps/cadence and record achieved times. Allow users to extend them. There is
no claim that the default duration establishes stationarity at every Kn.

Use disjoint measurement windows rather than treating successive cumulative
averages as independent blocks. Preserve at least one complete field snapshot
per block and enough subwindows to check drift. Use `seedNumber`, which was
found in the port, with stable distinct recorded seeds per case identity.
If parallel execution is enabled, verify distinct streams across ranks and
record the effective seed policy: this port does not automatically offset an
explicit `seedNumber` by rank. Reproducibility is conditional on the recorded
solver and MPI layout; a fixed seed does not establish independent samples or
equivalence across rank counts.

The per-case commands verify the plus solver and OpenFOAM environment, mesh,
initialize an empty domain, and run the prescribed inlet-to-vacuum expansion.
Fresh initialization refuses existing solver output. A resume preserves the
case inputs and sampling definition and uses available complete checkpoints;
it never silently starts over. The WSL paths found during reconnaissance are
documented setup examples, not hard-coded paths in portable runners.

## Breakdown diagnostic

Use the existing general Bird parameter with local VHS MFP by default,
`B = -(lambda U / c_bar) d(ln rho)/dr`, with threshold 0.05 configurable.
Also report the existing `radial_approx` form used in Stewart's Eq. (3) as a
named comparison. Do not mix those definitions or the hard-sphere convention.

Search from the fixed inlet through the configured outer radius for the first
outward crossing. A radial wrapper may adapt the existing search to the axis;
it must not duplicate the physical expansion or Bird equations.

Return distinct outcomes:

- `CROSSING_IN_DOMAIN`: inlet is below the threshold and an outward crossing
  is bracketed; report `r_B / D_e` and its state.
- `ALREADY_NONCONTINUUM_AT_INLET`: inlet already meets/exceeds the threshold;
  no upstream continuum interval is claimed.
- `NOT_REACHED_IN_DOMAIN`: no crossing within the configured domain. A
  domain-limited result from the existing interface selector does not make the
  outlet a breakdown surface.
- `INVALID_STATE`: missing/non-finite/nonpositive inputs or invalid model state;
  report the reason rather than constructing a radius.

For the default mapping, general Bird threshold 0.05 gives approximately
`r_B / D_e = 2.788` at Kn `0.001` and `14.992` at Kn `0.0001`. Kn `0.01` has
`B_in = 0.06243`, and the larger Kn cases are also already above the threshold
at this anchor. Do not force all seven cases to have a continuum region.

This location belongs to the radial-axis surrogate. A sphere at that radius
is the surface implied by this radial approximation, not a validated full
rocket-plume breakdown surface. No angular-core choice from the existing
three-dimensional interface study is smuggled into the radial result.

Overlay predicted breakdown on radial profiles. Where adequately sampled
results permit, also report the measured radial Bird profile and crossing,
reusing gradient helpers and providing validity/coverage. Missing or noisy
gradients leave that measured result unavailable; the prediction remains labeled.
Departure from the isentropic temperature downstream of predicted breakdown
is physical behavior to inspect, not automatically a failed refinement.

## User workflow and case identity

The proposed controls in `study.yaml` are the reference case, seven Kn values,
inlet/outer radius settings, cone-angle setting, predicted/measured mesh targets,
parcel targets, timestep policy, startup/sampling/block settings, seeds, and optional
resource/convergence limits. Reject unknown keys and invalid physical values.

Provide these family entry points (exact option spelling may be settled by the
builder; the behavior is required):

1. `generate_cases.py --dry-run`: report settings, source mapping, breakdown,
   planned cells/weights/populations/steps, and budget gaps. Run no solver.
2. `generate_cases.py`: generate only the seven prediction-based `r0_ppc20`
   initial cases by default. Filters select Kn; explicit overrides can request
   different initial mesh or parcel targets. There is no predefined finer-mesh
   matrix independent of measured need.
3. Each case has `Allmesh`, `Allrun`, and `Allpost`. A user can run a single
   case without running the rest of the study.
4. `refine_case.py <case> --dry-run` reads a complete measured assessment and
   reports needed mesh, parcel-weight, or timestep corrections. `--mesh` or
   `--parcels` can isolate a correction; their defaults hold the same target
   occupancy or the same mesh respectively. Applying a plan writes a new
   numerical revision, without changing the parent or launching a run. If all
   chosen criteria pass, the default command reports that refinement is
   unnecessary and creates no duplicate case. Missing measurements are reported
   as missing, rather than replaced by an unannounced prediction-based plan.
5. `analyse.py` inventories available cases and results, compares valid pairs,
   and reports remaining measured gaps. It does not require further runs after
   the chosen criteria are satisfied.

Use stable directories such as `Cases/Kn0p001/r0_ppc20` and a new `r1_ppc20`
after a needed correction. Revision numbers do not imply fixed edge ratios or
cell-count multipliers; independent mesh hashes distinguish mesh changes from
weight-only revisions. A case record includes
resolved physical/source/domain inputs, geometry identity, numerical settings,
parent refinement identity, measured source fields/windows and hashes where
applicable, generator version, and solver provenance. Study and case records
must survive moving the family tree.

Generation is idempotent. Matching cases are kept byte-for-byte; physical or
numerical conflicts are refused with a concrete difference report. Existing
meshes, checkpoints, and results are never overwritten by regeneration or
refinement. New physical settings require a new study/output tree. No implicit
clean command or deletion is part of the normal workflow.

## Measurements, comparisons, and outputs

The hyStrath adapter reads `rhoN_Ar`, `Ttra_Ar`, `UMean_Ar`, and
`dsmcNMean_Ar`, plus actual cell geometry and sampling windows. Supply direct
temperature to `rarefaction.core.Context`; a field-name substitution alone
would incorrectly retain the stock solver's temperature-from-moments dependency.
Enable the port's native MFP and collision-time measurements as cross-checks.
Canonical study MFP uses the existing Cai-compatible VHS convention on measured
positive density and temperature; native sentinels are not valid physical values.

Per case, write machine-readable settings/plan, a diagnostic summary and radial
CSV, and plots of:

- `n / n_e`, `U_x / U_e`, and `T_tra / T_e` against nominal `r / D_e`, with
  uncertainty where available and the continuum curve labeled as a reference;
- predicted and measured local mean free path;
- actual maximum-edge/each-MFP ratio, showing the measured acceptance target 0.8
  and the earlier 1 and 1/3 context;
- mean parcels per cell, showing 20;
- timestep/collision-time and transport ratios, with configured limits.

Preserve spatial resolution when forming profiles: include each mesh's native
radial layers. For a comparison, conservatively project the finer result onto
the coarser radial intervals using cell-volume overlap and average transverse
cells consistently. Do not let a fixed coarse plotting grid hide refinement
differences. Record comparable-domain and valid-data coverage.

Compare only cases with the same physical/domain/sampling identity. Show mesh
pairs at fixed parcel target and parcel pairs on a fixed mesh. Report relative
L2 differences in density, axial velocity, and temperature using the finer or
higher-parcel result as reference. Reuse the existing normalized metric and
uncertainty-aware verdict, with radial-bin volume weights. Six disjoint block
profiles supply standard errors; test temporal drift and flag when block
independence/stationarity is unsupported. Such cases remain inconclusive.

The normal completion status is **criteria satisfied** when reliable measured
coverage demonstrates the selected 0.8 mesh ratio, at least 20 mean parcels per
cell, and the timestep limits throughout the assessed domain. Insufficient data
leaves the relevant criterion unresolved. A passing initial or revised case
needs no additional mesh sweep or pairwise profile comparison to meet this
study's numerical requirements.

Optional profile-convergence tolerances are **unset by default**: when multiple
revisions exist, the report supplies differences, uncertainty, and measured
resolution/occupancy gaps for the user to assess.
If the user declares tolerances, report pass/fail/inconclusive per pair using
the existing metric contract. Missing blocks, inadequate occupancy, invalid
coverage, unsettled sampling, or absent comparison cases prevent a positive
optional convergence verdict, but absence of extra comparison cases is not a
failed mesh criterion. Report numerical criteria separately from empirical
profile convergence; neither establishes accuracy against the full Cai plume.

The study summary contains, for each generated case: Kn, inlet/outer radii,
revision and mesh identity, target occupancy, cell count, weight,
predicted/measured population,
worst and distributional mesh/occupancy/time ratios, violation counts and volume
fractions, breakdown outcome, sampling status, and available refinement pairs.
Include CSV plus a readable Markdown report and overlay plots. Identify the
needed correction, or explicitly state that the measured criteria are met,
without automatically launching another run.

## Verification and delivery

Implementation must verify physical input inheritance at all seven Kn values;
agreement with the existing mapping/Bird functions; immutable physical identity
under refinement; use of measured rather than inlet-only MFP in the acceptance
decision; planned versus built mesh geometry/volume; direct satisfaction of
measured edge bounds by a revised mesh; correct global-weight scaling and
predicted minimum occupancy; native hyStrath field/count/sampling interpretation; and
idempotent generation that preserves a case containing results.

Use bounded synthetic fields to check measured violations, missing coverage,
comparison on unequal radial meshes, distinct breakdown outcomes, and honest
inconclusive results when evidence is missing. Check that an already compliant
case generates no further revision, and that a local failure does not trigger a
predefined global bisection sequence. Reuse existing physics tests
rather than copying implementations into new tests.

After design approval and an implementation plan, run a tiny explicitly labeled
mechanics case in WSL through mesh, initialization, a short `dsmcFoam+` run,
field reading, and analysis. Confirm boundary/model compatibility with the port
and record the measured outcome. This short case is not a convergence result.
Also verify generation of all seven initial cases and the independent
correction paths using measured fields. Do not launch the scientific matrix or
submit scheduler jobs as part of the implementation checks.

Verify the shared utility separately on the unrelated Cartesian case, both
with field-name overrides and with its standalone hex refinement backend.

Deliver working documented commands, the new family and library support,
focused checks, and the actual verification record. Existing case inputs and
hyStrath sources remain intact.

## Calls left to the builder

Module organization, dictionary assembly details, radial-block grouping,
equivalent CLI spelling, plotting layout, and manifest serialization details
are implementation choices. Use existing repository conventions and avoid
duplicated physical formulas. The shared utility's independence from this
study and its concrete second-case verification are requirements.

The default domain extent, automatic narrow-cone rule, initial predicted-MFP
factor and geometric cap, duration/block defaults, and absence of default
convergence tolerances above are proposals made by the designer, not choices previously
made by the owner. They must be visible in the approval handoff.

## Sources

- Owner's design discussion in this chat, including selection of approach 1.
- `cases/1D`; `cases/cai2012-kn-extension` and its physical reference template.
- `plumetools/cai2012ext/interface.py`; `docs/cai2012-kn-interface.md`.
- Stewart and Lumpkin (2012), *Axisymmetric Plume Simulations with NASA's DSMC
  Analysis Code*, Eqs. (1)-(4), pp. 2 and 4-6.
  Owner-supplied file: `C:/Users/andy9/Downloads/Documents/Stewart2012.pdf`.
  Public copy: <https://ntrs.nasa.gov/citations/20120016934>.
- Read-only inspection of the owner's WSL hyStrath-port sources, including
  `dsmcVolFields`, `dsmcSpherical`, `dsmcCloud` (`seedNumber`), and runtime help.
- Read-only inspection of installed OpenFOAM v2512 `refineMesh` and
  `refineHexMesh` sources and help, including `refineHexMesh` field mapping and
  default extension of selected cells to maintain the 2:1 condition.

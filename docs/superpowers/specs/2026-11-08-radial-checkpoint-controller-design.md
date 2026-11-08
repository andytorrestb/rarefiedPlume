# Radial study: coarse starts and automatic checkpoint corrections

Date: 2026-11-08. Status: approved by the owner, including the proposed defaults.

## Intent and agreed behavior

Make the existing serial WSL study useful for quick numerical sizing. Start each
of the seven argon exit-Kn cases cheaply, measure what its cells require, and
spend additional computation only on useful sampling or measured corrections.
Keep the existing exit physics, mapped upstream inlet, and breakdown calculation.

The owner requested:

- A fixed 2.5-degree half-angle and outer distance of 10 nozzle-exit diameters.
- Approximately 100 total cells for the initial mesh, then measured refinement.
- Periodic checks of the numerical requirements in every cell.
- On a trusted failure: stop, clean outputs, adjust settings, and restart from
  initial conditions. Do not carry the previous parcel cloud into a changed case.
- A dedicated, concise results log that persists across all retries.
- The checkpoint controller approach around the installed native dsmcFoam+
  solver. Continue to run one case at a time in WSL.

Earlier requirements still apply: measured maximum edge <= 0.8 local mean free
path, mean occupancy >= 20 parcels per cell, separate collision/transport
timestep limits, no further refinement of passing cells simply for a sweep,
and a refinement utility reusable by other cases.

## Approaches

1. **External checkpoint controller (selected).** Native solver runs end at
   complete write boundaries. Python assesses them and either continues the
   same cloud, accepts the case, or prepares a clean corrected attempt. Reuses
   the existing diagnostics, plans and serial execution; no solver rebuild.
2. **Solver-side checking (not selected).** A new native model or function object
   would perform checks during solver execution. It reduces relaunches but adds
   OpenFOAM build and maintenance work and makes reuse harder.

## Scope, exclusions and existing code

This is a controller addition involving configuration, initial meshing, staged
native execution, incremental assessment, cleanup, and restart bookkeeping.
Reuse plumetools.refinement for physical cell checks and executable correction
plans; keep Cai/radial geometry in plumetools.radialstudy. Do not edit the
hyStrath port or add AMR, MPI, parallel queues, a different physical expansion
model, or an automatic full scientific run.

Today initial_mesh() enforces a predicted MFP/geometric bound and can shrink the
angle automatically. run_case() runs to a fixed final time. The assessment
requires all six configured sampling blocks before calling them complete.
run_study() skips cases only when their native end time has been reached.
Those assumptions need explicit controller-aware paths; do not weaken the
existing ordinary-run or immutable-input checks to make the new path work.

## Geometry and initial mesh

Use the existing square expanding streamtube and specular side boundaries.
The 2.5-degree half-angle applies to each coordinate-plane side of that tube;
its diagonal corner has a larger angle. This change does not make a circular
axisymmetric wedge. Retain the existing meaning of outer_over_D: the outlet
plane's axial/radial model coordinate is 10 D_e. Its flat-plane corners are
slightly farther from the origin; document that distinction.

Proposed default: 100 uniformly spaced radial layers and one cell across each
transverse direction, giving exactly 100 cells. Uniform radial spacing keeps
predicted molecules per cell reasonably balanced in this expansion. Add an
explicit initial-cell-count mesh option rather than misusing max_cells as a
coarsening instruction. Keep the existing predicted-target initial mesh mode
available for old configurations and for users who want it.

The new default initial mesh is an exploratory mesh. It is allowed to violate
the predicted and measured edge targets; it must not be labeled adequate merely
because it has 100 cells. Its input record reports predicted violations and
cost. Refinements retain the inlet, angle, outlet, exit state, and gas.
The initial-cell count is not a final-cell cap. Check the proposed refined mesh
against resource budgets before building or running it.

A read-only estimate using the current model, 100 uniform layers, 2.5 degrees
and 10 D_e predicts about 2,041 parcels for each initial case at a target of 20.
Dense cases can require much greater transverse subdivision than that initial
mesh. Preserve honest budget-limited outcomes instead of promising that all
seven final meshes can be small.

## Periodic execution and evidence

Keep the current physical write cadence and startup definition initially.
Proposed assessment interval: one native write window. Before startup ends,
continue without judging the still-filling domain. Post-startup, accumulate
disjoint reset-at-write windows on the unchanged mesh.

Proposed lightweight trust gate: at least four complete post-startup windows,
grouped into two blocks, with valid local density, temperature and mean parcel
counts in every cell and stable density/velocity/temperature between the two
block means. Use a configurable 5% relative change default, with explicit
absolute floors for quantities near zero; report the actual tolerance. This
is a sizing check, not a statistical-independence or scientific-convergence
certificate. Do not relabel the legacy six-block comparison analysis as proven
or complete when the controller stops early.

Unstable or insufficient evidence means continue the same unchanged attempt
to another write boundary. Empty or missing measurements never pass and never
produce guessed MFP-driven refinement. If adequate measurements cannot be
obtained within the configured maximum sampling/step budget, stop as unresolved.

Once the trust gate is satisfied, assess all cells with the existing MFP and
occupancy definitions plus the existing timestep checks:

- All requirements pass: record acceptance, retain final outputs and move to
  the next Kn case. The six-block maximum is no longer a mandatory run length
  for this sizing workflow.
- Trusted failure: prepare a measured correction, record its reason, then
  perform the clean restart described below.
- Missing coverage or unstable evidence: continue sampling without cleaning or
  changing the current cloud. Record unresolved if the maximum is reached.

Continue using the current separation of occupancy, mesh and timestep
corrections. In particular, increase parcels first when population is inadequate
for trusting a measured MFP. Reuse existing plans and radial refinement rather
than inventing a second mesh engine. Correct only demonstrated failures.
Proposed small retry margin: when increasing parcel population, aim for 25
mean parcels per cell while retaining 20 as the acceptance requirement. This
avoids repeated retries caused by landing exactly on a noisy threshold.

Native segments must end normally at complete write boundaries; do not SIGKILL
the solver for routine assessments. The installed solver has no arbitrary
alternate-controlDict command-line option. The builder may use a controlled
temporary endTime override with an exact recorded canonical/active file hash
and crash recovery. Validate immutable inputs before staging, record the
specific override before launching, and restore the canonical dictionary before
postprocessing. An unrelated edit must still refuse restart. A segment reaching
its temporary endTime is not proof that the case has been accepted.

## Attempts, cleanup and persistent summary

Retain the existing revision-directory organization. Each corrected attempt
has its own immutable settings/identity and fresh initialization. A correction
prepares the next revision and verifies its plan while the source evidence
still exists. Before any cleanup, durably record the assessment, source/plan
hashes, proposed changes, and the next attempt's identity.

One append-only human-readable study log lives above all Kn/attempt folders:

    Cases/adaptive-study.log

Append short timestamped events for start, checkpoint assessment, continuation,
retry, cleanup, acceptance, budget limit, and solver failure. Every assessment
includes Kn, attempt/revision, simulated time, elapsed wall time, steps, cell count,
particle weight, timestep, minimum mean parcels, maximum edge/MFP ratio,
mesh/occupancy/timestep failed and unresolved counts, and the action/reason.
Retry entries include old and new settings. Print concise progress to stderr.
Never truncate this log when restarting the shell command or cleaning an attempt.

Each Kn root keeps a small atomically updated controller-state.json with physical
identity, current/next attempt, active segment hashes, outcome and transition
status. This file supports recovery; it is not a substitute for the readable log.
Keep compact per-attempt settings, assessment and correction-plan records.
Retain the original dictionaries and radial mesh metadata. Clean only outputs
owned by an automatically superseded attempt: generated polyMesh/C/V, numeric
time outputs, initialized cloud, sampling inventories and known native runtime
output directories. Record any known native log's useful failure context before
discarding large logs. Preserve the dedicated summary and all compact records.

Never clean a passed attempt, a solver-crashed attempt, an unrelated file, an
unrecognized case tree, or the user's existing first-run cases automatically.
Cleanup is confined to resolved, controller-owned paths below the current Kn
root. The next revision must be fully prepared and the transition committed
before cleaning; Allrestart can finish an interrupted transition without
creating another revision or losing the summary.

## Serial scripts, recovery and budgets

Keep Allrun/Allrestart as the simple study entry points. For the new study
settings they dispatch the controller, one Kn case at a time, using the newest
active attempt rather than queuing cleaned superseded attempts for rerun.
Allrun starts untouched cases; Allrestart skips controller-accepted cases,
recovers controller transitions, resumes the same unchanged attempt where
possible, and starts untouched cases. Legacy configurations keep their current
fixed-duration behavior. Never identify an early-stopped accepted case merely
from the latest numeric time: verify its recorded acceptance and input/mesh/
measurement identity. Existing complete-checkpoint safeguards remain in force.

Proposed controller defaults: enabled for the new study settings, a maximum of
five automatic retries per Kn, and checks every write window. Keep existing
cell, predicted-parcel and timestep budgets; also account for solver steps
cumulatively across attempts so retries cannot silently multiply the step
budget. Check each proposed segment and correction before executing it.
When a retry, sampling, cell, parcel or cumulative-step limit is reached, log a
non-passing outcome and stop the serial study with a clear reason. A solver
failure likewise stops the queue and preserves its outputs for diagnosis.
Do not silently apply the old budget-exceedance flag to the retry limit.

Both scripts retain read-only --dry-run support. It reports geometry, initial
cells, current attempt/outcome, next action, retry count and resource estimates.
It never launches a solver, creates a log, cleans files, or changes state.

Changed geometry and sampling policy require a new study tree. Generate the
new defaults under a fresh output tree and provide migration instructions; do
not overwrite or automatically clean the owner's existing study run. Physical
source mapping and breakdown search are reused and evaluated on the new domain;
a crossing outside 10 D_e is reported as not reached in this domain.

## Validation and decisions left to the builder

Use focused tests for geometry/counts, read-only dry runs, accepted/failed/
unresolved evidence, selective corrections, retry limits, persistent logging,
cleanup containment, crash recovery, identity checks and serial failure stops.
Use a tiny native WSL test to verify staged stops, restart without reinitializing
an unchanged mesh, a clean corrected retry, and no mutation of a finished case.
Run the existing relevant suites and report pre-existing failures separately.
Do not run the full scientific matrix as part of implementation validation.

The named defaults for initial radial distribution, trust threshold, retry
margin and retry limit are proposals by the implementer, not prior owner
decisions. The builder chooses module layout, strict config key names, exact
native-staging mechanism, absolute stability floors, and log formatting while
preserving the behavior and invariants above. Any need to change the geometry,
physical inlet or acceptance requirements returns to the owner.

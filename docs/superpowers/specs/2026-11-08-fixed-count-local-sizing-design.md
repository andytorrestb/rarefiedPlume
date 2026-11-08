# Fixed-count local DSMC sizing study

Status: proposed for approval. Approach1 selected by the owner on2026-11-08.

## Purpose and agreed scope

Obtain useful cell-size and parcel guidance quickly at the seven existing exit
Kn values. Keep approximately100 cells by shrinking the simulated patch instead
of increasing its cell count. The owner permits changing both radial length and
half-angle, while keeping the upstream inlet location and exit-based Kn fixed.
Do not stop the current run or modify its case inputs, outputs or controller state.

This is a focused extension of the existing sizing controller. Reuse the source
model, measured diagnostics, native checkpoint execution, clean-retry mechanism,
parallel case runner, persistent summary and input/checkpoint safeguards.

## Approaches

1. **Individually sized local patches (selected).** Predict a suitable patch for
   each Kn from the existing expansion model, then shrink it if measured cell
   requirements fail. This retains useful Kn-dependent length-scale guidance.
2. **One common tiny patch.** Size all cases for the densest condition. Geometry
   is directly comparable but unnecessarily small for the rarefied cases.

## Geometry and corrections

- Exactly100 uniform radial cells, one cell in each transverse direction, on
  every attempt. Keep the existing square streamtube and boundary models.
- Preserve gas, physical reference, exit Kn and prescribed upstream inlet state
  and location. Adjust outlet position and half-angle for each patch.
- Use predicted local MFP to choose the initial radial span and angle. Proposed
  maximum radial span:1 D_e; proposed maximum half-angle:2.5 degrees. These are
  ceilings, not required sizes. Dense cases can be much shorter and narrower.
- Select both dimensions using actual slanted/transverse edge lengths and the
  existing0.8-lambda target. Check the existing model over the proposed patch,
  rather than assuming every cell has the inlet MFP.
- Once measurements are usable, a spatial failure contracts the patch with
  conservative bounds from those measurements. Rebuild100 cells, recalculate
  particle weight and timestep, and initialize a fresh cloud. Never respond by
  adding radial layers or transverse subdivisions in this mode.
- Occupancy or timestep failures retain their separate corrections. Missing or
  unstable measurements cannot cause a pass or an invented MFP-based correction.
- Domain changes are explicitly local-study attempts, not fixed-domain mesh
  convergence tests. Preserve the existing fixed-domain identity checks in the
  original controller; do not globally relax them.

## Runtime and proposed defaults

Recompute transit time, startup, write cadence and maximum sampling duration for
each new patch. Do not inherit the10D_e case's physical run duration. Retain the
current startup and sampling-transit factors and four-window5% stability gate.
Acceptance remains edge<=0.8 measured MFP, mean parcels>=20, and the existing
collision/transport timestep limits.

Proposed cost controls, chosen by the implementer for owner review:

- Aim for25 mean parcels/cell from the first attempt to reduce threshold retries.
- Maximum5,000 predicted parcels per case; check measured population when
  available. Population is budgeted separately from the100-cell count.
- Maximum50,000 cumulative solver steps and five retries per Kn.
- A five-minute wall-time budget checked at safe write boundaries. This is a
  checkpoint budget, not a hard real-time termination guarantee.
- If a budget prevents sufficient evidence or a passing result, stop and log
  the unmet requirement. Do not relax the numerical criteria to force success.

Record solver time separately from assessment/file-I/O time, so remaining
overhead can be identified rather than attributed entirely to cell count.

## Interpretation, outputs and coexistence

A passing result describes numerical adequacy within the small sampled patch.
It does not validate the full expansion or a downstream breakdown surface.
Retain the existing analytical breakdown calculation as a separately labeled
prediction over a reference interval extending to25D_e; explicitly report when
its crossing lies outside the simulated patch. Do not claim a measured crossing
where the patch provides no such evidence.

Use separate settings and a fresh `Cases-local` output tree. The current
`Cases-adaptive` tree and active processes remain unchanged. Use the existing
Allrun/Allrestart with an explicit --cases selection for the new tree.

Keep an append-only `Cases-local/adaptive-study.log` across retries, including
Kn, attempt, cell count, inlet/outlet, half-angle, MFP, edge/MFP ratio, occupancy,
particle weight, timestep, simulated duration, solver/assessment/I-O wall times,
correction reason and final outcome. Retain compact per-attempt records and
plans; preserve unknown files and accepted results during owned cleanup.

## Validation and builder choices

Verify constant100-cell counts, dimension contraction, preserved inlet/exit
reference, recomputed durations, population/step/time limits, persistent logs,
legacy behavior, and interrupted-retry recovery. Run a bounded native local
pilot in temporary storage; do not launch the full matrix or interfere with
the existing run. Report observed performance without promising a speedup that
has not been measured.

The builder chooses configuration/API names and the bounded contraction
algorithm, reusing the existing geometry/MFP calculations. Minimum representable
patch dimensions must fail explicitly rather than produce invalid meshes.

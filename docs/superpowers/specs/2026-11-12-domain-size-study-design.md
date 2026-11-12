# Domain-size study for the Kn_e 1e-3 breakdown surface — design

Date: 2026-11-12
Status: approved in conversation; awaiting written-spec review

## Purpose

The 3D fitted-ellipsoid run (`hpc-kn1em3-ellipsoid`, job 79594) stopped in
`prepare-dsmc` with `return_mass_fraction = 0.9873 (limit 0.05)`. Two causes
are mixed together:

1. **Topology (not fixed by domain size).** The ellipsoid's vertex touches the
   source plane at one point. 99.87 % of the source area is injected directly
   by DSMC, and that gas flows into the CFD ellipsoid upstream and out again
   downstream, so reverse ≈ forward Maxwellian mass. The `fitted-capsule`
   mode, whose cylinder covers the source and which screens candidates on
   return and subsonic fractions, addresses this.
2. **Domain size and resolution (this study).** The source-connected
   continuum region reached x = 6.3 D at the 1× threshold and 7.6 D at 2×, in
   an 8 D domain. The retained volume grew 5× between those thresholds. The
   raw boundary lay in graded outer cells beyond the 2 D × ±1 D uniform core
   and was ragged (8,827 disconnected pockets). Both the capsule's feasible
   set and its screen depend on this surface.

**Goal:** find the smallest computational domain for which (a) the breakdown
surface is domain-independent and fully resolved in uniform cells, and (b) a
`fitted-capsule` interface passes the return and subsonic screen with a shape
that is stable across domain sizes. The chosen rung then proceeds to
`prepare-dsmc` for the check on the actual mesh.

Out of scope: DSMC runs, resolution convergence (D/18 is accepted for
sizing), symmetry/quarter geometry, changes to breakdown thresholds or
suitability criteria.

## Fixed choices (from the brainstorming session)

| Decision | Choice |
|---|---|
| Mesh scaling | Fixed cell size; the breakdown surface must lie in uniform cells |
| Cell budget | ≤ 3 M cells per case, full 3D (no symmetry planes) |
| Core | D/18 (0.0556 D), uniform core x ∈ [0, 8.5 D], \|y\|,\|z\| ≤ 2.75 D (snapped to 2.78 D), identical in every rung. D/16 was rejected at pre-flight: its staircase source disk is +3.45 % in area (D/18: +0.60 %), and the CFD source is not rescaled. |
| Sweep | One proportional ladder, half-width = x_max / 2 |
| Convergence | Successive-rung change < 5 % of every 1×-threshold metric |
| Interface | `fitted-capsule` in every rung, CFD-only through `post-cfd` |

## Case matrix

Common to every rung: `source.knudsen: 1.0e-3`, `mesh.geometry: cartesian`,
`core_cell_over_D: 0.0556`, `core_x_over_D: 8.5`, `core_half_over_D: 2.75`,
`outer_cells: 16`, `max_outer_expansion: 250.0` (the x48 axial segment needs
214.6; with the default 8 the planner puts an adjacent-cell jump of up to 10×
at the core edge). This gives divisions 169 × 132 × 132 = **2,944,656 cells**
in every rung. The numerics copy `hpc-kn1em3-capsule` (Minmod,
`viscosity_cap {enabled: true, factor: 1.0}`,
`breakdown.retained_stability_limit: 0.02`, capsule fit
`clearance_cells: 0.25`), and `cfd.end_time_transits: 1.0`
(transit = x_max / U_e).

| Profile | x_max / D | half-width / D | Max adjacent growth (planner) | End time | Est. steps (Δt ≈ 6.0e-7 s) |
|---|---|---|---|---|---|
| `domstudy-kn1em3-x10` | 10 | 5 | 1.112 | ~2.8 ms | ~4.7 k |
| `domstudy-kn1em3-x16` | 16 | 8 | 1.248 | ~4.5 ms | ~7.6 k |
| `domstudy-kn1em3-x24` | 24 | 12 | 1.327 | ~6.8 ms | ~11 k |
| `domstudy-kn1em3-x32` | 32 | 16 | 1.373 | ~9.0 ms | ~15 k |
| `domstudy-kn1em3-x48` | 48 | 24 | 1.430 | ~13.6 ms | ~23 k |

Wall-time estimates span 0.5–2.3 h (ellipsoid-run rate, 1.25e-7 s per cell
per step) to 8–37 h (capsule-run rate, 2e-6 s per cell per step) on 192
ranks, within the 120 h `longq` walltime. The 16× gap between the two
measured rates is noted, not investigated, in this study.

**Known confounder:** outer grading coarsens as the domain grows. Any change
made to keep the largest rungs stable (for example the viscosity cap factor)
is recorded in the profile comment and the study report.

The running `hpc-kn1em3-capsule` case (30 × ±15 D, D/20 core of 3 × ±1.5 D)
is a side reference only. Its surface lies in graded cells, so it is not a
rung.

## Execution

1. **Profiles:** five entries in `cases/cai2012-hybrid/hybrid.yaml`, each
   written out in full (no YAML merge keys), labelled
   `"DOMAIN STUDY: Kn_e 0.001, x_max N D (D/18 core 8.5 x ±2.75 D, sizing only)"`.
2. **Wrapper:** `cases/cai2012-hybrid/run3DCase.sh` gains an optional
   `STAGES` variable: a list (space- or colon-separated) run in order, stopping at the
   first failure. When `STAGES` is unset, the existing `STAGE` behaviour is
   unchanged. The study uses `STAGES=setup:cfd:post-cfd`.
3. **Pre-flight (login node, Python only):** `Allrun <case> setup --profile
   <p>` for all five rungs. Each must report n_cells ≤ 3,000,000, an
   expansion within its `max_outer_expansion`, and a memory estimate under
   the guard. The real planner values replace the estimates in this spec's
   table and in the profile comments.
4. **Submission:** all five at once, with `msub`, never `qsub`:
   `msub -v PROFILE=domstudy-kn1em3-xNN,STAGES=setup:cfd:post-cfd run3DCase.sh`.
   Case directories default to `runs/<profile>`.

## Analysis: `plumetools/hybridplume/domainstudy.py`

A new module with a CLI subcommand,
`python -m plumetools.hybridplume domain-study <case> [<case> ...] --out <dir>`.
It reads only the files each rung's `post-cfd` already writes (`breakdown.npz`,
`breakdown.json`, the capsule fit and screen output, `setup.json`). The first
plan task confirms those files hold cell centres and bounds, per-factor
retained masks (or the inputs to recompute them) and the capsule screen rows.
If they do not, the module recomputes the masks from the stored indicators
using `breakdown`'s own functions rather than duplicating them.

Per rung and per threshold factor (0.5, 1, 2):

- **Surface metrics** on the source-connected retained region: axial extent on
  the centreline; r_b(x), the maximum retained radius in slabs centred at
  x = 1, 2, 4, 6, 8 D (slab width one core cell); maximum radius; enclosed
  volume.
- **Validity:** minimum distance from the 2× retained region to the graded
  zone (the core box boundary) and to the outer domain boundary. Pass:
  ≥ 1 D to the graded zone.
- **Capsule:** R, L, a, volume, screen status, return_mass_fraction,
  subsonic_area_fraction of the selected candidate.

Across rungs, ordered by x_max:

- Relative change of each 1×-threshold surface metric and of capsule volume
  versus the next-larger rung: |m_k − m_{k+1}| / |m_{k+1}|. For r_b(x) this
  is evaluated per station and the maximum is reported. A station where either
  rung has no retained cells counts as a failure, not as zero.
- **Verdict per rung:** `valid` (validity check), `converged` (all changes
  < 5 %; undefined for the largest rung), `capsule_passes` (screen status
  `passed`).

Outputs in `--out`: `domain_study.json` (all numbers and verdicts),
`domain_study.md` (tables) and `domain_study.png` (z = 0 overlay of each
rung's 1× raw boundary and selected capsule section, core box drawn).

Tests (pytest, TDD) use small synthetic voxel fields with known answers:
a retained cone whose extent and radii are analytic; a truncation case where
the region touches the graded zone (validity fails); two rungs differing by a
known percentage (convergence boundary at 5 %); and a missing station
(counts as failure). A config test asserts that every `domstudy-*` profile
plans to ≤ 3 M cells with an identical core.

## Decision rule

- **Chosen domain** = the smallest rung that is `valid`, `converged` against
  the next-larger rung, and `capsule_passes`. It proceeds to `prepare-dsmc`
  (separate follow-up, not in this plan).
- **Core containment fails** in the larger rungs → the core is too small at
  3 M full-3D cells; revisit with a larger core budget or quarter geometry.
- **Surface converges but no capsule passes** → domain size was not the
  limiting factor. The report states which criterion fails (return or
  subsonic) and by how much; that becomes a separate investigation.
- **No convergence by x48** → report the trend; the study does not extend
  automatically.

## Risks

- Far-field emptying and crashes in the largest rungs (seen in job 79485).
  Mitigation: viscosity cap from the capsule profile; resumable stages; at
  1 transit the run ends earlier than the 2.75-transit crash.
- Aggressive outer grading (growth ~1.35) may degrade far-field stability or
  reflect waves inward. Visible as a rung-to-rung trend in the surface
  metrics; recorded, not corrected.
- One transit may not settle the near field in every rung; the existing
  `retained_stability_limit` (2 % between the last two writes) flags this in
  `post-cfd`.

# Fixed-count Local Sizing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Keep100 cells while choosing/shrinking a local DSMC patch, reducing populations and simulated durations without touching the active study.

**Architecture:** Add an opt-in local policy and focused geometry/revision module. Reuse the existing generation/controller/parallel runner with a common local-study identity, while keeping each attempt's actual geometry and sampling identity distinct. Keep existing fixed-domain validation unchanged.

**Tech Stack:** Python, NumPy/SciPy, YAML, native OpenFOAM v2512 hyStrath, bash/WSL, pytest.

**Spec:** `docs/superpowers/specs/2026-11-08-fixed-count-local-sizing-design.md`.

## Global Constraints

- Owner approved separately sized patches per Kn; no wall-clock cutoff.
- Exactly100 radial cells,1x1 transverse, fixed physical inlet and exit reference.
- Cap initial radial span at1D_e and half-angle at2.5degrees; actual edges <=0.8 local MFP.
- Start near25parcels/cell, accept>=20;5000 predicted parcels,50000 cumulative steps,5 retries.
- Recompute durations when geometry changes; retain frozen duration on unchanged geometry.
- Separate Cases-local output/settings and persistent log; never touch current Cases-adaptive files/processes.
- Report analytical breakdown through25D_e separately; no full-plume convergence claim.

## Review Focus

1. Loading a resized child must retain its recorded source, geometry and sampling without resetting its dimensions.
2. Recovery and queue grouping must accept local-domain changes only for identical fixed inlet/gas/Kn/policy.
3. Retry budgets must use the child's shortened cadence, not the parent's old duration.
4. Domain contraction must never add cells or silently accept invalid/sparse measurements.
5. Old fixed-domain cases and active processes must retain their existing behavior.

### Task1: Local geometry and initial generation

**Files:** config.py; new local.py; generation.py; study-local.yaml; generate_local_cases.py; local unit tests.

**Interfaces:** `Local` settings enabled=False,max_span_over_D=1,initial_parcels=25,reference_outer_over_D=25. `initial_local_model(config,kn)` returns resolved RadialModel; `contract_local_model(model,mesh,target)` returns a smaller model and100-cell mesh. `study_identity(record)` in case.py returns a common local identity or the legacy physical identity.

- [ ] Add failing tests for seven100-cell predicted-target patches, bounded dimensions/populations/steps, no writes on dry-run, explicit invalid settings and legacy compatibility.
- [ ] Resolve initial dimensions from the existing model over the proposed patch with a small safety margin; check actual geometry; serialize resolved settings.
- [ ] Generate initial weight for25 parcels while keeping acceptance20. Write fresh Cases-local via separate settings/wrapper; do not change old defaults.
- [ ] Run tests and commit.

### Task2: Measured contractions and recoverable controller integration

**Files:** local.py, case.py, generation.py, adaptive.py, run.py, segments.py; unit and bounded native tests.

**Interfaces:** `revise_local_case(parent,plan,correction,delta_t_s,dry_run)` follows existing revision return shape, including child numerics/sampling; atomic publication includes audit records. Local study identity binds inlet state/location, fixed physical inputs/Kn, target and dimensionless sampling policy, not changing outlet/angle.

- [ ] Add failing tests for100-cell shrink retries, fresh sampling on geometry change, preserved inlet, population-only unchanged duration, different-Kn/inlet/policy refusal, grouping/recovery and unchanged active-tree hashes.
- [ ] Dispatch only enabled local cases; keep the legacy checks. Interpolate parent measured density to estimate child populations, recompute new-domain sampling, and use measured timestep guidance conservatively.
- [ ] Adapt root identity checks and child-budget cadence. Preserve publication-before-clean and exact source/ownership checks.
- [ ] Add dimension, MFP/population/duration and phase timing summaries. Separate native solver wall time from controller/assessment/I-O overhead; no time termination logic.
- [ ] Run focused tests and a bounded native local retry, record real costs, commit.

### Task3: Deliver a runnable separate study

**Files:** README/tutorial/CLI/API, mkdocs navigation, dedicated bounded verification script/report.

- [ ] Verify actual generation and parallel-wrapper dry runs in Cases-local; build all seven100-cell meshes without launching the matrix.
- [ ] Run relevant legacy/local suites and docs checks; record known baseline limitations rather than expanding scope.
- [ ] Obtain one independent focused review, address material findings, and commit changes atomically.
- [ ] Integrate compatible additions into the normal checkout, preserving unrelated edits. Generate/build Cases-local there; never stop current processes or change their Cases-adaptive inputs/state.

## Execution record

Implement inline in the existing isolated worktree. The source code for the
running study is in the normal checkout; do not modify it until compatibility
checks pass. No new full scientific matrix, wall-time cap, or extra worktree.

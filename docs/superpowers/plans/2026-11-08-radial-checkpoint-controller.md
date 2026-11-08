# Radial Checkpoint Controller Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the serial radial sizing study from 100-cell initial meshes with periodic measured corrections, clean retries, and a persistent study log.

**Architecture:** Keep the native solver and shared refinement engine. Add incremental sizing assessment, a recoverable native segment runner, and a serial controller with per-Kn state and one append-only study log. Existing configurations default to the legacy fixed-duration path.

**Tech Stack:** Python, NumPy/SciPy/PyYAML, OpenFOAM v2512 hyStrath dsmcFoam+, bash/WSL, pytest.

**Spec:** `docs/superpowers/specs/2026-11-08-radial-checkpoint-controller-design.md` (owner approved).

## Global Constraints

- 2.5-degree half-angle, outlet plane at 10 D_e, fixed existing upstream inlet.
- 100 uniform radial layers, one transverse cell in each direction initially.
- Maximum edge <= 0.8 measured local MFP; mean parcels >= 20; unchanged timestep limits.
- Four post-startup windows, two blocks, 5% per-cell change gate; no statistical-independence claim.
- Aim for 25 parcels when increasing occupancy; at most five automatic retries.
- One serial native solver at a time. No solver source edits or full scientific matrix.
- Preserve legacy runs, immutable inputs, complete-checkpoint checks, and portable refinement API.
- Never clean before the corrected child and durable retry record exist; never clean accepted/crashed/unowned runs.
- `Cases-adaptive/adaptive-study.log` is append-only and survives cleanup and command restarts.

## Review Focus

1. Interrupted temporary controlDict overrides must restore only an exact known file; unrelated edits refuse recovery (Task 3).
2. A segment reaching its temporary endTime must not be skipped as an accepted case (Tasks 3/4).
3. Cleanup encountering symlinks or unrelated files must not escape the owned attempt or delete user data (Task 4).
4. Restart during child-publication/cleanup must not create duplicate retries or truncate the log (Task 4).
5. Sparse or unstable measurements and exhausted budgets must never produce a numerical pass (Tasks 2/4).

## File responsibilities

- `config.py`, `mesh.py`, family `study.yaml`: optional initial cell count and controller settings; new defaults.
- `post.py`: share measured calculations; separate partial-run sizing trust from legacy six-block analysis.
- New `segments.py`: audited temporary controlDict, normal native stops, exact checkpoint validation, recovery.
- New `adaptive.py`: controller state, summary log, owned cleanup, retries and early acceptance.
- `run.py`, `__main__.py`, family wrappers: dispatch new/legacy workflows without weakening either.
- `generation.py`: expose numerical costs on dry revision plans; preserve immutable revision creation.
- Unit tests grouped by settings/assessment/segments/controller; native test exercises the actual scripts.
- Family README, tutorial, CLI/API/method documentation: new defaults, output tree, log and outcome semantics.

### Task 1: Cheap initial mesh and strict controller settings

**Files:** `plumetools/radialstudy/config.py`, `mesh.py`, family `study.yaml`, `generate_cases.py`, unit configuration/mesh tests.

**Interfaces:** Add `Mesh.initial_cells: int | None = None` and `Controller` on `RadialStudyConfig`, serialized by settings(). Controller fields: enabled=False, check_every_windows=1, min_windows=4, stability_relative=.05, density_floor=1., velocity_floor=1., temperature_floor=1., retry_parcels=25., max_retries=5. New family settings enable it and use `Cases-adaptive` by default; legacy files omit it.

- [ ] Write tests: enabled defaults generate exactly100 cells at2.5 degrees and10D_e; explicit predicted mode still enforces predicted targets; invalid counts/tolerances rejected; old case records load unchanged.
- [ ] Run those tests and observe expected failures before code changes.
- [ ] Add strict parsing and optional coarse-mesh branch to initial_mesh(model). Report exploratory predicted edge failures in case sizing. Pin existing test model helpers to the old defaults where they characterize legacy behavior.
- [ ] Run configuration/physics/mesh/case/generation tests; commit this independently useful mesh/config change.

### Task 2: Incremental numerical sizing assessment

**Files:** `plumetools/radialstudy/post.py`, new `tests/unit/test_radialstudy_incremental.py`.

**Interfaces:** `checkpoint_assessment(case) -> tuple[dict, RefinementPlan | None]`; internal `_compute(case, *, incremental=False)` shares moments and physical checks. Return recommended_delta_t_s and the exact correction plan when usable; missing evidence returns unresolved and no executable guess.

- [ ] Write tests with actual generated fields: fewer than4 windows unresolved; unstable per-cell n/U/T unresolved;4 stable windows usable; low occupancy requests parcels; edge/timestep failures selected; no new independence/convergence certificate.
- [ ] Run tests RED.
- [ ] Use the latest min_windows windows for the two-block stability gate, retain disjoint write provenance, and use measured moments on the current attempt for diagnostics. Floors are1 m^-3,1 m/s,1 K; compare per-cell n,T and velocity vector norms. Preserve legacy _statistics() and full-run comparisons.
- [ ] Run incremental, post, report and sampling tests GREEN; commit.

### Task 3: Recoverable native checkpoint segments

**Files:** new `plumetools/radialstudy/segments.py`, new segment unit tests, native mechanics tests.

**Interfaces:** `recover_control(case) -> None`, `run_segment(case, end_s, *, runner=None) -> dict`. Use existing native run/environment/mesh/checkpoint helpers. Segment results include start_s,end_s,steps,status.

- [ ] Write tests: fresh segment initializes once; continuation preserves the cloud; normal stop validates exact complete checkpoint; native failure restores canonical inputs; interrupted override recovers exact known bytes; unrelated edits refuse; temporary completion is not full-run acceptance.
- [ ] Observe RED.
- [ ] Store canonical control bytes plus active hash in an atomically published recovery manifest before replacing controlDict. Only endTime changes. Restore in finally and on next invocation after process loss. Verify all other immutable inputs. Record native segments and unique logs in existing run-state.json; check cadence alignment and actual checkpoint after success.
- [ ] Run segment/checkpoint/legacy-run tests and one bounded native staged stop; commit.

### Task 4: Automatic corrections, persistent log and bounded cleanup

**Files:** new `plumetools/radialstudy/adaptive.py`, `generation.py`, new controller unit tests.

**Interfaces:** `adaptive_preview(case, *, resume=False) -> dict`, `run_adaptive(case, *, resume=False, allow_budget_exceedance=False, progress=None) -> dict`. Per-Kn `controller-state.json`; global `Cases-adaptive/adaptive-study.log`. Preview is read-only.

- [ ] Write tests for accepted early stop, unresolved continuation, occupancy->mesh/timestep corrections, retry cap, cumulative charged step budget, dry-run immutability and append-only log contents.
- [ ] Add cleanup tests for unknown files/symlinks/outside paths, accepted/crashed/unowned cases, prepared-child-before-clean ordering, and restart between transition phases without duplicate child creation.
- [ ] Observe RED, then implement the state machine. State contains owner id, physical identity, attempt records, charged steps, pending segment/transition and final outcome. Charge requested native segment steps before launch conservatively; do not lose charges on crashes.
- [ ] At trusted failure use the existing measured plan/revise_case. Override only the needed occupancy recommendation to target25; keep acceptance20. Dry revision plans expose numerical estimates so budgets are checked before publication/build/run. Prepare child, persist summary/plan/transition, clean only allowlisted owned output paths, switch active revision and initialize fresh.
- [ ] Persist per-attempt summary and hashes. Append UTC timestamped concise events (Kn,attempt,time,wall time,steps,cells,weight,dt,minimum occupancy,maxedge/MFP,failed/unresolved counts,action/reason and settings changes). Flush/fsync the log before cleaning; atomically save state.
- [ ] Verify accepted measurement/input/mesh identity when resuming; retain accepted outputs. Budget/unresolved/solver failures stop the study with actionable reasons and retained evidence. Run tests GREEN; commit.

### Task 5: Serial script integration and migration

**Files:** `run.py`, `__main__.py`, family `Allrun`, `Allrestart`, `generate_cases.py`, suite tests and `tests/openfoam/test_radialstudy_openfoam.py`, relevant documentation.

**Interfaces:** run_study keeps its existing signature. Group controller-enabled revisions by Kn root and dispatch once to their active attempt. Legacy cases keep their existing per-revision queue. Reject mixed incompatible identities. Family wrappers default to fresh `Cases-adaptive`; explicit --cases selects other trees.

- [ ] Write tests for a mixed serial study's preflight, cleaned attempts not requeued, accepted early-stopped cases skipped, Allrestart recovery, failure stops, default output tree, and script dry runs from a different directory.
- [ ] Observe RED; add dispatch and fresh-tree defaults. Existing Cases remains untouched. New defaults' predicted breakdown is evaluated only within10D_e.
- [ ] Add native bounded integration covering multiple actual segments, unchanged-cloud resume, one clean corrected retry, persistent log, and accepted skip without mutation. Use small configuration overrides rather than full scientific cases.
- [ ] Update README/tutorial/CLI/API/method documentation, including log path, stop outcomes, stability-vs-convergence distinction, resource limitations and how to run existing legacy Cases explicitly.
- [ ] Run focused suites, bounded native tests, full non-OpenFOAM suite and docs checks. Record known baseline failures by name; do not expand scope to fix unrelated HPC/nozzle fixtures.
- [ ] Generate the seven new dictionary-only defaults into the fresh Cases-adaptive tree, inspect a dry-run, and commit integration/docs. Obtain one independent whole-change review, resolve material findings, and integrate locally into main without adding intermediate branches.

## Execution

Implement directly in the existing isolated worktree, keeping the owner's prior
preference for atomic commits and a clean branch list. The design is approved;
the implementation follows these steps without further routine permission gates.

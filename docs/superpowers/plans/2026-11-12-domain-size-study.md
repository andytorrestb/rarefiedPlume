# Domain-Size Study Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Five CFD-only rungs of increasing domain size (x_max = 10, 16, 24, 32, 48 D) at Kn_e 1e-3 with an identical D/18 core, plus a `domain-study` analysis that decides the smallest domain whose breakdown surface is domain-independent and whose capsule interface passes the screen.

**Architecture:** Profiles in `hybrid.yaml` reuse the existing `setup → cfd → post-cfd` stages unchanged; `run3DCase.sh` gains a `STAGES` list so jobs stop after `post-cfd`. A new module `plumetools/hybridplume/domainstudy.py` has a pure array layer (metrics, comparison, verdicts) and an IO layer that reads each rung's `results/` and CFD mesh, recomputing the 0.5× and 2× retained masks with `breakdown`'s own functions.

**Tech Stack:** Python 3.11 (repo `.venv`), numpy, scipy (via `breakdown`), matplotlib (optional, plots), pytest; bash; OpenFOAM v2512 `rhoCentralFoam`; Moab `msub`.

**Spec:** `docs/superpowers/specs/2026-11-12-domain-size-study-design.md`

## Global Constraints

- Run all Python with the repo venv: `source .venv/bin/activate` from the repository root (`/scratch/torresfa/OpenFOAM/rarefiedPlume`).
- Submit jobs with `msub`, never `qsub`.
- Every rung: `source.knudsen: 1.0e-3`; mesh `core_cell_over_D: 0.0556`, `core_x_over_D: 8.5`, `core_half_over_D: 2.75`, `outer_cells: 16`, `max_outer_expansion: 250.0`; planner result 169 × 132 × 132 = 2,944,656 cells, ≤ 3,000,000.
- Threshold factors 0.5, 1.0, 2.0; stations x = 1, 2, 4, 6, 8 D; convergence tolerance 0.05 (relative); validity requires the 2× region ≥ 1.0 D from the graded zone.
- Do not modify breakdown thresholds, `interface.criteria`, or existing profiles.
- CFD-only: no rung runs `prepare-dsmc` or later in this plan.

## Review Focus

1. **A rung whose `post-cfd` breakdown failed the stability check** (no `breakdown.npz`, only `breakdown.failed.json`). Expect status `unsettled` with the failure text, no crash, and the comparison treating it as missing. Test: Task 4 `test_load_rung_unsettled`.
2. **A rung whose CFD crashed or never ran** (no `results/breakdown*`). Expect status `missing`; the study still reports the other rungs. Test: Task 4 `test_load_rung_missing`, Task 5 `test_study_with_missing_rung`.
3. **Capsule fit not `fitted`** (`no_feasible_capsule`) or screen status `no_passing_candidate`. Expect `capsule_passes: false` with the status shown, not a KeyError. Test: Task 4 `test_capsule_summary_states`.
4. **Recomputed 1× mask ≠ saved `retained`** (config drift between breakdown and analysis). Expect a `WorkflowError`, never a silent report. Test: Task 4 `test_load_rung_rejects_mask_drift`.
5. **Empty retained region at 0.5×, or a station with no retained cells.** Expect `None` metrics and an infinite change (a failure), never a division by zero. Test: Task 3 `test_empty_region_metrics_are_none` and `test_missing_station_is_failure`.

---

### Task 1: Domain-study profiles

**Files:**
- Modify: `cases/cai2012-hybrid/hybrid.yaml` (append after `hpc-kn1em3-capsule`)
- Test: `tests/unit/test_hybridplume_domstudy_profiles.py`

**Interfaces:**
- Produces: profiles `domstudy-kn1em3-x10`, `-x16`, `-x24`, `-x32`, `-x48`.

- [ ] **Step 1: Write the failing test**

```python
NAMES = {10: 5.0, 16: 8.0, 24: 12.0, 32: 16.0, 48: 24.0}   # x_max -> half-width

@pytest.mark.parametrize("x", NAMES)
def test_domstudy_profile_plans(x):
    cfg = config.load(STUDY, f"domstudy-kn1em3-x{x}")
    src = config.build_source(cfg)
    p = mesh.plan(cfg, src)
    assert p.divisions == (169, 132, 132) and p.n_cells == 2_944_656
    assert abs(mesh.predicted_source(p, src).area_error) < 0.01
    assert p.max_expansion_ratio < 1.45
    assert cfg.get("domain.x_max_over_D") == x and cfg.get("domain.half_width_over_D") == NAMES[x]
    assert cfg.get("interface.mode") == "fitted-capsule"
    assert cfg.get("source.knudsen") == 1.0e-3
    assert cfg.get("cfd.end_time_transits") == 1.0
    assert cfg.get("breakdown.retained_stability_limit") == 0.02

def test_domstudy_cores_identical():
    cores = {(p.x_segments[0].hi, p.y_segments[1].lo, p.y_segments[1].hi)
             for p in (mesh.plan(c, config.build_source(c)) for c in
                       (config.load(STUDY, f"domstudy-kn1em3-x{x}") for x in NAMES))}
    assert len(cores) == 1
```

Use the `ROOT`/`STUDY` header from `tests/unit/test_hybridplume_capsule_profile.py`.

- [ ] **Step 2: Run, expect FAIL** (unknown profile)

Run: `pytest tests/unit/test_hybridplume_domstudy_profiles.py -v`

- [ ] **Step 3: Add the five profiles**

Each is written out in full (no merge keys), copying `hpc-kn1em3-capsule`'s `cfd`, `breakdown`, `interface` (`mode: fitted-capsule`, same `fit`), `dsmc.mesh` and `resources` blocks. Per rung, change only:
`label: "DOMAIN STUDY: Kn_e 0.001, x_max N D (D/18 core 8.5 x ±2.75 D, sizing only)"`,
`domain: {x_max_over_D: N, half_width_over_D: N/2}`,
`mesh: {geometry: cartesian, core_cell_over_D: 0.0556, core_x_over_D: 8.5, core_half_over_D: 2.75, outer_cells: 16, max_outer_expansion: 250.0}`,
`cfd: {end_time_transits: 1.0, viscosity_cap: {enabled: true, factor: 1.0}}`.
Precede the block with a comment citing the spec and noting D/16 was rejected for its +3.45 % source-area error.

- [ ] **Step 4: Run, expect PASS**; also run `pytest tests/unit/test_hybridplume_capsule_profile.py tests/unit/test_config.py -q` (no regressions).

- [ ] **Step 5: Commit** — `git add cases/cai2012-hybrid/hybrid.yaml tests/unit/test_hybridplume_domstudy_profiles.py && git commit -m "feat(cai2012-hybrid): domain-size study profiles x10..x48"`

---

### Task 2: `STAGES` list in `run3DCase.sh`

**Files:**
- Modify: `cases/cai2012-hybrid/run3DCase.sh` (the final `args=`/`Allrun` lines)
- Test: `tests/unit/test_hybridplume_run3dcase.py`

**Interfaces:**
- Produces: env var `STAGES` — stage names separated by spaces or colons (`setup:cfd:post-cfd`, which survives `msub -v` comma/quote handling), run in order with the same `args`; the first non-zero exit stops the job with that status. Unset → existing `STAGE` behaviour, byte-for-byte the same command line.

Note: the working tree already carries an uncommitted user change (`PROFILE` default → `hpc-kn1em3-capsule`). Keep it. Do not commit this file until the user confirms that change may go in with it.

- [ ] **Step 1: Write the failing test**

```python
def test_run3dcase_syntax():
    assert subprocess.run(["bash", "-n", str(SCRIPT)]).returncode == 0

def test_run3dcase_has_stage_list():
    text = SCRIPT.read_text()
    assert '${STAGES//:/ }' in text                # colon or space separated
    assert 'STAGE=${STAGE:-all}' in text          # single-stage default kept
```

- [ ] **Step 2: Run, expect FAIL** on the second test.

- [ ] **Step 3: Implement.** Build the argument list after the case directory with a placeholder stage; if `STAGES` is non-empty, loop `bash "$SCRIPT_DIR/Allrun" "$CASE_DIR" "$stage" "${rest[@]}" || exit $?`, printing `>>> stage <name>` before each; else run the existing single command. Update the `printf` banner to show `STAGES` when set.

- [ ] **Step 4: Run, expect PASS.**

- [ ] **Step 5: Commit** (after user confirmation per the note) — `git commit -m "feat(cai2012-hybrid): STAGES list in run3DCase.sh"`

---

### Task 3: Metrics, comparison and verdicts (pure arrays)

**Files:**
- Create: `plumetools/hybridplume/domainstudy.py`
- Test: `tests/unit/test_hybridplume_domainstudy.py`

**Interfaces:**
- Produces:
  - `FACTORS = (0.5, 1.0, 2.0)`, `STATIONS_OVER_D = (1.0, 2.0, 4.0, 6.0, 8.0)`, `TOLERANCE = 0.05`, `MIN_CORE_CLEARANCE_OVER_D = 1.0`
  - `@dataclass(frozen=True) class Box: lo: np.ndarray; hi: np.ndarray` (3-vectors, metres)
  - `surface_metrics(lower, upper, volumes, retained, D, stations_over_D=STATIONS_OVER_D) -> dict` with keys `axial_extent_over_D` (max `upper[:,0]` over retained cells whose y and z footprints contain 0), `r_b_over_D` (dict `str(station) -> float|None`: max cell-centre radius over retained cells whose x footprint contains the station), `max_radius_over_D`, `volume_m3`. Any metric over an empty set is `None`.
  - `clearance_over_D(lower, upper, retained, box: Box, D) -> float|None`: min over retained cells of distance from footprint to the box's x_max, y_min, y_max, z_min, z_max faces (the x_min face is the source plane and is excluded); negative if a footprint pokes outside.
  - `relative_change(small, large) -> float`: `abs(small-large)/abs(large)`; `inf` if either is `None` or `large == 0`.
  - `compare(small: dict, large: dict) -> dict`: changes for the 1× surface metrics (`axial_extent_over_D`, `max_radius_over_D`, `volume_m3`, and `r_b_max` = max over stations) and `capsule_volume_m3`; inputs are rung records as produced by Task 4 (`rung["surfaces"]["1.0"]`, `rung["capsule"]["volume_m3"]`).
  - `verdicts(rungs: list[dict]) -> list[dict]`: rungs sorted by `x_max_over_D` (duplicate x_max → `ValueError`); each output adds `valid` (`clearance_core_over_D >= 1.0`), `converged` (all `compare` changes vs next-larger `ok` rung `< TOLERANCE`; `None` for the largest), `capsule_passes`, `changes`. A rung with `status != "ok"` gets all three verdicts `False` (largest: `converged None`). Returns also `chosen`: `x_max_over_D` of the smallest rung with all three `True`, else `None` — as `{"rungs": [...], "chosen": ...}`.

- [ ] **Step 1: Write the failing tests**

Fixture: a uniform lattice, `h = 0.1`, `D = 1.0`, x ∈ [0, 2], y,z ∈ [−1, 1]; `retained = (r <= 0.5 + 0.25*x) & (x_centre <= 1.6)`.

```python
def test_cone_metrics(cone):
    m = domainstudy.surface_metrics(*cone, D=1.0, stations_over_D=(0.55, 1.25))
    assert m["axial_extent_over_D"] == pytest.approx(1.6)
    assert abs(m["r_b_over_D"]["1.25"] - 0.8125) <= 0.1*np.sqrt(2)
    assert m["volume_m3"] == pytest.approx(cone_retained.sum()*1e-3)

def test_empty_region_metrics_are_none(cone):     # Review Focus 5
    m = domainstudy.surface_metrics(lower, upper, vol, np.zeros_like(ret), D=1.0)
    assert m["axial_extent_over_D"] is None and m["volume_m3"] == 0.0
    assert all(v is None for v in m["r_b_over_D"].values())

def test_clearance_ignores_source_plane():
    box = Box(np.array([0., -1, -1]), np.array([2., 1, 1]))
    # cone touches x=0 by construction; nearest other face: y/z at r=0.9 near x=1.6
    assert domainstudy.clearance_over_D(lower, upper, ret, box, 1.0) == pytest.approx(0.1, abs=0.05)

def test_clearance_negative_when_outside():
    # a retained footprint extending 0.2 past box.hi in y → -0.2

def test_relative_change_boundary():
    assert domainstudy.relative_change(1.049, 1.0) < 0.05
    assert domainstudy.relative_change(1.051, 1.0) > 0.05

def test_missing_station_is_failure():            # Review Focus 5
    assert domainstudy.relative_change(None, 1.0) == float("inf")
    assert domainstudy.relative_change(1.0, 0.0) == float("inf")

def test_verdicts_choose_smallest_converged():
    # three ok rungs x=10,16,24; 10→16 changes 8 %, 16→24 changes 2 %; all valid, capsules pass
    out = domainstudy.verdicts([r24, r10, r16])          # unsorted input
    assert [r["x_max_over_D"] for r in out["rungs"]] == [10, 16, 24]
    assert out["chosen"] == 16
    assert out["rungs"][-1]["converged"] is None

def test_verdicts_invalid_rung_not_chosen():
    # r16 clearance_core_over_D = 0.5 → valid False → chosen None

def test_verdicts_duplicate_xmax_raises():
    with pytest.raises(ValueError): domainstudy.verdicts([r10, r10])
```

Build the small rung records with a helper `_rung(x, axial, rmax, vol, rb, cap_vol, clearance=2.0, passes=True, status="ok")` in the test file, shaped like Task 4's output.

- [ ] **Step 2: Run, expect FAIL** (module missing). `pytest tests/unit/test_hybridplume_domainstudy.py -v`
- [ ] **Step 3: Implement** the interfaces above in `domainstudy.py`. Numpy only; no IO.
- [ ] **Step 4: Run, expect PASS.**
- [ ] **Step 5: Commit** — `git commit -m "feat(hybridplume): domain-study metrics and verdicts"`

---

### Task 4: Load one rung from its case directory

**Files:**
- Modify: `plumetools/hybridplume/domainstudy.py`
- Test: `tests/unit/test_hybridplume_domainstudy_io.py`

**Interfaces:**
- Consumes: Task 3 functions; `workflow.load_case`, `workflow.results_dir`, `workflow.WorkflowError`, `config.build_source`, `meshgeom.read_mesh`, `MeshGeometry.cell_size()`, `breakdown.INDICATORS`, `breakdown.continuum_mask`, `breakdown.source_connected`, `mesh.PATCHES["source"]`, `spatial._selected`.
- Produces: `load_rung(case_dir: Path) -> dict` with keys `case`, `profile`, `x_max_over_D`, `half_width_over_D`, `status` (`"ok" | "unsettled" | "missing"`), `reason` (str, for non-ok), and for `ok`: `surfaces` (`{"0.5"|"1.0"|"2.0": surface_metrics}`), `clearance_core_over_D` (2× region vs core box), `clearance_domain_over_D` (2× region vs domain box), `retained_stability` (copied from `breakdown.json`), `capsule` (see below), `boundary_xy` (z = 0 slice: list of `[x, y]` retained-cell centres of the 1× region with |z| ≤ h/2, for plotting), `core_box`, `domain_box` (each `[[lo],[hi]]` in D).
  `capsule`: `{"fit_status", "screen_status", "R", "L", "a", "volume_m3", "return_mass_fraction", "subsonic_area_fraction", "passes": bool}`; numeric fields `None` unless `fit_status == "fitted"`; `passes` true only when `screen_status == "passed"`.
  Core box: bounding box of footprints of cells whose `cell_size()` equals `core_cell_over_D*D` on all three axes (rtol 1e-6). Domain box: bounding box of all footprints.
  `status` rules: `breakdown.npz` present → `ok`; else `breakdown.failed.json` present → `unsettled` with its `failures` joined as `reason`; else `missing`. `interface.json` absent on an `ok` rung → `capsule.fit_status = "missing"`.
  Recompute masks per factor from the saved indicators with `cfg.require("breakdown.thresholds")` and `cfg.get("breakdown.criterion", ["Kn_GLL", "B_exp"])`; if the 1× recomputed mask differs from saved `retained` in any cell, raise `WorkflowError("recomputed 1x retained mask differs from breakdown.npz; rerun post-cfd")`.

- [ ] **Step 1: Write the failing tests**

Fixture (pattern from `tests/unit/test_hybridplume_ellipsoid.py:112-120`): resolved config from `domstudy-kn1em3-x10` with an overlay making `core_cell_over_D` match a `box_mesh` whose middle edges are uniform at `h` and whose outer edges are wider (so the core box is identifiable); write `results/breakdown.npz` with `Kn_GLL`, `B_exp` set below threshold inside a cone and above outside, other indicators 0, `valid` all true, `retained` = the source-connected 1× mask computed in the fixture with `breakdown` functions, `mask`, `near_lip`, `time`.

```python
def test_load_rung_ok(rung_case):
    r = domainstudy.load_rung(rung_case)
    assert r["status"] == "ok"
    assert r["surfaces"]["2.0"]["volume_m3"] >= r["surfaces"]["1.0"]["volume_m3"] >= r["surfaces"]["0.5"]["volume_m3"]
    assert r["core_box"][1][0] == pytest.approx(<core x_hi in D>)

def test_load_rung_unsettled(rung_case):            # Review Focus 1
    (res/"breakdown.npz").unlink()
    write_json(res/"breakdown.failed.json", {"failures": ["retained volume fraction changed by 0.05"]})
    r = domainstudy.load_rung(rung_case)
    assert r["status"] == "unsettled" and "0.05" in r["reason"]

def test_load_rung_missing(tmp_path):               # Review Focus 2
    # resolved config only, no results/ → status "missing"

def test_load_rung_rejects_mask_drift(rung_case):    # Review Focus 4
    # rewrite breakdown.npz with one extra retained cell → WorkflowError

def test_capsule_summary_states(rung_case):          # Review Focus 3
    # interface.json fit {"status": "no_feasible_capsule", ...} → fit_status no_feasible_capsule, passes False, R None
    # fit fitted + screen status "no_passing_candidate" → numbers present, passes False
    # fit fitted + "passed" → passes True
```

- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement `load_rung`.**
- [ ] **Step 4: Run, expect PASS**, plus `pytest tests/unit/test_hybridplume_domainstudy.py -q`.
- [ ] **Step 5: Commit** — `git commit -m "feat(hybridplume): load domain-study rung from case results"`

---

### Task 5: Study report and CLI

**Files:**
- Modify: `plumetools/hybridplume/domainstudy.py`, `plumetools/hybridplume/__main__.py`
- Test: `tests/unit/test_hybridplume_domainstudy_io.py` (append)

**Interfaces:**
- Consumes: `load_rung`, `verdicts`, `workflow.write_json`.
- Produces:
  - `study(cases: list[Path], out: Path) -> dict` — loads every case, runs `verdicts`, writes `out/domain_study.json` (the verdicts doc plus `spec` path, `tolerance`, `factors`, `stations_over_D`), `out/domain_study.md`, and `out/domain_study.png` (z = 0 overlay: each ok rung's `boundary_xy` as a scatter in its own colour, core box outline, capsule section when fitted; skipped with `plot_status: "matplotlib unavailable"` in the JSON if import fails). Returns the JSON doc.
  - `domain_study.md` contains: one table row per rung — x_max, status, axial extent, max radius, volume (1×), clearance to core (2×), capsule status, return fraction, subsonic fraction, valid, converged, capsule passes — then a changes table (each metric vs next-larger rung, in %), then the line `Chosen domain: x_max = N D` or `Chosen domain: none`, and one line per non-ok rung with its `reason`.
  - CLI: `python -m plumetools.hybridplume domain-study CASE [CASE ...] --out DIR` → prints the JSON summary via `_print`; `WorkflowError` → exit 2 like other steps. Add as its own subparser branch (it takes several cases, unlike the generic steps).

- [ ] **Step 1: Write the failing tests**

```python
def test_study_writes_outputs(two_rung_cases, tmp_path):
    doc = domainstudy.study(two_rung_cases, tmp_path/"out")
    assert (tmp_path/"out/domain_study.json").is_file()
    md = (tmp_path/"out/domain_study.md").read_text()
    assert "Chosen domain:" in md and len(doc["rungs"]) == 2

def test_study_with_missing_rung(two_rung_cases, tmp_path):   # Review Focus 2
    shutil.rmtree(two_rung_cases[1]/"results")
    doc = domainstudy.study(two_rung_cases, tmp_path/"out")
    assert [r["status"] for r in doc["rungs"]] == ["ok", "missing"]
    assert doc["chosen"] is None

def test_cli_domain_study(two_rung_cases, tmp_path):
    rc = main(["domain-study", *map(str, two_rung_cases), "--out", str(tmp_path/"o")])
    assert rc == 0 and (tmp_path/"o/domain_study.json").is_file()
```

`two_rung_cases`: the Task 4 fixture built twice with overlays `domain.x_max_over_D` 10 and 16 (the mesh is the same small box; only the resolved x_max differs).

- [ ] **Step 2: Run, expect FAIL.**
- [ ] **Step 3: Implement `study` and the CLI branch.**
- [ ] **Step 4: Run the full unit suite, expect PASS:** `pytest tests/unit -q`
- [ ] **Step 5: Commit** — `git commit -m "feat(hybridplume): domain-study report and CLI"`

---

### Task 6: Pre-flight and submission (operational)

**Files:** none modified, except the spec/profile comments if the planner values differ from the Global Constraints.

- [ ] **Step 1: Set up all five rungs on the login node (Python only)**

```bash
cd /scratch/torresfa/OpenFOAM/rarefiedPlume
source cases/cai2012-kn-extension/hpc_env.sh
for x in 10 16 24 32 48; do
  bash cases/cai2012-hybrid/Allrun cases/cai2012-hybrid/runs/domstudy-kn1em3-x$x setup \
       --config cases/cai2012-hybrid/hybrid.yaml --profile domstudy-kn1em3-x$x -np 192
done
```

Expected for each: `results/setup.json` with `mesh.n_cells == 2944656`, `mesh.max_adjacent_cell_ratio <= 1.43`, a resources/memory estimate under `max_memory_GB` (1100), `source_discretisation.relative_area_error ≈ 0.006`, and a CFD end time ≈ x/10 × 2.83 ms. Any mismatch → stop and report.

- [ ] **Step 2: Check the `STAGES` path end-to-end without the solver:** `STAGES="setup:setup" PROFILE=domstudy-kn1em3-x10 bash cases/cai2012-hybrid/run3DCase.sh` on the login node; expect the banner, `setup` resuming ("set up from …"), exit 0.

- [ ] **Step 3: Ask the user to confirm submission** (five 192-core `longq` jobs). Do not submit without an explicit yes.

- [ ] **Step 4: Submit**

```bash
cd /scratch/torresfa/OpenFOAM/rarefiedPlume/cases/cai2012-hybrid
for x in 10 16 24 32 48; do
  msub -v PROFILE=domstudy-kn1em3-x$x,STAGES=setup:cfd:post-cfd run3DCase.sh
done
showq -u torresfa
```

Expected: five job IDs, listed by `showq`.

- [ ] **Step 5: Record job IDs** in `cases/cai2012-hybrid/README.md` under a new "Domain-size study" section (profiles, submit command, analysis command from Task 7); commit.

---

### Task 7: Analyse (after jobs finish)

- [ ] **Step 1: Run**

```bash
python -m plumetools.hybridplume domain-study \
  cases/cai2012-hybrid/runs/domstudy-kn1em3-x{10,16,24,32,48} \
  --out cases/cai2012-hybrid/runs/domain-study
```

Expected: `domain_study.{json,md,png}`; exit 0 even if some rungs are `unsettled`/`missing`.

- [ ] **Step 2: Report to the user** the chosen domain (or which decision-rule branch applies from the spec: core containment fails / no capsule passes / no convergence), with the `domain_study.md` tables. `prepare-dsmc` on the chosen rung is a separate follow-up, not part of this plan.

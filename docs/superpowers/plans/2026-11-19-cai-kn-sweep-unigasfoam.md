# Cai Knudsen sweep on uniGasFoam — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new case family, `cases/cai2012-unigas`, runs the seven-decade Cai
plume on a patched uniGasFoam (hybrid USP-DSMC). Its output is converted to
dsmcFoam field names, so the existing analysis, a post-run suitability verdict,
the scaling inventory and a dsmcFoam cross-check all run on it.

**Architecture:**
* A new package, `plumetools/cai2012ugf/`, reuses the frozen `plumetools.cai2012`
  and `plumetools.cai2012ext` for physics, mesh, weight, run window and analysis.
* It adds the uniGasFoam dictionaries, a lossless field converter, a post-run
  verdict, restart and leg control, scaling and report adapters.
* uniGasFoam itself gets one opt-in patch on a branch of the owner's fork.

**Tech Stack:** Python ≥ 3.9 (numpy, PyYAML, pytest), Bash runner scripts,
OpenFOAM v2512, uniGasFoam (C++, wmake), PBS.

**Spec:** `docs/superpowers/specs/2026-11-19-cai-kn-sweep-unigasfoam-design.md`
(approved 2026-11-19). Read it before any task; section numbers below (§n) refer
to it.

## Global Constraints

* **Do not edit frozen code.** That means anything under
  `plumetools/cai2012/` or `plumetools/cai2012ext/`, and any existing case
  family, manifest or results tree. Import them; wrap them; never copy a
  function out of them.
* **Physics, mesh and run window are inherited.** The generated `case.yaml`
  equals the kn-extension's for the same Kn in every key outside `meta.*`.
  Mesh dictionaries are byte-identical, and `nEquivalentParticles` is the
  kn-extension's weight to 1e-12 relative.
* **Cell weighting is off:** `cellWeightedSimulation false`;
  `cellWeightAdaptation false`.
* **Collision settings:** `collisionModel hybrid`, `bgkCollisionModel
  unifiedStochasticParticleSBGK`, `dsmcCollisionPartnerModel noTimeCounter`,
  `dsmcCollisionModel variableHardSphere`, `Tref 273`, `macroInterpolation
  true`, `theta 1e-3`.
* **Argon:** mass 6.63e-26 kg, diameter 4.17e-10 m, **omega 0.74** (never the
  tutorial's 0.81), 0 rotational and vibrational modes, alpha 1.0.
* **Adaptive settings:** `adaptiveSimulation true`, `subCellAdaptation true`,
  `timeStepAdaptation true`, `adaptationInterval 20`, `maxSubCellSizeMFPRatio
  0.5`, `maxTimeStepMCTRatio 0.2`, `maxCourantNumber 0.5`,
  `dsmcCellsTimeStepLimited true`, `maxSubCellLevels 10`.
* **Hybrid decomposition:** `localKnudsen`, `breakdownMax 0.05`, `theta 0.2`,
  `smoothingPasses 5`, `decompositionInterval 100`, `resetAtDecomposition on`,
  `resetAtDecompositionUntilTime` = `average_start_s`.
* **Run window:** `average_start_s` = 0.0501435 s; 4.5 transits of averaging;
  sampling writes every 3.536e-4 s (the case's `dsmc.write_interval_s`);
  startup writes every 1 domain transit.
* **Initial state:** `uniGasMeshFill` at 1e10 m⁻³, 300 K, at rest.
* **Every result carries its solver.** Fits never pool uniGasFoam with dsmcFoam.
* **No validation claims.** README and report text says "not validated", as the
  repository does everywhere.
* **uniGasFoam changes** go on branch `dsmc-time-step-limit` of
  `../uniGasFoam`. Patch keys default to the published behaviour. Do not push;
  the owner pushes.
* **Commits** end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  rarefiedPlume work stays on branch `unigas-kn-sweep`.
* **Test commands:** `python -m pytest -m "not needs_openfoam"` is the fast
  suite. `needs_openfoam` tests run in WSL with OpenFOAM v2512 sourced.

## Review Focus

1. **Resume after a walltime kill mid-write.** Some ranks hold a complete new
   write and have already deleted their previous cloud; others don't. uniGasFoam
   re-initialises from vacuum when a time has no `lagrangian/uniGas`. The
   resume plan must refuse with manual-restart instructions, never continue.
   *Test:* Task 7, `test_resume_refuses_time_without_cloud_on_every_rank`.
2. **Empty and near-empty cells.**
   * `rhoM_Ar = 0`; `uniGasRhoNMean_Ar` written as 0.001; the mean free path is
     0, inf or NaN.
   * The converter must keep them empty.
   * The verdict must exclude them from DSMC statistics, counted, without
     NaN-crashing.
   * *Tests:* Task 4, `test_empty_cells_stay_empty`; Task 5,
     `test_verdict_excludes_unoccupied_cells`.
3. **Converter run twice, or on frames from before the averaging start.**
   * A second run is a no-op.
   * Frames before `average_start_s` get only `rhoN`, never `*Mean`.
   * A differing existing target is refused.
   * *Tests:* Task 4, `test_second_run_is_noop`, `test_no_mean_fields_before_average_start`
     and `test_refuses_differing_target`.
4. **Edited inputs drifting from the kn-extension.** A `baseCase/case.yaml` or
   `study.yaml` edit must fail validation, naming the key.
   *Test:* Task 6, `test_template_difference_is_named`.
5. **Analysis run before the verdict exists, or a Kn ≥ 1 case with no USP
   cells.** The scripts show `post-run-verdict-pending` (never production), and
   the verdict handles a USP fraction of 0.
   *Tests:* Task 5, `test_context_pending_without_verdict` and
   `test_verdict_with_no_usp_cells`.

---

### Task 1: Build uniGasFoam on v2512 and add the opt-in patch (fork repo)

All paths in this task are in `../uniGasFoam` (relative to the rarefiedPlume
root).

**Files:**
- Modify: `src/lagrangian/uniGas/dynamicAdaptation/uniGasDynamicAdapter.H`. Add
  member `bool dsmcCellsTimeStepLimited_;` next to `cellWeightAdaptation_`.
- Modify: `src/lagrangian/uniGas/dynamicAdaptation/uniGasDynamicAdapter.C`:
  * the initialiser list (`maxSubCellLevels_(10)`, about line 49) and the
    member initialiser for the new bool (`false`);
  * the dictionary read (about line 202, beside `cellWeightAdaptation`);
  * `calculateTimeStep` (about line 347).
- Modify: `doc/GUIDE.md`, the adaptiveProperties key list.

**Interfaces:**
- Produces: `adaptiveProperties` keys `dsmcCellsTimeStepLimited` (bool,
  default `false`) and `maxSubCellLevels` (label ≥ 1, default 10). Tasks 3 and
  7 rely on these names. It also produces the fork commit hash that Task 6
  stores in `study.yaml` as `meta.unigas_commit`.

- [ ] **Step 1: Branch the fork.**
  `git -C ../uniGasFoam switch -c dsmc-time-step-limit`. Expected: `Switched to
  a new branch`.
- [ ] **Step 2: Build unpatched in WSL.**
  Run `source <OpenFOAM-v2512>/etc/bashrc && cd ../uniGasFoam && ./Allwmake -j`.
  Expected: `which uniGasFoam` prints a path under `$FOAM_USER_APPBIN`.
  If v2512 fails to compile, make the minimal fix and commit it separately as
  `fix: build on OpenFOAM v2512`, listing every change in the message.
- [ ] **Step 3: Record a baseline on the plume tutorial.**
  1. Copy `tutorials/uniGasFoam/plumeImpingement` to a scratch directory.
  2. Set `decompositionInterval 20` in `hybridDecompositionDict`, so DSMC cells
     appear early.
  3. Set `endTime 3e-5` in `controlDict`.
  4. Run `blockMesh`, `cp -r orig.0 0`, then `uniGasFoam > log.baseline`.

  Keep `grep -E "Time step|Number of particles" log.baseline` as the baseline.
- [ ] **Step 4: Implement both keys.**
  * The constructor reads them with `getOrDefault`.
  * `maxSubCellLevels < 1` is a `FatalIOError`.
  * In `calculateTimeStep`, the MCT condition becomes:

    ```cpp
    (dsmcCellsTimeStepLimited_ || cloud_.bgkCollisionModelName() != "unifiedStochasticParticleSBGK")
    && cloud_.cellCollModelId()[cell] != cloud_.bgkCollModelId()
    && instMaxTimeStepMCTRatio < timeStepMCTRatio_[cell]
    ```

  Document both keys in `doc/GUIDE.md`, giving each one's default and that the
  default reproduces published behaviour.
- [ ] **Step 5: Rebuild and check that the defaults are unchanged.**
  Rerun Step 3's case with neither key present, to `log.patched-default`.
  Expected: the `Time step` / `Number of particles` lines equal the baseline.
  If the run is not bit-deterministic, the `Time step` sequence alone must be
  equal.
- [ ] **Step 6: Check the time-step limit.**
  1. Rerun with `dsmcCellsTimeStepLimited true;` in `adaptiveProperties`.
  2. Check the last write with this command:

     ```bash
     T=$(foamListTimes -latestTime) python - <<'PY'
     import os, numpy as np
     from plumetools.foamio.volfield import read_internal_field as r
     t = os.environ["T"]
     dsmc = r(f"{t}/uniGasCollisionModelId") == 1
     print(t, dsmc.sum(), np.nanmax(r(f"{t}/timeStepMCTRatio_Ar")[dsmc]))
     PY
     ```

  Expected: the printed maximum is ≤ 0.2 × 1.05 (the field is time-averaged).
  Run the same check on `log.patched-default`'s case. If that one is also
  ≤ 0.2, record "patch not exercised by the tutorial" in the commit message;
  the Kn = 0.001 pilot (Task 10) exercises it.
- [ ] **Step 7: Check the sub-cell cap.**
  Rerun with `maxSubCellLevels 2;`. Expected: the maximum component of
  `uniGasSubCellLevels` at the last write is ≤ 2.
- [ ] **Step 8: Commit in the fork.**
  `git -C ../uniGasFoam commit -am "feat: opt-in collision-time limit for DSMC cells under USP; configurable sub-cell cap"`.
  Record the hash for Task 6.

### Task 2: `unigas.yaml` settings

The cai2012 `case.yaml` schema rejects new sections. The uniGasFoam settings
therefore live beside it in each case's `unigas.yaml`, as interface mode keeps
`interface.yaml`.

**Files:**
- Create: `plumetools/cai2012ugf/__init__.py`. Its docstring states the
  package's purpose and that it never edits the frozen packages.
- Create: `plumetools/cai2012ugf/settings.py`
- Test: `tests/unit/test_cai2012ugf_settings.py`

**Interfaces:**
- Produces:
  * `class UniGasConfigError(ValueError)`.
  * Frozen dataclasses:
    * `CollisionSettings(model, bgk_model, dsmc_partner_model, dsmc_model, macro_interpolation, theta)`
    * `AdaptiveSettings(adaptation_interval, max_sub_cell_size_mfp_ratio, max_time_step_mct_ratio, max_courant_number, dsmc_cells_time_step_limited, max_sub_cell_levels)`
    * `DecompositionSettings(breakdown_max, theta, smoothing_passes, decomposition_interval)`
    * `WriteSettings(startup_interval_transits, keep_lagrangian)`
    * `FieldSettings(mean_field, window_field)`
    * `VerdictSettings(max_effective_spacing_over_mfp: float | None, max_time_step_mct_ratio, min_particles_per_dsmc_cell: float | None, min_particles_per_sub_cell, occupancy_floor_n_over_n0)`
    * `UniGasSettings(cell_weighted, collision, adaptive, decomposition, writes, fields, verdict)`
  * `load_settings(case_dir: Path) -> UniGasSettings`
  * `settings_from_dict(data: dict, where: str) -> UniGasSettings`
  * The YAML keys equal the field names, one section per dataclass, with
    `cell_weighted` at top level.

- [ ] **Step 1: Write the failing tests.**

```python
def test_defaults_are_the_spec_values():
    s = settings_from_dict({}, "test")
    assert s.cell_weighted is False
    assert (s.collision.model, s.collision.bgk_model) == ("hybrid", "unifiedStochasticParticleSBGK")
    assert (s.collision.dsmc_partner_model, s.collision.dsmc_model) == ("noTimeCounter", "variableHardSphere")
    assert s.collision.macro_interpolation is True and s.collision.theta == 1e-3
    a = s.adaptive
    assert (a.adaptation_interval, a.max_sub_cell_size_mfp_ratio, a.max_time_step_mct_ratio,
            a.max_courant_number, a.dsmc_cells_time_step_limited, a.max_sub_cell_levels) == (20, 0.5, 0.2, 0.5, True, 10)
    d = s.decomposition
    assert (d.breakdown_max, d.theta, d.smoothing_passes, d.decomposition_interval) == (0.05, 0.2, 5, 100)
    assert (s.writes.startup_interval_transits, s.writes.keep_lagrangian) == (1.0, False)
    assert (s.fields.mean_field, s.fields.window_field) == ("Ar", "ArWindow")
    v = s.verdict
    assert v.max_effective_spacing_over_mfp is None and v.min_particles_per_dsmc_cell is None
    assert (v.max_time_step_mct_ratio, v.min_particles_per_sub_cell, v.occupancy_floor_n_over_n0) == (0.2, 2.0, 1e-3)

def test_cell_weighting_is_refused():           # decision 5
    with pytest.raises(UniGasConfigError, match="cell_weighted"):
        settings_from_dict({"cell_weighted": True}, "x")

def test_unknown_key_is_named():
    with pytest.raises(UniGasConfigError, match="adaptive.maxCourant"):
        settings_from_dict({"adaptive": {"maxCourant": 1}}, "x")

def test_non_usp_or_non_hybrid_is_refused():
    for bad in ({"collision": {"model": "dsmc"}}, {"collision": {"bgk_model": "stochasticParticleSBGK"}}):
        with pytest.raises(UniGasConfigError):
            settings_from_dict(bad, "x")

def test_time_step_limit_off_is_refused():      # decision 4
    with pytest.raises(UniGasConfigError, match="dsmc_cells_time_step_limited"):
        settings_from_dict({"adaptive": {"dsmc_cells_time_step_limited": False}}, "x")

def test_same_mean_and_window_field_is_refused():
    with pytest.raises(UniGasConfigError):
        settings_from_dict({"fields": {"window_field": "Ar"}}, "x")

def test_missing_file_names_its_path(tmp_path):
    with pytest.raises(UniGasConfigError, match="unigas.yaml"):
        load_settings(tmp_path)
```

- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_settings.py -q`. Expected:
  ImportError / FAIL.
- [ ] **Step 3: Implement `settings.py`.**
  Validation: positive numbers where physical; `max_sub_cell_levels ≥ 1`; and
  the refusals above, each naming the dotted key and citing the spec decision.
- [ ] **Step 4: Run the tests; they pass.** Same command. Expected: 7 passed.
- [ ] **Step 5: Commit.** `feat(cai2012ugf): load and validate unigas.yaml`

### Task 3: The uniGasFoam dictionaries

**Files:**
- Create: `plumetools/cai2012ugf/dictionaries.py`
- Test: `tests/unit/test_cai2012ugf_dictionaries.py`. Build `cfg`, `exit_state`
  and `run` from the kn-extension template at Kn = 0.001 via
  `plumetools.cai2012ext.direct_case.derive_config`, as
  `tests/unit/_cai2012.py` helpers do.

**Interfaces:**
- Consumes:
  * `UniGasSettings` (Task 2);
  * `cfg: CaiCaseConfig`, `exit_state: ExitState`, `run: RunSettings`
    (`plumetools.cai2012`);
  * the frozen `plumetools.cai2012.dictionaries.render_decompose_par_dict(cfg)`
    and `plumetools.foamio.primitives.write_lines(path, lines)`.
- Produces:
  * `render_uni_gas_properties(cfg, run, settings) -> str`
  * `render_boundaries_dict(cfg, exit_state, settings) -> str`
  * `render_field_properties_dict(cfg, run, settings) -> str`
  * `render_initialisation_dict(cfg) -> str`
  * `render_hybrid_decomposition_dict(run, settings) -> str`
  * `render_control_dict(cfg, run, settings, leg: str) -> str`, where `leg` is
    `"startup"` or `"sampling"`
  * `write_case_dictionaries(case_dir: Path, cfg, exit_state, run, settings) -> list[Path]`.
    It writes `constant/uniGasProperties` and these under `system/`:
    `boundariesDict`, `fieldPropertiesDict`, `uniGasInitialisationDict`,
    `hybridDecompositionDict`, `chemReactDict` (an empty `reactions ();`),
    `controlDict.startup`, `controlDict.sampling`, `controlDict` (a copy of
    `.startup`), `decomposeParDict`, `fvSchemes` and `fvSolution`. It also
    creates an empty `0/`.
  * Constant `LEGS = ("startup", "sampling")`.
  * Constant `GENERATOR = "plumetools.cai2012ugf.dictionaries"`, written in
    every header.

- [ ] **Step 1: Write the failing tests.**
  Read values back with a small `entry(text, key)` regex helper in the test
  file.

```python
def test_weight_is_the_kn_extension_weight(case):
    assert float(entry(props(case), "nEquivalentParticles")) == pytest.approx(case.run.n_equivalent_particles, rel=1e-12)
def test_cell_weighting_off_and_adaptive_keys(case):
    t = props(case)
    for key, value in [("cellWeightedSimulation", "false"), ("cellWeightAdaptation", "false"),
                       ("axisymmetricSimulation", "false"), ("subCellAdaptation", "true"),
                       ("timeStepAdaptation", "true"), ("adaptationInterval", "20"),
                       ("maxSubCellSizeMFPRatio", "0.5"), ("maxTimeStepMCTRatio", "0.2"),
                       ("maxCourantNumber", "0.5"), ("dsmcCellsTimeStepLimited", "true"),
                       ("maxSubCellLevels", "10")]:
        assert entry(t, key) == value
def test_argon_is_the_repository_vhs_model(case):
    t = props(case)
    assert float(entry(t, "omega")) == 0.74 and "0.81" not in t
    assert float(entry(t, "mass")) == 6.63e-26 and float(entry(t, "diameter")) == 4.17e-10
def test_collision_block(case):
    t = props(case)
    assert entry(t, "collisionModel") == "hybrid"
    assert entry(t, "bgkCollisionModel") == "unifiedStochasticParticleSBGK"
    assert entry(t, "Tref") == "273" and entry(t, "macroInterpolation") == "true"
def test_nozzle_injects_and_deletes_vacuum_only_deletes(case):
    t = render_boundaries_dict(case.cfg, case.exit_state, SETTINGS)
    assert patch_models(t, "nozzle") == {"uniGasDeletionPatch", "uniGasFreeStreamInflowPatch"}
    assert patch_models(t, "upstreamVacuum") == {"uniGasDeletionPatch"}
    assert patch_models(t, "vacuum") == {"uniGasDeletionPatch"}
    assert number_density(t, "Ar") == pytest.approx(case.exit_state.number_density_per_m3, rel=1e-12)
    assert inflow_velocity(t) == pytest.approx((case.exit_state.speed_m_per_s, 0.0, 0.0), rel=1e-12)
def test_averaging_window(case):
    t = render_field_properties_dict(case.cfg, case.run, SETTINGS)
    assert reset_until(t, "Ar") == pytest.approx(case.run.average_start_s)
    assert reset_until(t, "ArWindow") > case.run.end_time_s
    assert measure_mfp(t, "Ar") == "true"
def test_two_legs(case):
    s = render_control_dict(case.cfg, case.run, SETTINGS, "startup")
    p = render_control_dict(case.cfg, case.run, SETTINGS, "sampling")
    assert entry(s, "application") == entry(p, "application") == "uniGasFoam"
    assert float(entry(s, "endTime")) == pytest.approx(case.run.average_start_s)
    assert float(entry(s, "writeInterval")) == pytest.approx(case.run.domain_transit_s)
    assert float(entry(p, "endTime")) == pytest.approx(case.run.end_time_s)
    assert float(entry(p, "writeInterval")) == pytest.approx(case.run.write_interval_s)
    assert entry(p, "writeControl") == "adjustableRunTime" and entry(p, "purgeWrite") == "0"
    assert entry(p, "startFrom") == "latestTime" and entry(p, "writeFormat") == "ascii"
def test_decomposition(case):
    t = render_hybrid_decomposition_dict(case.run, SETTINGS)
    assert (entry(t, "breakdownMax"), entry(t, "smoothingPasses"), entry(t, "decompositionInterval")) == ("0.05", "5", "100")
    assert float(entry(t, "resetAtDecompositionUntilTime")) == pytest.approx(case.run.average_start_s)
def test_write_case_dictionaries_writes_the_set(case, tmp_path):
    names = {p.relative_to(tmp_path).as_posix() for p in write_case_dictionaries(tmp_path, case.cfg, case.exit_state, case.run, SETTINGS)}
    assert {"constant/uniGasProperties", "system/controlDict", "system/controlDict.startup",
            "system/controlDict.sampling", "system/chemReactDict", "system/decomposeParDict"} <= names
    assert (tmp_path / "system/controlDict").read_text() == (tmp_path / "system/controlDict.startup").read_text()
    assert (tmp_path / "0").is_dir()
```

  `ExitState` attribute names: use the ones `plumetools.cai2012.inflow.ExitState`
  actually defines. The names above are placeholders where they differ.
- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_dictionaries.py -q`.
- [ ] **Step 3: Implement.**
  * Model each layout on `tutorials/uniGasFoam/plumeImpingement` and
    `doc/GUIDE.md`.
  * Patch names come from `cfg.mesh.patch_names` (roles `nozzle`,
    `upstream_vacuum`, `outer`).
  * The inflow velocity is `U0 × cfg.geometry.plume_axis`. All four inflow
    temperatures are `T0`.
  * `measureErrors false`.
  * `ArWindow` gets `measureMeanFreePath false` and `resetAtOutputUntilTime`
    2 × `end_time_s`.
  * The startup `writeInterval` is `settings.writes.startup_interval_transits *
    run.domain_transit_s`.
  * `deltaT` is `run.delta_t_s`. It is only the initial step; uniGasFoam adapts
    it.
  * `timePrecision 8`; `nTerminalOutputs 10`; `runTimeModifiable yes`.
  * Comment `purgeWrite 0` in the file: every sampled frame is needed by
    `Allconverge.py`.
  * Copy `fvSchemes` / `fvSolution` from the tutorial's
    (`plumeImpingement/system`). Task 11 confirms them on v2512.
- [ ] **Step 4: Run the tests; they pass.** Expected: all pass.
- [ ] **Step 5: Commit.** `feat(cai2012ugf): render the uniGasFoam case dictionaries`

### Task 4: The converter

**Files:**
- Create: `plumetools/cai2012ugf/convert.py`, with a CLI:
  `python -m plumetools.cai2012ugf.convert <case_dir>`.
- Test: `tests/unit/test_cai2012ugf_convert.py`

**Interfaces:**
- Consumes:
  * `plumetools.foamio.volfield.read_internal_field(path) -> np.ndarray`;
  * `plumetools.foamio.primitives.write_lines`;
  * the frozen `plumetools.cai2012.post.read_case(case_dir, species_mass_kg, time=..., mesh_time=...)`
    (tests only);
  * `FieldSettings` (Task 2).
- Produces:
  * `class ConversionError(RuntimeError)`
  * `@dataclass ConversionResult(time: str, written: list[Path], kept: list[str], skipped: list[str])`
  * `convert_time(time_dir: Path, *, species_mass_kg: float, average_start_s: float, fields: FieldSettings, patches: dict) -> ConversionResult`
  * `convert_case(case_dir: Path) -> list[ConversionResult]`. It loads
    `case.yaml` and `unigas.yaml` and visits every reconstructed numeric time
    directory except `0`.
  * `MAPPING`: the table below, as data.

| target | class | dimensions | value |
|---|---|---|---|
| `rhoNMean` | volScalarField | `[0 -3 0 0 0 0 0]` | `rhoN_<f>` |
| `rhoMMean` | volScalarField | `[1 -3 0 0 0 0 0]` | `rhoM_<f>` |
| `momentumMean` | volVectorField | `[1 -2 -1 0 0 0 0]` | `rhoM_<f> * UMean_<f>` |
| `linearKEMean` | volScalarField | `[1 -1 -2 0 0 0 0]` | `0.5*rhoM*|U|^2 + 1.5*rhoN*k_B*translationalT_<f>`, with k_B = 1.380649e-23 |
| `dsmcRhoNMean` | volScalarField | `[0 -3 0 0 0 0 0]` | `uniGasRhoNMean_<f>`, set to 0 where `rhoN_<f> == 0` |
| `rhoN` | volScalarField | `[0 -3 0 0 0 0 0]` | `rhoN_<w>` |

`<f>` is `fields.mean_field` and `<w>` is `fields.window_field`. The `*Mean`
targets are written only when the time is ≥ `average_start_s` × (1 − 1e-9).
`rhoN` is written at every time that has `rhoN_<w>`.

- [ ] **Step 1: Write the failing tests.**
  Synthetic cases are built in `tmp_path`. Each has:
  * a copy of `cases/cai2012-kn-extension/baseCase/case.yaml` and an empty
    `unigas.yaml` (`{}`), so `convert_case` derives `average_start_s` =
    0.0501435 s;
  * a `C`, a `V`, and a `constant/polyMesh/boundary` with three patches;
  * uniGasFoam fields for 5 cells, one of them empty, at a startup time `0.03`
    and a sampled time `0.06`.

```python
def test_round_trip_is_exact(synthetic_case):
    convert_case(synthetic_case.dir)
    s = read_case(synthetic_case.dir, 6.63e-26, time=synthetic_case.sampled_time)
    occ = synthetic_case.rhoM > 0
    np.testing.assert_allclose(s.number_density, synthetic_case.rhoN, rtol=1e-12)
    np.testing.assert_allclose(s.velocity[occ], synthetic_case.U[occ], rtol=1e-12)
    np.testing.assert_allclose(s.temperature[occ], synthetic_case.T[occ], rtol=1e-12)
def test_empty_cells_stay_empty(synthetic_case):        # Review Focus 2
    convert_case(synthetic_case.dir)
    t = synthetic_case.dir / synthetic_case.sampled_time
    empty = synthetic_case.rhoM == 0
    assert (read_internal_field(t / "dsmcRhoNMean")[empty] == 0).all()   # source has 0.001 there
    assert (read_internal_field(t / "rhoMMean")[empty] == 0).all()
def test_no_mean_fields_before_average_start(synthetic_case):   # Review Focus 3
    convert_case(synthetic_case.dir)
    early = synthetic_case.dir / synthetic_case.startup_time
    assert (early / "rhoN").is_file() and not (early / "rhoNMean").exists()
def test_second_run_is_noop(synthetic_case):            # Review Focus 3
    convert_case(synthetic_case.dir)
    results = convert_case(synthetic_case.dir)
    assert all(not r.written for r in results)
def test_refuses_differing_target(synthetic_case):      # Review Focus 3
    convert_case(synthetic_case.dir)
    (synthetic_case.dir / synthetic_case.sampled_time / "rhoNMean").write_text("tampered")
    with pytest.raises(ConversionError, match="rhoNMean"):
        convert_case(synthetic_case.dir)
def test_never_touches_source_fields(synthetic_case):
    before = digests(synthetic_case.dir, "*_Ar")
    convert_case(synthetic_case.dir)
    assert digests(synthetic_case.dir, "*_Ar") == before
def test_header_names_converter_and_sources(synthetic_case):
    convert_case(synthetic_case.dir)
    text = (synthetic_case.dir / synthetic_case.sampled_time / "momentumMean").read_text()
    assert "plumetools.cai2012ugf.convert" in text and "rhoM_Ar" in text and "UMean_Ar" in text
def test_missing_window_field_is_reported_not_fatal(synthetic_case):
    (synthetic_case.dir / synthetic_case.sampled_time / "rhoN_ArWindow").unlink()
    r = {x.time: x for x in convert_case(synthetic_case.dir)}[synthetic_case.sampled_time]
    assert "rhoN" in r.skipped
```

- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_convert.py -q`.
- [ ] **Step 3: Implement.**
  * Use `read_internal_field` to read.
  * Write each field with `write_lines` and a standard FoamFile header. The
    header carries a `// converted by plumetools.cai2012ugf.convert from <sources>`
    line. Use `internalField nonuniform List<scalar|vector>` and a `calculated`
    zero `boundaryField` for each patch in `constant/polyMesh/boundary`, or
    `empty` for empty patches.
  * Build the patch dict the way `cases/cai2012-kn-extension/baseCase/runInflow.py`
    builds it.
  * Print at full precision (`repr`-equivalent, 17 significant digits), so the
    round trip is exact.
  * Compare an existing target by its parsed values (to 1e-15 relative) and
    header note. Equal means kept; different raises `ConversionError` naming
    the file.
- [ ] **Step 4: Run the tests; they pass.** Expected: 8 passed.
- [ ] **Step 5: Commit.** `feat(cai2012ugf): convert uniGasFoam fields to the dsmcFoam names plumetools reads`

### Task 5: The post-run verdict and the post context

**Files:**
- Create: `plumetools/cai2012ugf/verdict.py`, with a CLI:
  `python -m plumetools.cai2012ugf.verdict <case_dir>`, which writes
  `results/unigas_verdict.yaml`.
- Create: `plumetools/cai2012ugf/context.py`
- Test: `tests/unit/test_cai2012ugf_verdict.py`

**Interfaces:**
- Consumes:
  * `VerdictSettings` (Task 2);
  * `checks.max_cell_over_mfp` and `checks.min_particles_per_cell` from
    `case.yaml`, used when the verdict settings are `None`;
  * the frozen `plumetools.cai2012ext.direct_case.derive`,
    `plumetools.cai2012ext.interface_case.PostContext`, and
    `plumetools.cai2012ext.policy.source_limitations(cfg)`;
  * `plumetools.rarefaction.foam` per-cell extents (`cell_extents`, (n, 3) in
    m);
  * `results/stationarity.yaml`, written by `Allconverge.py`.
- Produces:
  * `effective_spacing_over_mfp(extents, sub_cell_levels, particles, mfp, *, min_particles_per_sub_cell=2.0) -> np.ndarray`
  * `assess_fields(fields: dict[str, np.ndarray], *, n0: float, settings: VerdictSettings, limits: dict, stationarity_verdict: str | None, known_limitations: list[str]) -> dict`
    (pure)
  * `assess_run(case_dir: Path, *, time: str | None = None) -> dict` (reads the
    files, calls `assess_fields`, writes the YAML)
  * `context.post_context(case_dir: Path) -> PostContext`
  * Verdict keys: `status` (`production-candidate` | `non-production`),
    `production_ready` (bool), `unresolved` (list[str]), plus `dsmc`, `usp`,
    `kn_gll`, `stationarity`, `known_limitations`, `time`.
  * The pending verdict when the file is absent: `{"status": "post-run-verdict-pending", "production_ready": False, "unresolved": ["post-run verdict not computed: run Allverdict.py"]}`.

The effective-spacing estimator is spec §7, made exact:

```text
S = sx*sy*sz ;  f = min(1, N / (p*S)) ** (1/3)         # p = min_particles_per_sub_cell
divisions_d = max(1, s_d * f)
effective   = max_d( extent_d / divisions_d ) / mfp     # NaN where mfp is not finite and > 0
```

The gates apply to DSMC cells (`uniGasCollisionModelId == 1`) that are occupied
(`uniGasRhoNMean > 0` after conversion, i.e. `dsmcRhoNMean > 0`):
* `effective ≤ limit`;
* `timeStepMCTRatio ≤ 0.2 × (1 + 1e-6)`;
* particles ≥ `min_particles_per_dsmc_cell`, counted only where
  `rhoN/n0 ≥ occupancy_floor_n_over_n0` (1e-3, Cai's lowest contour).

`production_ready` additionally requires a stationarity verdict of
`stable-window-demonstrated` and no known limitation. Statistics are reported as
max, p50, p95, p99 and the fraction over the limit, with counts of excluded
cells.

- [ ] **Step 1: Write the failing tests.**

```python
def test_dense_sub_cells_resolve():
    r = effective_spacing_over_mfp(np.array([[1.0, 1, 1]]), np.array([[5, 5, 5]]), np.array([1000.0]), np.array([1.0]))
    assert r[0] == pytest.approx(0.2)
def test_sparse_sub_cells_fall_back_toward_cell_size():
    r = effective_spacing_over_mfp(np.array([[1.0, 1, 1]]), np.array([[10, 10, 10]]), np.array([40.0]), np.array([0.1]))
    assert r[0] == pytest.approx(1 / (10 * 0.02 ** (1 / 3)) / 0.1, rel=1e-12)   # 3.684
def test_one_particle_means_cell_size():
    r = effective_spacing_over_mfp(np.array([[2.0, 1, 1]]), np.array([[10, 10, 10]]), np.array([1.0]), np.array([1.0]))
    assert r[0] == pytest.approx(2.0)
def test_verdict_excludes_unoccupied_cells():            # Review Focus 2
    v = assess_fields(fields_with(empty_dsmc_cell=True, mfp_in_empty=0.0), **COMMON)
    assert v["dsmc"]["excluded_unoccupied"] == 1 and np.isfinite(v["dsmc"]["effective_spacing_over_mfp"]["max"])
def test_verdict_with_no_usp_cells():                    # Review Focus 5
    v = assess_fields(fields_with(all_dsmc=True), **COMMON)
    assert v["usp"]["fraction"] == 0.0 and v["status"] in ("production-candidate", "non-production")
def test_each_gate_names_its_failure():
    for broken, phrase in [("spacing", "collision spacing"), ("mct", "time step"),
                           ("occupancy", "particles"), ("stationarity", "stationarity")]:
        v = assess_fields(fields_with(broken=broken), **COMMON)
        assert v["production_ready"] is False and any(phrase in u for u in v["unresolved"])
def test_known_limitation_blocks_production():
    v = assess_fields(fields_with(), **{**COMMON, "known_limitations": ["direct-exit inflow surface"]})
    assert v["status"] == "non-production"
def test_all_gates_pass():
    v = assess_fields(fields_with(), **COMMON)
    assert v["status"] == "production-candidate" and v["production_ready"] is True and v["unresolved"] == []
def test_context_pending_without_verdict(kn_case_dir):   # Review Focus 5
    c = post_context(kn_case_dir)
    assert c.verdict["status"] == "post-run-verdict-pending" and c.verdict["production_ready"] is False
    assert c.mode == "direct" and c.refusal is None
def test_context_reads_written_verdict(kn_case_dir):
    write_yaml(kn_case_dir / "results/unigas_verdict.yaml", {"status": "non-production", "production_ready": False, "unresolved": ["x"]})
    assert post_context(kn_case_dir).verdict["unresolved"] == ["x"]
```

  `kn_case_dir` is `tmp_path` holding a copy of
  `cases/cai2012-kn-extension/baseCase/case.yaml` with `exit.knudsen: 1.0`.
  Task 6 has not run yet, so do not depend on its generated family.
  `COMMON` is the `n0`, `settings`, `limits`, `stationarity_verdict` and
  `known_limitations` keyword arguments with passing values. `fields_with`
  builds a small consistent dict of arrays (keys `collision_model_id`,
  `sub_cell_levels`, `particles`, `mfp`, `extents`, `time_step_mct_ratio`,
  `kn_gll`, `number_density`, `centres`).
- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_verdict.py -q`.
- [ ] **Step 3: Implement.**
  * `assess_run` reads from the latest reconstructed time:
    `uniGasCollisionModelId`, `uniGasSubCellLevels`, `dsmcRhoNMean` (converted),
    `variableHardSphereMeanFreePath_Ar`, `timeStepMCTRatio_Ar`, `KnGLL`,
    `rhoNMean`, `C` and the mesh extents.
  * Known limitations come from `policy.source_limitations(cfg)`; Kn ≤ 0.001
    gives the direct-exit inflow note.
  * `post_context` returns the frozen `PostContext("direct", case.cfg,
    case.geom, case.exit_state, case.run, verdict,
    float(case.nozzle.area_error), None, None)`.
- [ ] **Step 4: Run the tests; they pass.** Expected: 10 passed.
- [ ] **Step 5: Commit.** `feat(cai2012ugf): post-run suitability verdict from DSMC cells`

### Task 6: Case family and idempotent generation

**Files:**
- Create: `cases/cai2012-unigas/study.yaml`. Copy the kn-extension's
  `study.yaml`, then change only:
  * `cases_dir: Cases`;
  * `manifest_name: manifest-unigas.yaml`;
  * `meta.study`, `meta.solver: "uniGasFoam hybrid USP-DSMC"`;
  * `meta.unigas_commit: <Task 1 hash>`;
  * `meta.unigas_branch: dsmc-time-step-limit`.

  Keep `meta.knudsen_suite: full_decades`.
- Create: `cases/cai2012-unigas/baseCase/case.yaml`, a copy of the
  kn-extension's with only its header comment and `meta.study` changed.
- Create: `cases/cai2012-unigas/baseCase/unigas.yaml`. Every key is explicit, at
  the Task 2 defaults, each commented with its spec section.
- Create: `cases/cai2012-unigas/baseCase/diagnostics.yaml`, a byte copy of the
  kn-extension's.
- Create: `plumetools/cai2012ugf/policy.py`, `plumetools/cai2012ugf/generation.py`.
- Create: `cases/cai2012-unigas/generate_cases.py`, `cases/cai2012-unigas/validate_cases.py`.
  They are thin CLIs; the logic is in the package.
- Create: `cases/cai2012-unigas/README.md`. It covers:
  * purpose and that nothing is validated;
  * commands;
  * the pending-verdict note;
  * the per-write-interval convergence series note (spec §5);
  * the cloud-restart warning (spec §6).
- Test: `tests/unit/test_cai2012ugf_generation.py`

**Interfaces:**
- Consumes:
  * the frozen `plumetools.cai2012.study.load_study`, `apply_case_parameters`;
  * `plumetools.cai2012ext.policy.flatten_config`, `check_study`;
  * `plumetools.cai2012ext.generation.stage`, `compare`, `decide`, `apply`,
    `build_entry`, `write_manifest`, `has_solver_output`;
  * `plumetools.cai2012ext.direct_case.derive_config`;
  * `plumetools.cai2012.mesh.write_mesh_setup`.
- Produces:
  * `KN_EXTENSION = Path("cases/cai2012-kn-extension")`, relative to the repo
    root
  * `policy.template_differences(ugf_cfg, ext_cfg) -> list[str]`: dotted keys
    that differ, ignoring `meta.*`
  * `policy.study_differences(ugf_study, ext_study) -> list[str]`: compares
    `kn_cases` (name, Kn, enabled, overrides), `particle_levels`,
    `run_time_multiplier` and `output_frequency_multiplier`
  * `generation.generate(family_root: Path, *, dry_run=False, refresh=False, regenerate=False, overwrite=False) -> list[CaseState]`.
    It has the same semantics and statuses as the kn-extension's
    `generate_cases.py`, and its staleness comparison covers `case.yaml` **and**
    `unigas.yaml`.
  * `generation.stage_pilot(family_root, kn_name: str) -> Path`: writes
    `Pilot/<case>` with a `PILOT` marker file; never in the manifest.
  * `generation.validate(family_root) -> list[str]`, where empty means pass.
    It returns failures for template or study drift, a mesh-dictionary byte
    mismatch, a weight mismatch, or `unigas.yaml` invalid in any case.

- [ ] **Step 1: Write the failing tests.**

```python
def test_template_matches_kn_extension_outside_meta():
    assert template_differences(load_case_config(UGF / "baseCase"), load_case_config(EXT / "baseCase")) == []
def test_template_difference_is_named(tmp_path):             # Review Focus 4
    edited = copy_base(tmp_path, set_key=("dsmc.courant_target", 0.3))
    assert template_differences(load_case_config(edited), load_case_config(EXT / "baseCase")) == ["dsmc.courant_target"]
def test_study_matches_kn_extension():
    assert study_differences(load_study(UGF / "study.yaml"), load_study(EXT / "study.yaml")) == []
def test_seven_cases_with_the_suite_names(generated):
    assert sorted(p.name for p in (generated / "Cases").iterdir()) == sorted(
        f"{k}_np1x" for k in ("Kn100", "Kn10", "Kn1", "Kn0p1", "Kn0p01", "Kn0p001", "Kn0p0001"))
def test_each_case_has_unigas_yaml(generated):
    for case in (generated / "Cases").iterdir():
        assert load_settings(case).cell_weighted is False
def test_mesh_dictionaries_byte_identical_and_weight_equal(generated, kn_extension_staged):
    assert validate(generated) == []          # includes both comparisons, all seven Kn
def test_second_generation_keeps_everything(generated):
    states = generate(generated)
    assert {s.status for s in states} == {"identical"}
def test_unigas_yaml_edit_is_stale(generated):
    edit(generated / "baseCase/unigas.yaml", "max_courant_number: 0.5", "max_courant_number: 0.4")
    assert {s.status for s in generate(generated, dry_run=True)} == {"stale"}
def test_case_with_results_is_never_replaced(generated):
    case = generated / "Cases/Kn1_np1x"; (case / "0.06").mkdir(); (case / "0.06/rhoN_Ar").write_text("x")
    edit(generated / "baseCase/unigas.yaml", "theta: 0.2", "theta: 0.3")
    states = {s.name: s for s in generate(generated, regenerate=True)}
    assert states["Kn1_np1x"].status == "conflict"               # refused, as the kn-extension does
    assert (case / "0.06/rhoN_Ar").read_text() == "x"
def test_pilot_is_outside_the_manifest(generated):
    pilot = stage_pilot(generated, "Kn0p001")
    assert (pilot / "PILOT").is_file() and "Pilot" not in (generated / "manifest-unigas.yaml").read_text()
```

  `generated` copies `cases/cai2012-unigas` to `tmp_path` and runs `generate`.
  `kn_extension_staged` stages the kn-extension's seven cases in `tmp_path`
  with the frozen generator.
- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_generation.py -q`.
- [ ] **Step 3: Implement.**
  * Delegate to the frozen generation functions wherever they accept this
    family's root.
  * Where one hard-codes the kn-extension (template inheritance against
    cai2012-health, `INPUT_FILES`), wrap it rather than copy it.
  * If a needed frozen function cannot be wrapped, stop and report back with
    the function and the reason (spec §10's rule applies here too).
  * The manifest's entries add `solver: uniGasFoam` and `unigas_commit`.
  * Validation stages the kn-extension cases into a temporary directory. It then
    compares `blockMeshDict`, `topoSetDict`, `createPatchDict`,
    `meshQualityDict` and `decomposeParDict` byte for byte, and
    `run.n_equivalent_particles` to 1e-12 relative.
- [ ] **Step 4: Run the tests; they pass.** Expected: all pass. Then run
  `python cases/cai2012-unigas/generate_cases.py --dry-run`. Expected: seven
  `new` cases; Kn = 0.001 and 1e-4 carry the inflow-surface limitation.
- [ ] **Step 5: Commit.** `feat(cai2012-unigas): generate the seven-Kn uniGasFoam family`

### Task 7: Mesh, run legs, restart safety and the PBS jobs

**Files:**
- Create: `plumetools/cai2012ugf/hpc.py`, with a CLI:
  `python -m plumetools.cai2012ugf.hpc {resume-plan,finished,record,check-build} <case_dir> [...]`.
- Create, in `cases/cai2012-unigas/baseCase/`:
  * `Allmesh`: the kn-extension `Allmesh` with dictionary generation through
    `plumetools.cai2012ugf.dictionaries.write_case_dictionaries`. The mesh steps
    are unchanged.
  * `Allrun`
  * `Allclean`
  * `runMesh.sh`, `runCase.sh`
- Create: `cases/cai2012-unigas/runMeshAll.sh`, the kn-extension's with this
  family's generator, validator and manifest names.
- Create: `cases/cai2012-unigas/runBuild.sh`. It is a PBS job that:
  1. sources `../cai2012-kn-extension/hpc_env.sh`;
  2. checks out `meta.unigas_commit` in `$UNIGAS_DIR` (default
     `/home2/torresfa/scratch/OpenFOAM/uniGasFoam`);
  3. runs `./Allwmake`;
  4. writes the commit to `$FOAM_USER_APPBIN/uniGasFoam.commit`.
- Test: `tests/unit/test_cai2012ugf_hpc.py`

**Interfaces:**
- Consumes:
  * the frozen `plumetools.foamio.checkpoints.plan_restart(case_dir, n_subdomains) -> RestartPlan`
    and `RestartPlan.shell()`;
  * `run.average_start_s`, `run.end_time_s`;
  * `settings.writes.keep_lagrangian`.
- Produces:
  * `leg_for(time_s: float | None, average_start_s: float, end_time_s: float) -> str`,
    returning `"startup"`, `"sampling"` or `"done"`. Use a relative tolerance
    of 1e-9 at each boundary.
  * `cloud_complete(case_dir: Path, time: str) -> bool`. It requires
    `lagrangian/uniGas` at `time` on every `processor*`, with one file set, and
    `uniform/time` present. Serial: in the case directory.
  * `resume_plan(case_dir: Path, n_subdomains: int) -> UgfRestartPlan`, with
    fields `base: RestartPlan`, `leg: str`, `refusal: str | None`, and
    `shell()` adding `LEG=` and `REFUSAL=0/1`.
  * `finished(case_dir) -> bool`
  * `record(case_dir: Path, seconds: float) -> Path`. It appends to
    `<family>/results/cost-model.yaml` under a file lock, with keys `case`,
    `path`, `Kn`, `seconds`, `cells`, `solver: uniGasFoam`, `unigas_commit`,
    `production_status`.
  * `check_build(case_dir) -> str | None`: an error if
    `$FOAM_USER_APPBIN/uniGasFoam.commit` ≠ the manifest's `unigas_commit`.
  * `write_provenance(case_dir: Path) -> Path`. It writes
    `results/provenance-unigas.yaml` (spec §11) with:
    * `solver`;
    * `unigas_commit`, the built commit, and the patch keys' values;
    * the OpenFOAM version (`$WM_PROJECT_VERSION`);
    * the plumetools git provenance
      (`plumetools.cai2012ext.provenance.git_provenance`);
    * the converter module and its version;
    * every `uniGasProperties` and `hybridDecompositionDict` value as written;
    * the initial-state protocol.

    `Allrun` calls it before the first leg, and `Allverdict.py` adds the pre-
    and post-run verdicts.

**`Allrun` flow:**
1. `check-build`.
2. Run `eval "$(python -m plumetools.cai2012ugf.hpc resume-plan . --np $NP)"`.
   Refuse when `REFUSAL=1`, printing the plan's text.
3. Decompose or reconstruct as `RestartPlan` says (as in the kn-extension
   `Allrun`). There is no initialiser step: uniGasFoam fills the mesh itself on a
   fresh start.
4. Loop while `LEG != done`:
   1. `cp system/controlDict.$LEG system/controlDict`;
   2. `mpirun -np $NP uniGasFoam -parallel [-keep-lagrangian] >> log.uniGasFoam`;
   3. recompute the plan.
5. `reconstructPar -newTimes`.
6. `python -m plumetools.cai2012ugf.convert .`

`--serial` and `--stop-at-transits X` exist for the smoke test and the pilot.
`--stop-at-transits X` caps every leg's `endTime` at X × `domain_transit_s`
(via `foamDictionary` on the copied `controlDict`). The loop then also ends
once the newest time is ≥ that cap, so a pilot stopped inside the startup leg
does not loop.
`runCase.sh` keeps the kn-extension's PBS header and flow, using
`cai2012ugf.hpc` in place of `cai2012ext.hpc`.

- [ ] **Step 1: Write the failing tests.**

```python
def test_legs():
    assert leg_for(None, 0.05, 0.063) == "startup"
    assert leg_for(0.049, 0.05, 0.063) == "startup"
    assert leg_for(0.05 * (1 - 1e-12), 0.05, 0.063) == "sampling"
    assert leg_for(0.063, 0.05, 0.063) == "done"
def test_resume_refuses_time_without_cloud_on_every_rank(tmp_path):    # Review Focus 1
    make_ranks(tmp_path, n=2, times=["0.04", "0.05"], cloud_at={"0.05": [0]}, previous_cloud_deleted=True)
    plan = resume_plan(tmp_path, 2)
    assert plan.refusal and "lagrangian/uniGas" in plan.refusal and "manual restart" in plan.refusal
def test_resume_uses_newest_time_with_complete_clouds(tmp_path):
    make_ranks(tmp_path, n=2, times=["0.04", "0.05"], cloud_at={"0.05": [0, 1]})
    plan = resume_plan(tmp_path, 2)
    assert plan.refusal is None and plan.base.resume == "0.05" and plan.leg == "sampling"
def test_fresh_case_starts_in_startup(tmp_path):
    assert resume_plan(make_fresh(tmp_path), 2).leg == "startup"
def test_check_build_mismatch_is_reported(tmp_path, monkeypatch):
    monkeypatch.setenv("FOAM_USER_APPBIN", str(write_commit(tmp_path, "abc")))
    assert "abc" in check_build(case_with_commit(tmp_path, "def"))
def test_record_appends_with_solver(tmp_path):
    path = record(make_case(tmp_path), 123.0)
    assert yaml.safe_load(path.read_text())["measurements"][-1]["solver"] == "uniGasFoam"
def test_provenance_records_solver_commit_and_patch_keys(tmp_path, monkeypatch):
    monkeypatch.setenv("WM_PROJECT_VERSION", "v2512")
    p = yaml.safe_load(write_provenance(make_case(tmp_path, with_dictionaries=True)).read_text())
    assert p["solver"] == "uniGasFoam hybrid USP-DSMC" and p["openfoam_version"] == "v2512"
    assert p["uniGasProperties"]["dsmcCellsTimeStepLimited"] == "true"
    assert p["uniGasProperties"]["cellWeightedSimulation"] == "false"
```

  `make_ranks` builds `processor*/<t>/uniform/time` and
  `processor*/<t>/lagrangian/uniGas/{positions,...}` with the frozen
  checkpoints module's divider convention.
- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_hpc.py -q`.
- [ ] **Step 3: Implement `hpc.py` and the scripts.**
  The refusal text names the time, the ranks lacking a cloud, and the manual
  restart: remove the incomplete time from every `processor*` and resubmit. The
  newest earlier write with a cloud on every rank is the restart point; if none
  exists, the case restarts from scratch.
- [ ] **Step 4: Run the tests; they pass.** Also run `bash -n` on every new
  script. Expected: no syntax errors.
- [ ] **Step 5: Commit.** `feat(cai2012-unigas): two-leg runner with cloud-safe restarts and PBS jobs`

### Task 8: The post-processing chain

**Files:**
- Create: `cases/cai2012-unigas/baseCase/Allcompare.py`, `Allconverge.py`. These
  are copies of the kn-extension's in which **only** the import line
  `from plumetools.cai2012ext import figures, interface_case` becomes two lines:
  `from plumetools.cai2012ext import figures` and
  `from plumetools.cai2012ugf import context as interface_case`.
- Create: `cases/cai2012-unigas/baseCase/Allpost`, `Alldiagnose.py`,
  `Allaudit.py`, `postProcess.py`, `runPost.sh`, `runViz.sh`. These are
  byte-identical copies of the kn-extension's.
- Create: `cases/cai2012-unigas/baseCase/Allverdict.py`, which calls
  `verdict.assess_run`.
- Create: `cases/cai2012-unigas/baseCase/Allpostugf`. It is the per-case
  sequence: `convert`, `Allpost`, `Alldiagnose.py`, `Allaudit.py`,
  `Allconverge.py`, `Allverdict.py`, `Allcompare.py`. `runPost.sh` calls it.
- Create: `cases/cai2012-unigas/AllpostCases`, `AllreconstructCases`,
  `AllcleanCases`, `AllmeshCases`. These are the kn-extension's, with the
  manifest name and the per-case post command changed.
- Modify: `plumetools/cai2012ugf/generation.py` `validate`, adding the copy
  check below.
- Test: `tests/unit/test_cai2012ugf_postchain.py`

**Interfaces:**
- Consumes: Task 4's `convert_case`, Task 5's `assess_run` / `post_context`.
- Produces: `generation.copy_differences(family_root) -> list[str]`, which is
  empty when every copied script matches its original. For `Allcompare.py` and
  `Allconverge.py` the allowed difference is exactly the import line replaced
  as above.

- [ ] **Step 1: Write the failing tests.**

```python
def test_copies_match_originals_except_the_context_import():
    assert copy_differences(UGF) == []
def test_a_stray_edit_in_a_copy_is_named(tmp_path):
    fam = copy_family(tmp_path); append(fam / "baseCase/Allcompare.py", "# edit\n")
    assert any("Allcompare.py" in d for d in copy_differences(fam))
def test_post_order():
    steps = read_steps(UGF / "baseCase/Allpostugf")
    assert steps.index("convert") < steps.index("Allpost") < steps.index("Allconverge.py") < steps.index("Allverdict.py") < steps.index("Allcompare.py")
```

- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_postchain.py -q`.
- [ ] **Step 3: Create the copies and scripts.** Use `cp`, then the single edit.
- [ ] **Step 4: Run the tests; they pass.** Also run `bash -n` on the shell
  scripts. Expected: pass.
- [ ] **Step 5: Commit.** `feat(cai2012-unigas): post chain on converted fields with the post-run verdict`

### Task 9: Scaling inventory as a separate solver stratum

**Files:**
- Create: `plumetools/cai2012ugf/scaling.py`, with
  `main(argv=None, *, study_root=None) -> int` and the subcommands `inventory`,
  `fit` and `report`, in the style of the frozen `scaling_workflow.main`.
- Create: `cases/cai2012-unigas/Allscale.py` (a CLI) and
  `cases/cai2012-unigas/scaling.yaml`, a copy of the kn-extension's.
- Test: `tests/unit/test_cai2012ugf_scaling.py`

**Interfaces:**
- Consumes:
  * the frozen `plumetools.cai2012ext.scaling_inventory.inventory(study_root)`
    and `load_case(row)`;
  * `scaling_workflow.load_config`, `load_available`, `require_loaded`,
    `run_assessment`, `write_summary`;
  * each case's `results/unigas_verdict.yaml`.
- Produces: `inventory(study_root: Path) -> dict`. This is the frozen inventory
  with each row re-labelled:
  * `row["solver"] = "uniGasFoam"`;
  * `row["stratum"] = f"{stratum}+uniGasFoam"` (unchanged when `None`);
  * `row["suitability"] = {"suitable", "status", "reasons", "criteria": "cai2012ugf post-run verdict (spec §7)"}`,
    where `reasons` are the verdict's `unresolved` plus the row's
    `known_source_limitations`;
  * `row["suitable"] = verdict.production_ready and row["compatible"] and row["fit_ready"]`;
  * `row["warnings"]` with the frozen pre-run feasibility reasons removed and
    the new reasons added.

  Outputs go under `cases/cai2012-unigas/results/scaling/`.

- [ ] **Step 1: Write the failing tests.**

```python
def test_rows_carry_the_solver_and_a_separate_stratum(family_with_two_results):
    rows = inventory(family_with_two_results)["cases"]
    assert {r["solver"] for r in rows} == {"uniGasFoam"}
    assert all(r["stratum"].endswith("+uniGasFoam") for r in rows if r["stratum"])
def test_suitability_comes_from_the_post_run_verdict(family_with_two_results):
    rows = {r["name"]: r for r in inventory(family_with_two_results)["cases"]}
    assert rows["Kn1_np1x"]["suitable"] is True          # verdict production_ready, fit-ready fields
    assert rows["Kn0p001_np1x"]["suitable"] is False     # inflow-surface limitation
    assert not any("Pre-run" in w for w in rows["Kn1_np1x"]["warnings"])
def test_strata_never_mix_solvers():
    ugf = {"stratum": "abc+uniGasFoam"}; dsmc = {"stratum": "abc"}
    assert ugf["stratum"] != dsmc["stratum"]
```

  `family_with_two_results` is a generated family in `tmp_path` with synthetic
  `C`, `V`, `rhoNMean` and the verdict, diagnostics, occupancy and stationarity
  YAMLs for two cases.
- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_scaling.py -q`.
- [ ] **Step 3: Implement.** If the frozen `inventory` refuses this family's
  root (registration paths, schema), stop and report back with the error. Do not
  copy the inventory code (spec §10).
- [ ] **Step 4: Run the tests; they pass.**
- [ ] **Step 5: Commit.** `feat(cai2012ugf): scaling inventory as a uniGasFoam stratum`

### Task 10: Pilot and cross-check reports

**Files:**
- Create: `plumetools/cai2012ugf/pilot.py`, with a CLI:
  `python -m plumetools.cai2012ugf.pilot <case_dir>`, which writes
  `results/pilot_report.yaml` and `.md`.
- Create: `plumetools/cai2012ugf/crosscheck.py`, with a CLI:
  `python -m plumetools.cai2012ugf.crosscheck`, which writes
  `cases/cai2012-unigas/results/crosscheck.{yaml,md}`.
- Create: `cases/cai2012-unigas/runPilot.sh`. It is a PBS job that runs, in
  order:
  1. `generate_cases.py --pilot Kn0p001`;
  2. `Allmesh`;
  3. `Allrun --stop-at-transits 2`;
  4. `reconstructPar -latestTime`;
  5. `convert`;
  6. `pilot`.
- Test: `tests/unit/test_cai2012ugf_reports.py`

**Interfaces:**
- Consumes: the uniGasFoam log line formats:
  * `Time = <t>`
  * `    Time step                       = <dt>`
  * `    Number of particles             = <n>`
  * `    Maximum imbalance               = <pct>`
  * `ExecutionTime = <s> s  ClockTime = <s> s`

  It also consumes Task 5's `effective_spacing_over_mfp` and field readers, and
  each family's `results/comparison.yaml`.
- Produces:
  * `parse_log(path: Path) -> dict[str, np.ndarray]`, with keys `time`, `dt`,
    `particles`, `imbalance`, `clock`
  * `projected_hours(t_now: float, end: float, recent_clock_per_sim_s: float, elapsed_s: float) -> float`
  * `report(case_dir: Path) -> dict`, with:
    * the Δt history summary and which cells set it (DSMC MCT, or Courant);
    * DSMC/USP fractions and the x/D extents of each;
    * the effective spacing distribution;
    * the count of cells at `max_sub_cell_levels`;
    * the particle count;
    * the clock per step;
    * imbalance (max, p95);
    * bytes in the newest write across ranks;
    * `projected_wall_hours`
  * `crosscheck.report(ugf_root: Path, dsmc_root: Path, cases=("Kn0p1_np1x", "Kn1_np1x")) -> dict`.
    It puts every numeric metric in both `comparison.yaml` files side by side,
    with the difference and a `numerics` note (adaptive Δt at Co 0.5 on
    max(u, v0) with 0.2 τ, against fixed Co 0.2 on U0 + 3σ). A missing side is
    reported, never fabricated.

`projected_wall_hours` is the clock per simulated second over the last 10 % of
the log, times (`end_time_s` − t_now), plus elapsed.

- [ ] **Step 1: Write the failing tests.**
  Import `from plumetools.cai2012ugf.crosscheck import report as crosscheck_report`.

```python
def test_parse_log_reads_the_solver_lines(tmp_path):
    log = write_log(tmp_path, steps=[(1e-7, 1e-7, 1000, 3.0, 10.0), (2e-7, 1e-7, 1100, 4.0, 20.0)])
    out = parse_log(log)
    assert list(out["particles"]) == [1000, 1100] and out["clock"][-1] == 20.0
def test_projection_uses_recent_rate():
    assert projected_hours(t_now=0.01, end=0.02, recent_clock_per_sim_s=3.6e6, elapsed_s=0) == pytest.approx(10.0)
def test_crosscheck_pairs_metrics_and_reports_missing(tmp_path):
    ugf, dsmc = two_families(tmp_path, ugf_cases=["Kn1_np1x"], dsmc_cases=["Kn1_np1x", "Kn0p1_np1x"])
    r = crosscheck_report(ugf, dsmc)
    assert r["cases"]["Kn1_np1x"]["metrics"] and r["cases"]["Kn0p1_np1x"]["missing"] == ["uniGasFoam"]
    assert "Co 0.5" in r["numerics"]
```

- [ ] **Step 2: Run the tests; they fail.**
  `python -m pytest tests/unit/test_cai2012ugf_reports.py -q`.
- [ ] **Step 3: Implement.** Add `--pilot <KnName>` to
  `cases/cai2012-unigas/generate_cases.py`; it calls `generation.stage_pilot`.
- [ ] **Step 4: Run the tests; they pass.** Also `bash -n runPilot.sh`.
- [ ] **Step 5: Commit.** `feat(cai2012ugf): pilot report and dsmcFoam cross-check`

### Task 11: Local end-to-end smoke and documentation

**Files:**
- Create: `tests/openfoam/test_cai2012ugf_smoke.py`, marked `needs_openfoam`.
  Skip it unless `uniGasFoam` is on PATH.
- Modify: `plumetools/cai2012ugf/dictionaries.py`, only if the smoke shows
  uniGasFoam needs `0/` files, `fvSchemes` entries or `chemReactDict` content
  that Task 3 did not write. Record each finding in the commit message.
- Create: `docs/tutorials/cai2012-unigas.md`. It covers the build, the pilot,
  go/no-go, production order, post, cross-check and scaling, with the commands.
  It links the spec.
- Modify:
  * `mkdocs.yml` (nav);
  * `README.md` (a case-families row);
  * `docs/reference/case-catalog.md` (a row);
  * `docs/development/known-gaps.md` (the entries below).

  The known-gaps entries are:
  * the unrun pilot;
  * the patch not upstream;
  * the cloud-restart manual case;
  * the convergence series semantics;
  * pending-verdict banners in `stationarity.yaml`.

**Interfaces:**
- Consumes: everything above.

- [ ] **Step 1: Write the smoke test.** It runs these steps:
  1. Stage `Kn0p01_np1x` into `tmp_path` with a mesh override coarsening it to
     about 60,000 cells (`core_cell_size_m: 0.04`, `outer_x_cells: 10`,
     `outer_lateral_cells: 10`, `outer_expansion: 1.0`).
  2. Run `./Allmesh`, then `./Allrun --serial --stop-at-transits 0.05`. In the
     test only, set `decompositionInterval` to 20, so both cell types appear.
  3. Run `Allpostugf`.

  Then assert:

```python
assert (case / t / "rhoNMean").is_file() and (case / t / "momentumMean").is_file()
assert set(np.unique(read_internal_field(case / t / "uniGasCollisionModelId"))) == {0.0, 1.0}
occ = yaml.safe_load((case / "results/occupancy.yaml").read_text())
assert occ["weight_check"]["consistent"] is True                 # cell weighting off (decision 5)
v = yaml.safe_load((case / "results/unigas_verdict.yaml").read_text())
assert v["status"] in ("production-candidate", "non-production")
assert (case / "results/comparison.yaml").is_file()
```

  `t` must be ≥ `average_start_s`. In the test only, override `average_start_s`
  through `case.yaml`'s `dsmc.average_start_s` to 0.4 of the stop time, so
  `*Mean` fields exist.
- [ ] **Step 2: Run it in WSL.**
  `python -m pytest tests/openfoam/test_cai2012ugf_smoke.py -m needs_openfoam -q`.
  Expected: PASS. Fix the dictionaries (Task 3's module) for every uniGasFoam
  startup error. Rerun until it passes.
- [ ] **Step 3: Run the fast suite.** `python -m pytest -m "not needs_openfoam" -q`.
  Expected: everything passes, apart from failures listed in
  `docs/development/known-gaps.md` (the ffmpeg-dependent viz test).
- [ ] **Step 4: Write the docs.** Run `.venv-docs/bin/python -m mkdocs build --strict`.
  Expected: build succeeds.
- [ ] **Step 5: Commit.** `test,docs(cai2012-unigas): local end-to-end smoke and the family tutorial`

---

## After the plan (owner actions, not tasks)

1. Push the fork branch. On the cluster, `qsub cases/cai2012-unigas/runBuild.sh`.
2. `qsub cases/cai2012-unigas/runPilot.sh`, then read
   `Pilot/Kn0p001_np1x/results/pilot_report.md`. Go/no-go is the owner's
   decision (spec §8).
3. On go: `generate_cases.py`, `qsub runMeshAll.sh`, then `runCase.sh` for
   Kn = 0.1 and 1 first, then the rest. After each, run `runPost.sh`.
4. `python -m plumetools.cai2012ugf.crosscheck` and `python Allscale.py inventory | fit | report`.

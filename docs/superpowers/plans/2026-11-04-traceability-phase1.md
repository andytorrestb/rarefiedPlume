# Traceability and a learning path — Phase 1 implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Phase 1 of the approved design: a checked reading trace and a per-case run trace for the `cai2012` and `nozzleplume` (r-4d) families, the Learn section with two teaching cases, and the family guides moved to `docs/families/`, all without changing any existing output.

**Architecture:** A new reporting package, `plumetools/trace/`, holds one trace map per family. The maps are checked against what the family's own writers produce, rendered into committed pages under `docs/trace/`, and read by `python -m plumetools.trace record`, which writes the `derive` lines of `log.trace`. A canonical bash block, embedded verbatim in each traced runner script, writes the `decide`/`step`/header/end lines. Teaching cases are plain case directories under `cases/learn/` whose scripts are byte-identical to their family's `baseCase`.

**Tech Stack:** Python ≥ 3.9 (numpy, scipy, PyYAML only), bash, pytest, OpenFOAM v2512 (only for Task 15), MkDocs Material (pinned in `docs/requirements.txt`).

**Spec:** [`docs/superpowers/specs/2026-11-04-traceability-and-learning-design.md`](../specs/2026-11-04-traceability-and-learning-design.md) (branch `docs/traceability-learning-spec`, commit `12fc011`). Phases 2 and 3 get their own plans once this phase's interfaces exist.

## Global Constraints

- **Invariant 1:** every generated dictionary and `0/` field of every Phase 1 `baseCase` is byte-identical before and after.
- **Invariant 2:** every existing script keeps its options, exit statuses and console output (stdout and stderr). The trace writes only to `log.trace`.
- **Invariant 3:** tracing never fails a run. If the trace can't be written or the derivation can't be recorded, the step proceeds exactly as before.
- **Invariant 4:** the unit suite passes unchanged, and `check_api.py`, `check_links.py`, `run_examples.py` and `mkdocs build --strict` pass.
- No edit to `plumetools/cai2012/`, `plumetools/markelov1999/`, any archived case directory (`1d`, `2d-planar`, `2d-wedge`, `caseFoamEx`, `wake-cylinder`), or the legacy ISS / particle-resolution cases. No family package imports `plumetools.trace`.
- `plumetools` runs on Python ≥ 3.9. Every new module starts with `from __future__ import annotations`, with no `match`, no `dataclass(slots=...)`, and no runtime use of `X | Y` types. No new runtime dependency.
- Teaching budget: serial, 1 core, < 2 GB, solver ≤ 3 minutes wall-clock. Cost guard: ≤ 20,000 cells, ≤ 5,000 steps, ≤ 5 × 10⁴ parcels (estimated, steady state).
- `log.trace` line grammar, verbatim from the spec:
  ```
  == <script>  <ISO-8601 time>  [<options>]  host=<hostname>  plumetools=<short sha>|unknown  OpenFOAM=<version>|none
     derive   <entry>  <value> [<unit>]  <formula>  <producer>
     decide   <text>
     step     <tool>  ok|FAILED(<status>)  <elapsed> s  <log file>
     note     <text>
  == end <script>  exit=<status>  <elapsed> s
  ```
- New files use LF line endings. The repository has no `* text=auto` rule, so on Windows checkouts stage with LF.
- Work on branch `feature/traceability-phase1`, created from `docs/traceability-learning-spec`. Run bash-dependent tests under WSL or Linux (they skip on native Windows, as `test_allrun_restart.py` does). Commit messages end with the session's attribution trailer lines.
- Phase 1 traced scripts (the `TRACED_SCRIPTS` list used by the tests):
  - `cases/{cai2012,cai2012-health,cai2012-kn-extension}/baseCase/{Allmesh,Allrun,Allpost}`
  - `cases/r-4d/{baseCase,baseCaseHPC,baseCaseFreePlumeHPC}/{Allmesh,Allrun,Allpost}`
  - the teaching cases' copies

## Decisions made while planning (refinements of the spec)

Each one follows from reading the code. The owner should confirm them during plan review.

1. **The derivation is recorded by a separate call.** `Allmesh` runs `python -m plumetools.trace record . <family>` after its dictionaries are written, instead of handing over in-process objects (spec §6.4). r-4d's `Allmesh` calls `python -m plumetools.nozzleplume setup .` and has no inline Python block. Passing objects in-process would mean editing `nozzleplume`, and the spec forbids family packages importing the trace. `record` re-derives with the map's `objects()`, and the value check proves those objects match what the family writes.
2. **`TraceEntry` gains a `writer` field:** the function that puts the value into the file, alongside `producer`, the function that computes it. The pages show both.
3. **`Allclean` removes `log.trace` only with `--mesh`.** The trace starts with `Allmesh`, whose outputs `--mesh` removes. `Allclean` writes no trace of its own, because an end line would re-create the file it just removed.
4. **Lesson 4's first leg uses `run_time_multiplier: 0.75`, not 0.5.** At 0.5, `derive_run_settings` refuses the case: the end time equals the averaging start of 2 transits. At 0.75 the first leg ends at 1,239 steps with 3 writes every 413 steps, and the resume reaches 1,652 steps and 4 writes. Both were verified with `--dict-only` at `3612457`.
5. **The completeness check covers numeric entries** (a number or a numeric vector, optionally `uniform`). Word entries that carry physics (collision, wall and inflow models) are mapped, but completeness doesn't enforce them.
6. **`0/` values are checked without a mesh.** The value check writes `0/` with the family's own field writer and synthetic `PatchInfo` patches, as `tests/unit/test_cai2012_foamio.py` does.
7. **Teaching cases default to one subdomain** (`n_subdomains: 1`), so a plain `./Allrun` is serial. Lesson 4 uses `./Allrun -np 2`.
8. **`StaticEntry` (spec §6.1) is deferred to the Phase 3 plan.** Nothing in Phase 1 uses it.
9. **The canonical bash block lives in `plumetools/trace/runner_block.sh`.** Every traced script embeds it verbatim, and a test enforces that.
10. **`run_examples.py`'s existing `r-4d Allmesh --dict-only` check gains `cases/r-4d/baseCase/log.trace` in its declared writes.** Without it, that check fails once the trace exists.

## Review Focus

The five conditions most likely to bite a user that no other test exercises. Each test is added to the task named.

1. **A case path containing a space.** The run trace still writes, and console output is unchanged. Test: `test_a_case_path_with_a_space_is_traced` (Task 8).
2. **Read-only invocations must create no file.** That covers `--help`, an unknown option, and r-4d's `Allmesh --describe`. Tests: `test_help_and_bad_options_write_no_trace` (Task 8), `test_describe_writes_no_trace` (Task 9).
3. **A failing OpenFOAM step.** It is recorded as `FAILED(<n>)` with `exit=1`, and the console error is unchanged. Test: `test_a_failing_step_is_recorded_and_the_console_is_unchanged` (Task 8).
4. **A resumed or repeated run appends a block.** It never truncates earlier ones. Test: `test_a_second_run_appends_a_second_block` (Task 8).
5. **A case outside a git checkout, or a machine without `git` or `hostname`.** These give `plumetools=unknown` / `host=unknown`, with nothing on stderr. Test: `test_outside_a_checkout_the_commit_is_unknown_and_silent` (Task 8).

---

### Task 1: Baseline capture of every Phase 1 output (before any change)

**Files:**
- Create: `tests/unit/_runner_stubs.py`
- Create: `tests/invariants/capture_outputs.py`
- Test: `tests/unit/test_capture_outputs.py`

**Interfaces:**
- Produces:
  - `_runner_stubs.STUB_TOOLS: tuple[str, ...]`: `dsmcFoam dsmcInitialise decomposePar reconstructPar mpirun foamDictionary blockMesh snappyHexMesh topoSet createPatch checkMesh postProcess`.
  - `make_stub_bin(root: Path, failing: dict[str, int] | None = None) -> Path`: stubs that append `"<tool> <args>"` to `root/"calls"`. A tool in `failing` exits with the given status. A `python` shim execs `sys.executable`.
  - `script_env(bin_dir: Path, *, path_extra: str = "") -> dict[str, str]`: as in `test_allrun_restart.py`. Stubs come first on `PATH`, `PYTHONPATH` gets the repo, and the `FOAM_*` variables are removed.
  - `make_resume_case(case: Path, script: Path, *, ranks: int = 2, family: str) -> Path`: the `test_allrun_restart.py` fixture generalised to any `Allrun`. For `family="cai2012"` it writes a stub `runInflow.py`; for `"nozzleplume"`, nothing extra.
  - `checkpoint(base: Path, name: str) -> None`: copied from `test_allrun_restart.py`.
  - `run_script(case: Path, script: str, *args: str, env: dict) -> subprocess.CompletedProcess`.
  - `capture_outputs.capture(out_dir: Path) -> list[str]` and `capture_outputs.compare(base: Path, after: Path) -> list[str]`.
  - CLI: `python tests/invariants/capture_outputs.py capture OUT_DIR` and `python tests/invariants/capture_outputs.py compare BASE AFTER`, which exits 1 when there are problems.

- [ ] **Step 1: Write the failing tests** in `tests/unit/test_capture_outputs.py`. Mark the module to skip without bash or on native Windows, as `test_allrun_restart.py` does.
  ```python
  def test_compare_ignores_log_trace_and_nothing_else(tmp_path):
      base, after = tmp_path / "base", tmp_path / "after"
      for root in (base, after):
          (root / "s").mkdir(parents=True); (root / "s" / "controlDict").write_text("deltaT 1;\n")
      (after / "s" / "log.trace").write_text("== Allmesh\n")
      assert capture_outputs.compare(base, after) == []
      (after / "s" / "controlDict").write_text("deltaT 2;\n")
      assert capture_outputs.compare(base, after) == ["s/controlDict: differs"]

  def test_a_scenario_new_since_the_baseline_is_a_note_not_a_problem(tmp_path): ...
      # after/ has an extra top-level scenario directory -> compare returns []

  def test_capture_records_every_phase1_scenario(tmp_path):
      problems = capture_outputs.capture(tmp_path)
      assert problems == []
      for name in capture_outputs.SCENARIOS:
          for f in ("stdout", "stderr", "exit"):
              assert (tmp_path / name / f).is_file()
  ```
- [ ] **Step 2: Run them to see them fail.** `python -m pytest tests/unit/test_capture_outputs.py -q` should fail with `ModuleNotFoundError: capture_outputs`.
- [ ] **Step 3: Implement `_runner_stubs.py` and `capture_outputs.py`.** `capture_outputs` inserts `tests/unit` into `sys.path` to import `_runner_stubs`. Each scenario copies its case into `OUT/<scenario>/case`, runs the command there, and saves `stdout`, `stderr`, `exit` and every file the command wrote. `SCENARIOS`:
  - `Allmesh --dict-only` for each of the six Phase 1 `baseCase`s;
  - r-4d `Allmesh --describe` for each of the three r-4d `baseCase`s;
  - `0/` written by the family field writer with synthetic patches. For Cai that's `foamfields.write_all` with the `PATCHES` shape from `test_cai2012_foamio.py`; for nozzleplume it's `foamfields.write_all` with patches from `cfg.mesh.patch_names`, adding a `wall` plate when `geom.plate` and a `symmetry` patch when there are symmetry planes;
  - with stub tools, for both families' `Allrun`: resume-in-place (`-np 2`), already-at-endTime, and the preflight failure "no constant/polyMesh";
  - `Allpost` with no time directories, which is the die path;
  - `Allclean --dry-run --mesh` after `--dict-only`, deleting any `log.trace` first, so the listing compares like with like. The unit tests of Task 8 cover `log.trace`'s removal.

  `compare` walks `base`, skips any file named `log.trace`, and reports `"<rel>: differs"` or `"<rel>: missing"`. Top-level scenarios only in `after` are printed as notes.
- [ ] **Step 4: Run the tests.** `python -m pytest tests/unit/test_capture_outputs.py -q` should pass.
- [ ] **Step 5: Capture the baseline from the unchanged scripts.** Put it outside the repository: `python tests/invariants/capture_outputs.py capture ~/rarefiedPlume-baseline-phase1`. Expected: exit 0, and one directory per scenario.
- [ ] **Step 6: Commit.** `git add tests/unit/_runner_stubs.py tests/invariants/capture_outputs.py tests/unit/test_capture_outputs.py`, then `git commit -m "test: capture Phase 1 outputs for the before/after comparison"`.

### Task 2: Trace model and dictionary read-back

**Files:**
- Create: `plumetools/trace/__init__.py`, `plumetools/trace/model.py`, `plumetools/trace/readback.py`
- Modify: `docs/_tools/hooks.py` (add `"trace": "trace"` to `FAMILIES`)
- Create: `docs/reference/api/trace/index.md` (the module table; later tasks add a row each)
- Test: `tests/unit/test_trace_model.py`, `tests/unit/test_trace_readback.py`

**Interfaces:**
- Produces, in `plumetools.trace.model` (all `@dataclass(frozen=True)`):
  - `Stage(name: str, module: str, does: str)`
  - `TraceEntry(file: str, key: str, quantity: str, unit: str, formula: str, producer: str, writer: str, value: str, inputs: tuple, tag: str, component: int | None = None, optional: bool = False, explained_in: str | None = None)`. `file` and `key` may contain `{species}`, `{nozzle}` and `{plate}`. `key` is a slash path (`moleculeProperties/{species}/mass`, `boundaryField/{nozzle}/value`). `unit` is `"-"` when dimensionless.
  - `DictionaryEntry(file: str, contents: str, plan: str, writer: str)`
  - `FamilyMap(family: str, title: str, cases: tuple, stages: tuple, entries: tuple, dictionaries: tuple, boilerplate: frozenset, objects: Callable, context: Callable, write_outputs: Callable)`, with methods `entry(file: str, key: str) -> TraceEntry` (raises `KeyError`) and `traced_files(ctx: dict) -> list[str]`, which returns the formatted `file` of every entry, in first-seen order.
  - `TAGS = frozenset({"PAPER", "ASSUMPTION", "DERIVED", "NUMERICAL", "REPO", "BIRD", "EXTENSION", "HEALTH", "LEGACY", "GIVEN", "TEACHING"})`
  - `resolve_value(objects: Mapping[str, Any], path: str) -> Any`. The first segment indexes `objects`, the rest are attributes. It raises `LookupError` naming the full path.
- Produces, in `plumetools.trace.readback`:
  - `read_entry(path: Path, key: str) -> str`: the raw value text without its `;`. Raises `KeyError("<key> not in <path>")`.
  - `numeric_entries(path: Path) -> dict[str, str]`: every entry, at any depth outside `FoamFile`, whose value `parse_numeric` accepts.
  - `parse_numeric(text: str) -> float | tuple | None`: accepts `uniform` + a number, a number, or `( numbers )`. Returns `None` for anything else, including `[dims]`, quoted strings and words.
  - `tolerance(text: str) -> float`: half a unit in the last printed significant digit of the largest-magnitude number in `text`.

- [ ] **Step 1: Write the failing tests.** Use fixture files written from these literal snippets, taken from a generated Cai case: a `dsmcProperties` with the `FoamFile` header, `//` comments and nested `moleculeProperties`; a `controlDict` with `functions { fieldAverage1 { timeStart 0.0501435; fields ( rhoN { mean on; } ); } }`; a `decomposeParDict` with `coeffs { n (1 2 2); }`; a `0/boundaryU` with `dimensions [0 1 -1 0 0 0 0];`, `internalField uniform (0 0 0);` and `boundaryField { nozzle { type fixedValue; value uniform (707.0 0 0); } vacuum { type zeroGradient; } }`.
  ```python
  assert read_entry(props, "nEquivalentParticles") == "3.310039e+09"
  assert read_entry(props, "moleculeProperties/Ar/mass") == "6.63e-26"
  assert read_entry(control, "functions/fieldAverage1/timeStart") == "0.0501435"
  assert read_entry(u, "boundaryField/nozzle/value") == "uniform (707.0 0 0)"
  with pytest.raises(KeyError, match="nope not in"): read_entry(props, "nope")
  ents = numeric_entries(decompose); assert ents["coeffs/n"] == "(1 2 2)" and "method" not in ents
  assert "FoamFile/version" not in numeric_entries(props)
  assert "dimensions" not in numeric_entries(u) and numeric_entries(u)["internalField"] == "uniform (0 0 0)"
  assert parse_numeric("uniform (707 0 0)") == (707.0, 0.0, 0.0)
  assert parse_numeric("1e+10") == 1e10 and parse_numeric("VariableHardSphere") is None
  assert parse_numeric('( "libplumeDsmcBoundaryModels.so" )') is None
  assert tolerance("3.310039e+09") == pytest.approx(500.0)
  assert tolerance("0.0586413") == pytest.approx(5e-8)
  assert tolerance("4") == 0.5
  # model
  assert resolve_value({"run": SimpleNamespace(dt=2.0)}, "run.dt") == 2.0
  with pytest.raises(LookupError, match="run.nope"): resolve_value({"run": SimpleNamespace()}, "run.nope")
  ```
  Also check that a `FamilyMap` built in the test returns the right entry from `entry()`, raises `KeyError` for an unknown one, and that `traced_files({"species": "Ar"})` formats placeholders.
- [ ] **Step 2: Run them to see them fail** (`python -m pytest tests/unit/test_trace_readback.py tests/unit/test_trace_model.py -q`). Expect `ModuleNotFoundError: plumetools.trace`.
- [ ] **Step 3: Implement.** `readback` tokenises the file after stripping `//` and `/* */` comments, tracks `{}` nesting and `()` lists, and keys entries by their brace path. Entries inside a `( … )` list that hold `{}` blocks (the `fields ( rhoN { … } )` form) are keyed `<list key>/<name>/<entry>`. `__init__.py` holds only a docstring pointing to `docs/trace/index.md`.
- [ ] **Step 4: Run the tests** (expect PASS), then `python docs/_tools/check_api.py` (expect exit 0, with the trace index listing `plumetools.trace`, `plumetools.trace.model` and `plumetools.trace.readback` and their principal interfaces).
- [ ] **Step 5: Commit** the five files and the hooks change: `git commit -m "feat(trace): trace model and OpenFOAM dictionary read-back"`.

### Task 3: The checks

**Files:**
- Create: `plumetools/trace/checks.py`
- Modify: `docs/reference/api/trace/index.md` (add the `checks` row)
- Test: `tests/unit/test_trace_checks.py`

**Interfaces:**
- Consumes: `FamilyMap`, `TraceEntry`, `resolve_value` (Task 2); `read_entry`, `numeric_entries`, `parse_numeric`, `tolerance` (Task 2).
- Produces, in `plumetools.trace.checks`. Every function returns a list of problem strings, each formatted `"<family> <case>: <file> : <key>: <message>"`, and an empty list means pass:
  - `accepted_keys(cfgs: Iterable) -> set[str]`: `"section.key"` for every field of every dataclass section of the loaded configs, plus the top-level non-dataclass fields. It takes the union over all instances, so a section that's `None` in one case (no plate) is still known from another.
  - `check_names(m: FamilyMap) -> list[str]`: every `producer` and `writer` imports and resolves. Message: `does not resolve`.
  - `check_inputs(m: FamilyMap, repo: Path) -> list[str]`: every input is in `accepted_keys` of all `m.cases`. Message: `not a key the loader accepts`.
  - `check_values(m: FamilyMap, repo: Path, work: Path) -> list[str]`: for each case, `objs = m.objects(repo/case)`, `ctx = m.context(objs)`, `m.write_outputs(objs, work/<case>)`, then every entry is read back and compared.
    - Numbers must agree within `tolerance` (per component for vectors; `component` selects one).
    - Words must be equal.
    - Message: `written <text> but the map says <value>`.
    - An `optional` entry is skipped when its value resolves to `None`, fails to resolve, or its key is absent. A non-optional absent key gives `missing`.
  - `check_completeness(m: FamilyMap, repo: Path, work: Path) -> list[str]`: in each traced file, every `numeric_entries` key that is neither a formatted entry key nor in the formatted `boilerplate` (`"<file>:<key>"`). Message: `not traced and not boilerplate`.
  - `check_tags(m: FamilyMap) -> list[str]`: every tag must be in `TAGS`.
  - `run_checks(m: FamilyMap, repo: Path, work: Path) -> list[str]`: all of the above, concatenated.

- [ ] **Step 1: Write the failing tests** against a toy family built in the test. Its `objects` returns `{"cfg": <dataclass with section dsmc: (n: int, dt: float)>, "run": SimpleNamespace(dt=0.5, v=(3.0, 0.0, 0.0))}`, and `write_outputs` writes a `system/controlDict` with `deltaT 0.5; startTime 0; extra 7;` and a `0/boundaryU` with `boundaryField { nozzle { value uniform (3 0 0); } }`.
  ```python
  def test_an_unresolvable_producer_is_reported(...): assert "does not resolve" in one(check_names(bad))
  def test_an_input_the_loader_does_not_accept_is_reported(...): assert "not a key the loader accepts" in one(...)
  def test_a_wrong_value_is_reported_with_both_numbers(...): assert "written 0.5 but the map says 0.25" in one(...)
  def test_a_value_within_the_printed_precision_passes(...)        # map 0.50000004 vs written "0.5"
  def test_a_vector_component_is_compared(...)                      # component=0 against (3 0 0)
  def test_a_missing_optional_entry_is_skipped_and_a_required_one_is_not(...)
  def test_an_untraced_number_is_reported_unless_boilerplate(...):  # "extra" reported; "startTime" allowlisted
  def test_an_unknown_tag_is_reported(...)
  def test_accepted_keys_union_covers_a_section_absent_in_one_case(...)
  ```
- [ ] **Step 2: Run them to see them fail.** Expect `ModuleNotFoundError: plumetools.trace.checks`.
- [ ] **Step 3: Implement** `checks.py` with the signatures above.
- [ ] **Step 4: Run** `python -m pytest tests/unit/test_trace_checks.py -q` (expect PASS) and `python docs/_tools/check_api.py` (expect exit 0).
- [ ] **Step 5: Commit:** `git commit -m "feat(trace): producer, input, value, completeness and tag checks"`.

### Task 4: The `cai2012` trace map

**Files:**
- Create: `plumetools/trace/maps/__init__.py`, `plumetools/trace/maps/cai2012.py`
- Modify: `docs/reference/api/trace/index.md` (add the `maps` and `maps.cai2012` rows)
- Test: `tests/unit/test_trace_maps.py`

**Interfaces:**
- Consumes: the model and checks (Tasks 2–3), plus the public names of `plumetools.cai2012`: `config.load_case_config`, `geometry.from_config`, `inflow.from_config`, `inflow.derive_run_settings`, `mesh.plan`, `mesh.predicted_nozzle_faces`, `mesh.write_mesh_setup`, `dictionaries.write_case_dictionaries`, `dictionaries.simple_split`, `foamfields.write_all`, and `plumetools.mesh.boundary.PatchInfo`.
- Produces:
  - `plumetools.trace.maps.MAPS: dict[str, str]`: family to module, `{"cai2012": "plumetools.trace.maps.cai2012"}` at this task.
  - `load_map(family: str) -> FamilyMap`, raising `KeyError` listing the known families.
  - `families() -> list[str]`.
  - `plumetools.trace.maps.cai2012.MAP: FamilyMap`, with `cases = ("cases/cai2012/baseCase", "cases/cai2012-health/baseCase", "cases/cai2012-kn-extension/baseCase")`. All three run the same `Allmesh`, so all three get checked. Task 10 appends the teaching case.

`MAP.objects(case_dir)` repeats `cases/cai2012/baseCase/Allmesh`'s Python block, in its order:
- `cfg = load_case_config(case_dir)`
- `geom = geometry.from_config(cfg)`
- `exit_state = inflow.from_config(cfg)`
- `plan = mesh.plan(cfg, geom, exit_state)`
- `run = inflow.derive_run_settings(cfg, exit_state, geom, min_cell_size_m=plan.min_cell_size_m, exit_cell_volume_m3=plan.exit_cell_volume_m3)`
- `nozzle = mesh.predicted_nozzle_faces(plan, geom)`
- `split = dictionaries.simple_split(cfg.dsmc.n_subdomains)`

It returns these as a dict.

`MAP.context(objs)` returns `{"species": cfg.gas.species_name, "nozzle": cfg.mesh.patch_names["nozzle"]}`.

`MAP.write_outputs(objs, out)` calls `mesh.write_mesh_setup`, `dictionaries.write_case_dictionaries` and `foamfields.write_all`. The patches are `{nozzle: PatchInfo(nozzle, "patch", nozzle.n_faces, 0)}`, plus `upstream_vacuum` and `outer` as `"patch"`, plus `symmetry` as `"symmetry"` when `cfg.geometry.symmetry_mode == "half_y"`.

Entries. For every entry, `inputs` are the `case.yaml` keys the producer reads, so read the producer to list them. The tag is the one `case.yaml`'s comment gives that input, and `DERIVED` for computed values. Formulas are copied as written here.

| file | key | value | producer | writer | unit | formula |
|---|---|---|---|---|---|---|
| `constant/dsmcProperties` | `nEquivalentParticles` | `run.n_equivalent_particles` | `plumetools.cai2012.inflow.derive_run_settings` | `plumetools.cai2012.dictionaries.render_dsmc_properties` | - | `W = n0 · V_exit / N_exit ÷ numerical_particle_multiplier; dsmc.n_equivalent_particles overrides` |
| `constant/dsmcProperties` | `BinaryCollisionModel` | `cfg.dsmc.binary_collision_model` | `plumetools.cai2012.config.load_case_config` | same writer | - | `case.yaml value` |
| `constant/dsmcProperties` | `InflowBoundaryModel` | `cfg.output.inflow_model` | loader | same | - | `case.yaml value` |
| `constant/dsmcProperties` | `VariableHardSphereCoeffs/Tref` | `cfg.gas.t_ref_K` | loader | same | K | `case.yaml value` |
| `constant/dsmcProperties` | `moleculeProperties/{species}/mass` | `cfg.gas.mass_kg` | loader | same | kg | `case.yaml value` |
| `constant/dsmcProperties` | `moleculeProperties/{species}/diameter` | `cfg.gas.diameter_m` | loader | same | m | `case.yaml value` |
| `constant/dsmcProperties` | `moleculeProperties/{species}/omega` | `cfg.gas.omega` | loader | same | - | `case.yaml value` |
| `constant/dsmcProperties` | `moleculeProperties/{species}/internalDegreesOfFreedom` | `cfg.gas.internal_degrees_of_freedom` | loader | same | - | `case.yaml value` |
| `system/controlDict` | `deltaT` | `run.delta_t_s` | `…inflow.derive_run_settings` | `plumetools.cai2012.dictionaries.render_control_dict` | s | `Δt = C · Δx_min / (U0 + 3σ); dsmc.delta_t_s overrides` |
| `system/controlDict` | `endTime` | `run.end_time_s` | same | same | s | `t_end = run_time_multiplier · (t_avg + N_sample · L_x / U0), whole steps and write intervals` |
| `system/controlDict` | `writeInterval` | `run.write_interval_steps` | same | same | steps | `steps per write, so the last write lands on t_end` |
| `system/controlDict` | `functions/fieldAverage1/timeStart` | `run.average_start_s` | same | same | s | `t_avg = 10⁴ t0 (basis cai) or N_transient · L_x / U0 (basis transits); dsmc.average_start_s overrides` |
| `system/dsmcInitialiseDict` | `numberDensities/{species}` | `cfg.dsmc.initial_number_density_per_m3` | loader | `…dictionaries.render_dsmc_initialise_dict` | 1/m^3 | `case.yaml value` |
| `system/dsmcInitialiseDict` | `temperature` | `cfg.dsmc.initial_temperature_K` | loader | same | K | `case.yaml value` |
| `system/decomposeParDict` | `numberOfSubdomains` | `cfg.dsmc.n_subdomains` | loader | `…dictionaries.render_decompose_par_dict` | - | `case.yaml value` |
| `system/decomposeParDict` | `coeffs/n` | `split` | `plumetools.cai2012.dictionaries.simple_split` | same | - | `(nx ny nz), nx·ny·nz = numberOfSubdomains, spread over the graded shell` |
| `0/boundaryU` | `boundaryField/{nozzle}/value` (component 0) | `exit_state.velocity_m_per_s` | `plumetools.cai2012.gas.speed_from_speed_ratio` | `plumetools.cai2012.foamfields.write_inflow_fields` | m/s | `U0 = S0 · √(2 R T0)` |
| `0/boundaryT` | `boundaryField/{nozzle}/value` | `exit_state.T0_K` | loader | same | K | `case.yaml value` |
| `0/boundaryNumberDensity_{species}` | `boundaryField/{nozzle}/value` | `exit_state.number_density_per_m3` | `plumetools.cai2012.gas.number_density_from_kn` | same | 1/m^3 | `n0 from λ0 = Kn · L (L = D), VHS convention` |

`dictionaries`, the dictionary-level entries:
- `system/blockMeshDict`: plan `plan`, writer `plumetools.cai2012.mesh.render_block_mesh_dict`
- `system/topoSetDict`: writer `…mesh.render_topo_set_dict`
- `system/createPatchDict`: writer `…mesh.render_create_patch_dict`
- `system/meshQualityDict`: writer `…mesh.render_mesh_quality_dict`
- `system/fvSchemes`, `system/fvSolution`: writers `…dictionaries.render_fv_schemes` / `render_fv_solution`
- the nine `0/` measurement fields: writer `…foamfields.write_measurement_fields`

`boilerplate`:
- `system/controlDict:startTime`, `purgeWrite`, `writePrecision`, `timePrecision`
- `system/dsmcInitialiseDict:velocity`
- `internalField` of the three inflow fields

Add any other entry only if the completeness check names it and it is a literal constant in the writer.

`stages` uses these eight rows, verbatim:

| name | module | does |
|---|---|---|
| load | `plumetools.cai2012.config.load_case_config` | validate `case.yaml`; unknown keys refused |
| derive | `plumetools.cai2012.inflow`, `.geometry`, `.mesh.plan` | exit state, geometry, mesh plan, run settings |
| assess | `plumetools.cai2012.checks.run_checks` | DSMC quality checks; a failure stops `Allmesh` |
| write dictionaries | `plumetools.cai2012.mesh.write_mesh_setup`, `.dictionaries.write_case_dictionaries` | `system/`, `constant/dsmcProperties` |
| mesh | `Allmesh`: `blockMesh`, `topoSet`, `createPatch`, `checkMesh` | `constant/polyMesh`, nozzle carved and verified |
| initialise | `runInflow.py` → `plumetools.cai2012.foamfields.write_all`; `dsmcInitialise` | `0/` |
| solve | `Allrun`: `dsmcFoam` (with `decomposePar`, `mpirun`, `reconstructPar` in parallel) | time directories |
| post | `Allpost`: `postProcess`, `postProcess.py` → `plumetools.cai2012.post` | `results/` |

- [ ] **Step 1: Write the failing test.**
  ```python
  @pytest.mark.parametrize("family", families())
  def test_every_map_passes_every_check(family, tmp_path):
      assert run_checks(load_map(family), REPO, tmp_path) == []
  def test_load_map_names_the_known_families_for_an_unknown_one():
      with pytest.raises(KeyError, match="cai2012"): load_map("nope")
  ```
- [ ] **Step 2: Run it to see it fail** (`python -m pytest tests/unit/test_trace_maps.py -q`). Expect an import error.
- [ ] **Step 3: Implement** `maps/__init__.py` and `maps/cai2012.py` with the tables above.
- [ ] **Step 4: Run** the test and `check_api.py`. Expect PASS and exit 0. A completeness problem means either a missing entry or an unlisted literal constant. Resolve it by the rule above and record which in the commit message.
- [ ] **Step 5: Commit:** `git commit -m "feat(trace): the cai2012 trace map"`.

### Task 5: The `nozzleplume` trace map

**Files:**
- Create: `plumetools/trace/maps/nozzleplume.py`
- Modify: `plumetools/trace/maps/__init__.py` (register), `docs/reference/api/trace/index.md`
- Test: `tests/unit/test_trace_maps.py` (the parametrised test now covers it)

**Interfaces:**
- Consumes: `plumetools.nozzleplume.workflow.load` (it returns `Case` with `cfg geom exit_state plan nozzle run`), plus `mesh.write_mesh_setup(case_dir, cfg, geom, plan)`, `dictionaries.write_case_dictionaries(case_dir, cfg, geom, exit_state, run)`, `dictionaries.simple_split`, `foamfields.write_all(case_dir, cfg, geom, exit_state, patches)` and `PatchInfo`.
- Produces:
  - `plumetools.trace.maps.nozzleplume.MAP: FamilyMap`, with `cases = ("cases/r-4d/baseCase", "cases/r-4d/baseCaseHPC", "cases/r-4d/baseCaseFreePlumeHPC")`. Task 11 appends the teaching case.
  - `objects`: the `Case` fields as a dict, plus `split = simple_split(n)` when `cfg.dsmc.decomposition_method == "simple"`, else `None`.
  - `context`: `species` = `cfg.gas.species_name`, `nozzle` and `plate` from `cfg.mesh.patch_names`.
  - `write_outputs`: patches from `cfg.mesh.patch_names` (nozzle, upstream_vacuum, outer as `patch`), plus `plate` as `wall` when `geom.plate`, plus `symmetry` as `symmetry` when `geom.symmetry_planes`.

Entries. Same columns and rules as Task 4. Producers are `plumetools.nozzleplume.run.derive` for run values, `plumetools.nozzleplume.exitstate.from_config` for exit values, and `plumetools.nozzleplume.config.load_case_config` for straight `case.yaml` values. Writers are `plumetools.nozzleplume.dictionaries.render_dsmc_properties`, `render_control_dict`, `render_dsmc_initialise_dict` and `render_decompose_par_dict`, and `plumetools.nozzleplume.foamfields.write_inflow_fields`.

| file | key | value | unit | formula | optional |
|---|---|---|---|---|---|
| `constant/dsmcProperties` | `nEquivalentParticles` | `run.n_equivalent_particles` | - | `dsmc.n_equivalent_particles, or n · V_exit / target_particles_per_exit_cell` | |
| `constant/dsmcProperties` | `BinaryCollisionModel` | `cfg.dsmc.collision_model` | - | `case.yaml value` | |
| `constant/dsmcProperties` | `WallInteractionModel` | `cfg.dsmc.wall_interaction_model` | - | `case.yaml value` | |
| `constant/dsmcProperties` | `VariableHardSphereCoeffs/Tref` | `cfg.gas.t_ref_K` | K | `case.yaml value` | yes |
| `constant/dsmcProperties` | `moleculeProperties/{species}/mass` | `cfg.gas.molecular_mass_kg` | kg | `case.yaml value` | |
| `constant/dsmcProperties` | `moleculeProperties/{species}/diameter` | `cfg.gas.molecular_diameter_m` | m | `case.yaml value` | |
| `constant/dsmcProperties` | `moleculeProperties/{species}/omega` | `cfg.gas.omega` | - | `case.yaml value` | |
| `constant/dsmcProperties` | `moleculeProperties/{species}/internalDegreesOfFreedom` | `cfg.gas.internal_degrees_of_freedom` | - | `case.yaml value` | |
| `system/controlDict` | `deltaT` | `run.delta_t_s` | s | `Δt = C · Δx_min / (U + 3σ); dsmc.delta_t_s overrides` | |
| `system/controlDict` | `endTime` | `run.end_time_s` | s | `t_end = t_avg + N_sample · x_max / U, whole steps and write intervals` | |
| `system/controlDict` | `writeInterval` | `run.write_interval_steps` | steps | `steps per write, so the last write lands on t_end` | |
| `system/controlDict` | `functions/fieldAverage1/timeStart` | `run.average_start_s` | s | `t_avg = N_transient · x_max / U; dsmc.average_start_s overrides` | |
| `system/dsmcInitialiseDict` | `numberDensities/{species}` | `cfg.dsmc.initial_number_density_per_m3` | 1/m^3 | `case.yaml value` | |
| `system/dsmcInitialiseDict` | `temperature` | `cfg.dsmc.initial_temperature_K` | K | `case.yaml value` | |
| `system/decomposeParDict` | `numberOfSubdomains` | `cfg.dsmc.n_subdomains` | - | `case.yaml value` | |
| `system/decomposeParDict` | `coeffs/n` | `split` | - | `(nx ny nz), nx·ny·nz = numberOfSubdomains` | yes |
| `0/boundaryU` | `boundaryField/{nozzle}/value` (component 0) | `exit_state.speed_m_per_s` | m/s | `case.yaml exit speed (or derived from the given form)` | |
| `0/boundaryT` | `boundaryField/{nozzle}/value` | `exit_state.temperature_K` | K | `case.yaml value` | |
| `0/boundaryT` | `boundaryField/{plate}/value` | `geom.plate.temperature_K` | K | `case.yaml plate.temperature_K (diffuse wall)` | yes |
| `0/boundaryNumberDensity_{species}` | `boundaryField/{nozzle}/value` | `exit_state.number_density_per_m3` | 1/m^3 | `case.yaml exit density (or derived from the given form)` | |

The read-back check confirms the `functions/fieldAverage1` and `VariableHardSphereCoeffs` paths. If the writer names a block differently, use the writer's name. Dictionary-level entries cover `blockMeshDict`, `snappyHexMeshDict`, `topoSetDict`, `createPatchDict` and `meshQualityDict`, with the matching `plumetools.nozzleplume.mesh.render_*` writer, plus `fvSchemes`, `fvSolution` and the `0/` measurement fields. Start `boilerplate` from Task 4's list; the plate's `boundaryU` value is a likely addition. The `stages` rows are Task 4's, with the family's own names: `workflow.load`, `checks.run_checks`, `setup`, `Allmesh` with `snappyHexMesh`, `Allrun`'s `python -m plumetools.nozzleplume fields`, and `Allpost`'s `post`/`diagnostics`.

- [ ] **Step 1:** Register `"nozzleplume"` in `MAPS`, then run `python -m pytest tests/unit/test_trace_maps.py -q`. Expect FAIL with an import error for the new module.
- [ ] **Step 2: Implement** `maps/nozzleplume.py`.
- [ ] **Step 3: Run** the test (expect PASS for both families, all three r-4d cases) and `check_api.py`.
- [ ] **Step 4: Commit:** `git commit -m "feat(trace): the nozzleplume trace map"`.

### Task 6: Trace pages, the drift check and the `pages`/`check` commands

**Files:**
- Create: `plumetools/trace/pages.py`, `plumetools/trace/__main__.py`, `docs/trace/index.md`
- Generate and commit: `docs/trace/cai2012.md`, `docs/trace/nozzleplume.md`
- Modify: `mkdocs.yml` (a top-level `Trace` section: `trace/index.md`, then `cai2012: trace/cai2012.md`, `nozzleplume (r-4d): trace/nozzleplume.md`), `.github/workflows/docs.yml` (in the `examples` job, after the install step, add the step `python -m plumetools.trace pages --check`), `docs/reference/api/trace/index.md`
- Test: `tests/unit/test_trace_pages.py`

**Interfaces:**
- Consumes: `load_map`, `families` (Task 4); `FamilyMap`, `resolve_value` (Task 2).
- Produces:
  - `page_path(repo: Path, family: str) -> Path`: `repo/"docs"/"trace"/f"{family}.md"`.
  - `render_page(m: FamilyMap, repo: Path) -> str`.
  - `write_pages(repo: Path) -> list[Path]`, written with `newline="\n"`.
  - `stale_pages(repo: Path, pages_dir: Path | None = None) -> list[Path]`.
  - CLI `python -m plumetools.trace pages [--check]` and `python -m plumetools.trace check [FAMILY ...]`. `check` runs `run_checks` for each family in a temporary directory, prints problems, and exits 1 on any. `record` is added in Task 7.

**Page layout,** in this order:
- `# Trace: <title>`
- the line `<!-- generated by python -m plumetools.trace pages from plumetools/trace/maps/<family>.py; do not edit -->`
- a paragraph naming the example case, `m.cases[0]`, and linking `index.md`
- `## How a case runs in this family`: the stages table, with columns *Stage · Module · Does*
- one `## <file>` section per traced file, in `traced_files` order, each a table with columns *Entry · Quantity · Formula · Computed by · Written by · From `case.yaml` · Tag · Value in `<case>`*. Values use `:.6g` plus a non-`-` unit; strings in code font; a vector shows the selected component, or the whole tuple. An optional entry that doesn't resolve for the example case shows `—`.
- `## Mesh dictionaries and other generated files`: columns *File · Contents · Planned by · Written by*
- `## Not traced: fixed by the writer`: the boilerplate keys as a list

**Drift message,** exact and printed to stderr by `pages --check`, one per stale page, exiting 1:
```
docs/trace/<family>.md is out of date with plumetools/trace/maps/<family>.py.
Regenerate it:  python -m plumetools.trace pages
```

`docs/trace/index.md` is hand-written:
- *What a trace page is* (the columns, explained in a sentence each);
- *What `log.trace` is* (the line grammar above, verbatim, plus an annotated example block);
- *Pages* (links to the two family pages; Phase 2 and 3 families are listed as "not yet traced").

- [ ] **Step 1: Write the failing tests.**
  ```python
  def test_committed_pages_match_the_maps():
      assert stale_pages(REPO) == []
  def test_a_stale_page_is_reported(tmp_path):
      shutil.copytree(REPO / "docs" / "trace", tmp_path / "trace")
      (tmp_path / "trace" / "cai2012.md").write_text("old\n")
      assert [p.name for p in stale_pages(REPO, tmp_path / "trace")] == ["cai2012.md"]
  def test_a_page_names_each_entry_its_producer_and_its_writer():
      text = render_page(load_map("cai2012"), REPO)
      for s in ("nEquivalentParticles", "`plumetools.cai2012.inflow.derive_run_settings`",
                "`plumetools.cai2012.dictionaries.render_dsmc_properties`", "## How a case runs in this family"):
          assert s in text
  def test_pages_are_written_with_lf(tmp_path): ...
  def test_pages_check_cli_prints_the_regenerate_command(tmp_path, capsys): ...   # main(["pages", "--check"]) on a stale copy -> 1
  ```
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement** `pages.py` and `__main__.py`, then write the pages with `python -m plumetools.trace pages`.
- [ ] **Step 4: Verify.**
  - `python -m pytest tests/unit/test_trace_pages.py -q` should PASS.
  - `python -m plumetools.trace pages --check` should exit 0.
  - `python -m plumetools.trace check` should exit 0.
  - `python docs/_tools/check_api.py` should pass.
  - `python docs/_tools/check_links.py --markdown-only` should pass.
  - Read both generated pages once in full, and confirm every formula against its producer's code. The value check can't verify formula text.
- [ ] **Step 5: Commit** the code, pages, index, nav and workflow step: `git commit -m "feat(trace): generated trace pages with a drift check"`.

### Task 7: `python -m plumetools.trace record`, the `derive` lines

**Files:**
- Create: `plumetools/trace/record.py`
- Modify: `plumetools/trace/__main__.py` (add the `record` subcommand), `docs/reference/api/trace/index.md`
- Test: `tests/unit/test_trace_record.py`

**Interfaces:**
- Consumes: `load_map`, `resolve_value`.
- Produces:
  - `derive_lines(m: FamilyMap, case_dir: Path) -> list[str]`: one line per `TraceEntry` whose value resolves to something other than `None`, in map order. The format is `f"   derive   {entry:<44} {value:<16} {formula}   {producer}"`:
    - `entry` is `"<file basename>:<formatted key>"`, for example `dsmcProperties:nEquivalentParticles` or `boundaryT:boundaryField/nozzle/value`;
    - `value` is a number as `f"{v:.3g}"`, followed by a space and the unit unless the unit is `"-"`, or a word verbatim, or the selected vector component;
    - `producer` is the dotted name without its `plumetools.` prefix.
  - `record(case_dir: Path, family: str) -> None`: appends the lines to `case_dir/"log.trace"`, opened with `newline="\n"`.
  - CLI `python -m plumetools.trace record CASE FAMILY`: exits 0 on success. On any exception it writes one line to stderr, exits 1, and leaves the file untouched.

- [ ] **Step 1: Write the failing tests.**
  ```python
  DERIVE = re.compile(r"^   derive   \S+\s+\S.*\S\s{3}\S+$")
  def test_record_appends_one_derive_line_per_resolvable_entry(tmp_path):
      case = copy_case("cases/cai2012/baseCase", tmp_path)
      record(case, "cai2012"); lines = (case / "log.trace").read_text().splitlines()
      assert len(lines) == len(load_map("cai2012").entries) and all(DERIVE.match(l) for l in lines)
      assert any(l.split()[1] == "dsmcProperties:nEquivalentParticles" for l in lines)
  def test_a_second_record_appends_rather_than_replaces(tmp_path): ...   # 2x the lines
  def test_an_optional_entry_that_does_not_resolve_is_left_out(tmp_path):  # baseCaseFreePlumeHPC: no boundaryT:…/plate/value line
  def test_the_cli_fails_quietly_and_writes_nothing_for_a_missing_case(tmp_path, capsys):
      assert main(["record", str(tmp_path / "nope"), "cai2012"]) == 1
      assert not (tmp_path / "nope" / "log.trace").exists() and len(capsys.readouterr().err.splitlines()) == 1
  ```
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement** `record.py` and the subcommand.
- [ ] **Step 4: Run** the tests and `check_api.py`. Expect PASS.
- [ ] **Step 5: Commit:** `git commit -m "feat(trace): record the derivation into log.trace"`.

### Task 8: The run-trace block, and the Cai runner scripts

**Files:**
- Create: `plumetools/trace/runner_block.sh`
- Modify: `pyproject.toml` (`[tool.setuptools.package-data]`: `plumetools = ["viz/*.yaml", "trace/*.sh"]`)
- Modify: `cases/cai2012/baseCase/{Allmesh,Allrun,Allpost,Allclean}`, then copy each byte-for-byte over `cases/cai2012-health/baseCase/` and `cases/cai2012-kn-extension/baseCase/`
- Test: `tests/unit/test_run_trace.py`

**Interfaces:**
- Consumes: `_runner_stubs` (Task 1); the `record` CLI (Task 7).
- Produces:
  - `TRACED_SCRIPTS` in `tests/unit/test_run_trace.py`: the paths listed in Global Constraints, extended by Tasks 9–11.
  - `block_of(script: Path) -> str` in the same module: the text from the `# >>> run trace` line through the `# <<< run trace` line, inclusive, with a trailing newline.
  - The bash functions `trace`, `trace_clock`, `trace_since`, `trace_begin`, `trace_decide`, `trace_step`, `trace_cmd` and `trace_derivation`, available in every traced script.

`plumetools/trace/runner_block.sh`, exact text. Each traced script embeds these lines verbatim, from the `# >>>` line through the `# <<<` line:

```bash
# >>> run trace ------------------------------------------------------------- #
# Appends one block per run of this script to log.trace (docs/trace/index.md).
# It never prints and never changes an exit status: every write is guarded and
# its errors are discarded. Canonical copy: plumetools/trace/runner_block.sh;
# a unit test keeps every traced script's copy identical to it.
TRACE_FILE=log.trace
TRACE_ARGS="$*"
TRACE_T0=""
trace() { { printf '%s\n' "$*" >> "$TRACE_FILE"; } 2>/dev/null || true; }
trace_clock() {
    local t
    t=$(date +%s.%N 2>/dev/null)
    case "$t" in ''|*N*) date +%s ;; *) echo "$t" ;; esac
}
trace_since() {
    awk -v a="$1" -v b="$(trace_clock)" 'BEGIN { printf "%.1f", b - a }' 2>/dev/null || echo "?"
}
trace_begin() {
    local sha host
    TRACE_T0=$(trace_clock)
    sha=$(git rev-parse --short HEAD 2>/dev/null) || sha=""
    host=$(hostname 2>/dev/null || uname -n 2>/dev/null) || host=""
    trace "== ${SELF}  $(date +%Y-%m-%dT%H:%M:%S)  ${TRACE_ARGS:+$TRACE_ARGS  }host=${host:-unknown}  plumetools=${sha:-unknown}  OpenFOAM=${WM_PROJECT_VERSION:-none}"
    trap 'TRACE_STATUS=$?; trace "== end ${SELF}  exit=${TRACE_STATUS}  $(trace_since "$TRACE_T0") s"' EXIT
}
trace_decide() { trace "   decide   $*"; }
trace_step() {   # NAME STATUS SECONDS LOG
    local state=ok
    [ "$2" -eq 0 ] || state="FAILED($2)"
    trace "$(printf '   step     %-22s %-10s %7s s   %s' "$1" "$state" "$3" "$4")"
}
trace_cmd() {    # NAME COMMAND...: run it untouched, then record a step line
    local name=$1 t0 status
    shift
    t0=$(trace_clock)
    "$@"
    status=$?
    trace_step "$name" "$status" "$(trace_since "$t0")" "-"
    return "$status"
}
trace_derivation() {   # FAMILY: the derive lines, from python -m plumetools.trace record
    "$PYTHON" -m plumetools.trace record . "$1" >/dev/null 2>&1 \
        || trace "   note     derivation not recorded: python -m plumetools.trace record . $1 failed"
}
# <<< run trace ------------------------------------------------------------- #
```

**`run()` in `Allmesh` and `Allrun`.** Replace only the line `if ! "$@" > "$log" 2>&1; then` with the lines below. The body of that `if` stays byte-for-byte, including `die "'$1' failed; full output in $PWD/$log"`.
```bash
    local trace_t0 status
    trace_t0=$(trace_clock)
    "$@" > "$log" 2>&1
    status=$?
    trace_step "${log#log.}" "$status" "$(trace_since "$trace_t0")" "$log"
    if [ "$status" -ne 0 ]; then
```

**Insertion points** (anchors are existing lines):

| Script | After / at | Insert |
|---|---|---|
| all three | the `cd "$(dirname "$0")" \|\| exit 1 …` line | the block |
| `Allmesh` | `have() { … }` | `trace_begin` |
| `Allmesh` | in the results-exist refusal, before `exit 1` | `trace_decide "refused to re-mesh: the case holds results (latest time $LATEST)"` |
| `Allmesh` | `$PYTHON - <<'PY' \|\| die "dictionary generation failed"` | prefix the command: `trace_cmd dictionaries $PYTHON - <<'PY' \|\| die …` |
| `Allmesh` | the `PY` line closing that block | `trace_derivation cai2012` |
| `Allmesh` | in `if [ "$DICT_ONLY" -eq 1 ]`, before `exit 0` | `trace_decide "--dict-only: dictionaries written; the mesh was not touched"` |
| `Allmesh` | the verification `$PYTHON - <<'PY' \|\| die …` | prefix `trace_cmd verify` |
| `Allrun` | `have() { … }` | `trace_begin` |
| `Allrun` | `eval "$PLAN"` | `trace_decide "checkpoints plan: RESUME=${RESUME:-none} DECOMPOSE=$DECOMPOSE RECONSTRUCT_FIRST=$RECONSTRUCT_FIRST np=$NP"` |
| `Allrun` | in the already-at-endTime block, before `exit 0` | `trace_decide "already at endTime (latest $RESUME, endTime $END_TIME): nothing to run"` |
| `Allrun` | `echo "==> continuing from t = $RESUME; 0/ is not regenerated"` | `trace_decide "continuing from t = $RESUME: 0/ not regenerated, dsmcInitialise skipped"` |
| `Allrun` | `$PYTHON runInflow.py \|\| die …` | prefix `trace_cmd runInflow.py` |
| `Allrun` | in the `--no-solve` block, before `exit 0` | `trace_decide "--no-solve: stopped after writing 0/"` |
| `Allrun` | `DECOMP_ARGS="-decomposeParDict $OVERRIDE"` | `trace_decide "-np $NP overrides decomposeParDict ($OVERRIDE)"` |
| `Allpost` | `have() { … }` | `trace_begin` |
| `Allpost` | the two `postProcess -func … > log.… 2>&1 \|\| die …` lines and `$PYTHON postProcess.py …` | prefix `trace_cmd writeCellCentres`, `trace_cmd writeCellVolumes`, `trace_cmd postProcess.py` |
| `Allclean` | `remove log.blockMesh log.topoSet log.createPatch log.checkMesh` | append ` log.trace` to that line (`--mesh` only); no block |

- [ ] **Step 1: Write the failing tests** in `tests/unit/test_run_trace.py`. Skip on native Windows or without bash. Each test runs a script copied into `tmp_path`, using `_runner_stubs`. Define `GRAMMAR` as one regex per line kind:
  - header `^== (Allmesh|Allrun|Allpost)  \d{4}-\d\d-\d\dT\d\d:\d\d:\d\d  (.+  )?host=\S+  plumetools=\S+  OpenFOAM=\S+$`
  - end `^== end (Allmesh|Allrun|Allpost)  exit=\d+  (\d+\.\d|\?) s$`
  - step `^   step     \S+\s+(ok|FAILED\(\d+\))\s+(\d+\.\d|\?) s   \S+$`
  - decide `^   decide   .+$`
  - derive `^   derive   .+$`
  - note `^   note     .+$`

  Every line of every `log.trace` must match one of them.
  ```python
  def test_every_traced_script_embeds_the_canonical_block():
      canon = (REPO / "plumetools/trace/runner_block.sh").read_text()
      for s in TRACED_SCRIPTS: assert block_of(s) == canon, s
  def test_the_three_cai_basecases_carry_identical_scripts():   # Allmesh Allrun Allpost Allclean runInflow.py postProcess.py
  def test_a_resumed_allrun_records_the_plan_and_its_steps(tmp_path):
      # make_resume_case(... ranks=2, family="cai2012"); run Allrun -np 2
      # lines: header; "decide   checkpoints plan: RESUME=0.2 DECOMPOSE=none"; step dsmcFoam ok; step reconstructPar ok; "== end Allrun  exit=0"
  def test_allmesh_dict_only_records_the_derivation_and_the_decision(tmp_path):
      # real python, copy of cases/cai2012/baseCase; bash Allmesh --dict-only
      # contains "derive   dsmcProperties:nEquivalentParticles", "step     dictionaries", "decide   --dict-only", end exit=0
  def test_an_unwritable_trace_changes_nothing(tmp_path):
      # same resume scenario twice; second with log.trace as a directory -> identical stdout, stderr, exit status
  def test_a_failed_derivation_record_leaves_a_note_and_the_console_alone(tmp_path):
      # python shim exits 1 when argv contains "plumetools.trace"; Allmesh --dict-only stdout equals a normal run's;
      # log.trace has "note     derivation not recorded"
  def test_a_failing_step_is_recorded_and_the_console_is_unchanged(tmp_path):   # Review Focus 3
      # make_stub_bin(failing={"dsmcFoam": 3}); stderr has "'dsmcFoam' failed; full output in"; "FAILED(3)"; "exit=1"
  def test_a_second_run_appends_a_second_block(tmp_path):                       # Review Focus 4
  def test_help_and_bad_options_write_no_trace(tmp_path):                       # Review Focus 2: --help, --bogus on all three
  def test_a_case_path_with_a_space_is_traced(tmp_path):                        # Review Focus 1: tmp_path / "with space" / "case"
  def test_outside_a_checkout_the_commit_is_unknown_and_silent(tmp_path):       # Review Focus 5
      # stub `git` (exit 128) and `hostname` (exit 1) first on PATH; tmp_path is outside any checkout
      # -> header has "plumetools=unknown" and a host from `uname -n`; stderr equals a normal run's
  def test_allclean_mesh_removes_the_trace_and_plain_allclean_keeps_it(tmp_path):
  ```
- [ ] **Step 2: Run them to see them fail** (`python -m pytest tests/unit/test_run_trace.py -q`): `runner_block.sh` doesn't exist and there's no `log.trace`.
- [ ] **Step 3: Implement.** Create the block file, edit the four `cai2012/baseCase` scripts at the anchors above, then `cp` them over the `cai2012-health` and `cai2012-kn-extension` copies.
- [ ] **Step 4: Run.**
  - `python -m pytest tests/unit/test_run_trace.py tests/unit/test_allrun_restart.py -q`: expect PASS. The existing restart tests must pass unchanged.
  - Invariant check: `python tests/invariants/capture_outputs.py capture ~/rarefiedPlume-after-task8 && python tests/invariants/capture_outputs.py compare ~/rarefiedPlume-baseline-phase1 ~/rarefiedPlume-after-task8`. Expect exit 0, with no problems in the Cai scenarios.
- [ ] **Step 5: Commit:** `git commit -m "feat(trace): run trace in the Cai runner scripts"`.

### Task 9: The r-4d runner scripts

**Files:**
- Modify: `cases/r-4d/{baseCase,baseCaseHPC,baseCaseFreePlumeHPC}/{Allmesh,Allrun,Allpost,Allclean}`
- Modify: `docs/_tools/run_examples.py`: the `r-4d Allmesh --dict-only` command's `writes` gains `"cases/r-4d/baseCase/log.trace"`.
- Test: `tests/unit/test_run_trace.py` (extend `TRACED_SCRIPTS`; parametrise the Allrun and Allmesh tests over the r-4d variants)

**Interfaces:**
- Consumes: the block and `run()` edit from Task 8.

**Insertion points.** As in Task 8, with these differences:
- `Allmesh`: put `trace_begin` after the `--describe` block's closing `fi` (`exec $PYTHON -m plumetools.nozzleplume describe .`), not after `have()`, so `--describe` stays write-nothing.
- `Allmesh`: prefix `$PYTHON -m plumetools.nozzleplume setup .` with `trace_cmd setup`, and follow its `|| die …` statement with `trace_derivation nozzleplume`.
- `Allmesh`: prefix `$PYTHON -m plumetools.nozzleplume verify .` with `trace_cmd verify`.
- `Allmesh`: the `--dict-only` decide line as in Task 8.
- `Allrun`: prefix `$PYTHON -m plumetools.nozzleplume fields .` with `trace_cmd fields`.
- `Allpost`: prefix each step command that is a plain statement with `trace_cmd <name>`. Those are `postProcess` (named `writeCellCentres`/`writeCellVolumes`), `python -m plumetools.nozzleplume post` (`post`), `… diagnostics` (`diagnostics`) and the `pvpython` render (`render_slices`). Never wrap a command inside `$( … )` or a pipeline.
- `Allclean`: append ` log.trace` to the `--mesh` `remove log.blockMesh …` line.

The variants' `Allmesh`/`Allpost` differ only in comments, so the anchors exist in all three.

- [ ] **Step 1: Add the failing tests.** Add the r-4d scripts to `TRACED_SCRIPTS`, parametrise the resume, failing-step, unwritable and append tests with `family="nozzleplume"`, and add:
  ```python
  @pytest.mark.parametrize("variant", ["baseCase", "baseCaseHPC", "baseCaseFreePlumeHPC"])
  def test_describe_writes_no_trace(variant, tmp_path):                 # Review Focus 2
      case = copy_case(f"cases/r-4d/{variant}", tmp_path); run_script(case, "Allmesh", "--describe", env=...)
      assert not (case / "log.trace").exists()
  def test_r4d_allmesh_dict_only_records_the_derivation(tmp_path): ...  # "derive   dsmcProperties:nEquivalentParticles"
  ```
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Implement** the twelve script edits.
- [ ] **Step 4: Run.**
  - `python -m pytest tests/unit/test_run_trace.py -q` should PASS.
  - `python docs/_tools/run_examples.py --only "r-4d Allmesh --dict-only (bash, no OpenFOAM)"` should PASS.
  - Capture and compare against the baseline, as in Task 8: exit 0.
- [ ] **Step 5: Commit:** `git commit -m "feat(trace): run trace in the r-4d runner scripts"`.

### Task 10: Teaching case `01-cai-collisionless`

**Files:**
- Create: `cases/learn/README.md`, `cases/learn/01-cai-collisionless/case.yaml`
- Copy byte-for-byte from `cases/cai2012/baseCase/`: `Allmesh Allrun Allpost Allclean runInflow.py postProcess.py open.foam`
- Modify: `.gitignore` (after the cai2012 block: the same generated-file entries rewritten for `cases/learn/01-cai-collisionless/`, plus `cases/learn/runs/`)
- Modify: `plumetools/trace/maps/cai2012.py`: `cases` gains `"cases/learn/01-cai-collisionless"`
- Modify: `docs/_tools/run_examples.py` (the commands below)
- Modify: `tests/unit/test_run_trace.py` (`TRACED_SCRIPTS`)
- Test: `tests/unit/test_learn_cases.py`

**Interfaces:**
- Consumes: the `cai2012` map (Task 4), `plumetools.cai2012.mesh.estimate_particles(plan, geom, exit_state, n_equivalent_particles)["total_particles"]`, `plan.n_cells`, `run.n_steps`.
- Produces: `LEARN_CASES = {"cases/learn/01-cai-collisionless": ("cai2012", "cases/cai2012/baseCase")}` in the test module. Task 11 adds its case.

**`case.yaml`.** Copy `cases/cai2012/baseCase/case.yaml`, keeping every comment. Change exactly the values below, each carrying a `[TEACHING]` comment that gives the production value and the reason, for example `# [TEACHING] production 10.0 -- 3 D keeps the solver under 3 minutes`. Leave `exit.knudsen: 100.0` and every other value as is.

| key | value | key | value |
|---|---|---|---|
| `geometry.symmetry_mode` | `half_y` | `dsmc.transient_basis` | `transits` |
| `geometry.x_max_over_D` | `3.0` | `dsmc.transient_domain_transits` | `2.0` |
| `geometry.y_half_over_D` | `1.5` | `dsmc.sampling_domain_transits` | `2.0` |
| `geometry.z_half_over_D` | `1.5` | `dsmc.n_subdomains` | `1` |
| `mesh.max_core_cell_over_D` | `0.08` | `resolution.target_particles_per_cell` | `10.0` |
| `mesh.core_x_over_D` | `1.0` | `post.centerline_x_over_D_max` | `3.0` |
| `mesh.core_half_over_D` | `0.75` | `post.plane_half_over_D` | `1.5` |
| `mesh.outer_x_cells` | `8` | `post.plane_points_x` | `120` |
| `mesh.outer_lateral_cells` | `4` | `post.plane_points_lateral` | `120` |
| `mesh.max_cells` | `60000` | | |

**`run_examples.py` commands:**
- `"learn 01 Allmesh --dict-only (bash, no OpenFOAM)"`: `["bash", "cases/learn/01-cai-collisionless/Allmesh", "--dict-only"]`, needing `["cases/learn"]`, `shell=True`, `writes=("cases/learn/01-cai-collisionless/system/", "cases/learn/01-cai-collisionless/constant/", "cases/learn/01-cai-collisionless/log.trace")`, `expect_stdout=r"= 7,938"`.
- `"learn 01 a D/10 nozzle is refused"`:
  - `before=[[PY, "-c", <edit mesh.max_core_cell_over_D to 0.1 with yaml.safe_load / safe_dump>]]`
  - the same `Allmesh --dict-only`
  - `expect_exit=1`, `expect_stdout=r"FAIL nozzle_area"`
  - `writes` = `case.yaml` and `log.trace`
- `"learn 04 first leg (run_time_multiplier 0.75)"`:
  - `before` edits `dsmc.run_time_multiplier` to `0.75`
  - the same `Allmesh --dict-only`
  - `expect_stdout=r"3 in total"`
  - `writes` as in the first command, plus `case.yaml`

- [ ] **Step 1: Write the failing tests** in `tests/unit/test_learn_cases.py`.
  ```python
  @pytest.mark.parametrize("case", LEARN_CASES)
  def test_the_scripts_are_byte_identical_to_the_family_basecase(case): ...
  @pytest.mark.parametrize("case", LEARN_CASES)
  def test_the_cost_guard(case):
      cells, steps, parcels = estimate(case)        # per family, from the map's objects()
      assert cells <= 20_000 and steps <= 5_000 and parcels <= 50_000
  @pytest.mark.parametrize("case", LEARN_CASES)
  def test_every_family_check_passes(case): ...     # cai2012.checks.run_checks(...).failures == [] / nozzleplume case.checks().failures == []
  @pytest.mark.parametrize("case", LEARN_CASES)
  def test_every_value_changed_from_the_basecase_is_tagged_teaching(case):
      # leaves (section.key; a list counts as one leaf) differing between yaml.safe_load(case) and the base;
      # the line "<key>:" inside its "<section>:" block carries "[TEACHING]"
  def test_lesson1_sizes_are_the_planned_ones():
      assert estimate("cases/learn/01-cai-collisionless") == (7938, 1652, pytest.approx(2.2721e4, rel=1e-3))
  ```
- [ ] **Step 2: Run them to see them fail** (the case doesn't exist).
- [ ] **Step 3: Create the case.** Add the `.gitignore` entries, extend `MAP.cases`, add the `run_examples` commands, and add the case's scripts to `TRACED_SCRIPTS`. Regenerate pages with `python -m plumetools.trace pages`. The pages are unaffected, because the example case is `cases[0]`, but run it anyway.
- [ ] **Step 4: Run.**
  - `python -m pytest tests/unit/test_learn_cases.py tests/unit/test_trace_maps.py tests/unit/test_run_trace.py tests/unit/test_trace_pages.py -q` should PASS.
  - `python docs/_tools/run_examples.py --only "learn 01 Allmesh --dict-only (bash, no OpenFOAM)" "learn 01 a D/10 nozzle is refused" "learn 04 first leg (run_time_multiplier 0.75)"` should PASS.
  - In a copy, `cases/learn/runs/check-01`, `bash Allmesh --dict-only` should print `21 x 14 x 27 = 7,938`, `= 1652 steps`, `estimated parcels        2.2721e+04` and no `FAIL`.
- [ ] **Step 5: Commit:** `git commit -m "feat(learn): teaching case 01, a collisionless Cai plume"`.

### Task 11: Teaching case `03-own-plate`

**Files:**
- Create: `cases/learn/03-own-plate/case.yaml`
- Copy byte-for-byte from `cases/r-4d/baseCase/`: `Allmesh Allrun Allpost Allclean Allviz viz.yaml open.foam`
- Modify: `.gitignore` (the r-4d entries rewritten for `cases/learn/03-own-plate/`)
- Modify: `plumetools/trace/maps/nozzleplume.py` (`cases` gains the case), `docs/_tools/run_examples.py`, `tests/unit/test_learn_cases.py` (`LEARN_CASES`), `tests/unit/test_run_trace.py` (`TRACED_SCRIPTS`)

**Interfaces:**
- Consumes: the `nozzleplume` map; `Case.plan.estimated_cells`, `Case.run.n_steps`, `Case.run.estimated_parcels`.

**`case.yaml`.** Block style, one key per line, structured like `cases/r-4d/baseCase/case.yaml`, with its `checks`, `post` and `diagnostics` sections copied unchanged. The listing below is compact; the file is not. Every other value is below, tagged `[TEACHING] invented for this lesson` (physical values) or `[TEACHING] production <x>` (numerical values):
- `execution_profile: smoke`
- `meta: {case_name: learn-03-own-plate, description: "teaching case: a small argon thruster firing at a plate"}`
- `nozzle.exit_diameter_m: 0.05`
- `exit: {temperature_K: 300.0, speed_m_per_s: 700.0, number_density_per_m3: 6.0e+19}`
- `gas: {species_name: Ar, molecular_mass_kg: 6.63e-26, molecular_diameter_m: 4.17e-10, omega: 0.74, internal_degrees_of_freedom: 0}`
- `plate: {standoff_m: 0.25, thickness_m: 0.02, size_y_m: 0.3, size_z_m: 0.3, centre_y_m: 0.0, centre_z_m: 0.0, temperature_K: 300.0}`
- `domain: {x_max_m: 0.32, y_half_m: 0.2, z_half_m: 0.2, symmetry_planes: [y, z]}`
- `mesh: {background_cell_size_m: 0.01, plate_refinement_level: 0, refinement_regions: [{name: exitCore, level: 1, x_max_m: 0.05, radius1_m: 0.04, radius2_m: 0.05}], n_cells_between_levels: 3, max_cells: 60000}`
- `dsmc: {n_equivalent_particles: null, target_particles_per_exit_cell: 10.0, collision_model: VariableHardSphere, wall_interaction_model: MaxwellianThermal, delta_t_s: null, courant_target: 0.2, end_time_s: null, average_start_s: null, write_interval_s: null, transient_domain_transits: 3.0, sampling_domain_transits: 2.0, n_subdomains: 1, decomposition_method: scotch}`

**`run_examples.py`:** `"learn 03 describe"`, running `[PY, "-m", "plumetools.nozzleplume", "describe", "cases/learn/03-own-plate"]`, needing `["cases/learn"]`, with `expect_stdout=r"estimated cells\s+12,909"`.

- [ ] **Step 1: Add the failing tests.** Add the case to `LEARN_CASES` and add `test_lesson3_sizes_are_the_planned_ones`, asserting `(12909, 3320, approx(2.05e4, rel=0.01))`.
- [ ] **Step 2: Run them to see them fail.**
- [ ] **Step 3: Create the case** and make the other modifications.
- [ ] **Step 4: Run.**
  - The four test modules from Task 10 should PASS.
  - `python docs/_tools/run_examples.py --only "learn 03 describe"` should PASS.
  - `python -m plumetools.nozzleplume describe cases/learn/03-own-plate` should show `Knudsen number           4.4134e-01`, `20 faces`, `area error               +1.86%`, `3,320 steps` and no `FAIL`.
- [ ] **Step 5: Commit:** `git commit -m "feat(learn): teaching case 03, your own nozzle on a plate"`.

### Task 12: Move the family guides to `docs/families/`

**Files:**
- Move: `git mv docs/tutorials docs/families`
- Modify: `mkdocs.yml` (rename the nav section to `Family guides` with `families/` paths); every file that links into `docs/tutorials/`
- Modify: each family guide's at-a-glance box (add a `Trace` row)

**Interfaces:**
- Consumes: the trace pages (Task 6).

- [ ] **Step 1: Find every link.** `git grep -n "tutorials/" -- docs README.md mkdocs.yml cases/*/README.md cases/ARCHIVE.md plumetools docs/_tools`, excluding the archived case directories. Record the list.
- [ ] **Step 2: Move and relink.**
  - `git mv docs/tutorials docs/families`.
  - Rewrite each link target's `tutorials/` path component to `families/`, link text unchanged.
  - Retitle `docs/families/index.md` from `# Tutorials` to `# Family guides`, and add a first sentence: `Operational guides to each case family. New here? Start with [Learn](../learn/index.md).`
  - If a file inside an archived case directory links to a moved page, leave that file alone and add a one-line stub at the old path linking the new one. Expect none: `check_links.py` doesn't check archived directories.
- [ ] **Step 3: Add the `Trace` row** to every at-a-glance box (`| **Trace** | … |`, after `Start in`):
  - `cai2012`: `[cai2012 trace](../trace/cai2012.md)`
  - `cai2012-health` and `cai2012-kn-extension` (direct): `per-case dictionaries: [cai2012 trace](../trace/cai2012.md); study generation: not yet traced (Phase 2)`
  - `r-4d`: `[nozzleplume trace](../trace/nozzleplume.md)`
  - interfaces, scaling, virtual collisions, hybrid, markelov1999 and 3d-inflow: `not yet traced (Phase 2)`
  - ISS and particle resolution: `not yet traced (Phase 3)`
- [ ] **Step 4: Verify.**
  - `python docs/_tools/check_links.py --markdown-only` should exit 0.
  - `python docs/_tools/check_api.py` should exit 0.
  - `python docs/_tools/run_examples.py --list` should show every page example still found.
  - With the docs environment: `python -m mkdocs build --strict && python docs/_tools/check_links.py --site site` should exit 0.
- [ ] **Step 5: Commit:** `git commit -m "docs: move the per-family tutorials to Family guides"`.

### Task 13: The Learn section and the extension guide

**Files:**
- Create:
  - `docs/learn/index.md`
  - `docs/learn/01-collisionless-plume.md`
  - `docs/learn/02-trace-one-number.md`
  - `docs/learn/03-your-own-case.md`
  - `docs/learn/04-parallel-and-resume.md`
  - `docs/explanation/how-a-case-runs.md` (lesson 0)
  - `docs/guides/add-a-diagnostic.md`
- Modify:
  - `mkdocs.yml`: a `Learn` section after Getting started, listing `learn/index.md`, `How a case runs: explanation/how-a-case-runs.md` and the four lessons; plus `Add a diagnostic: guides/add-a-diagnostic.md` under Guides; plus `How a case runs` under Explanation
  - `docs/getting-started/first-run.md`: make Learn the first item in "Where to go next"
  - `README.md`: under Quick start add `New here? Start with the [Learn](docs/learn/index.md) path.`
  - `docs/_tools/run_examples.py`: register the example `learn-trace-one-number`, needing `["cases/learn", "cases/cai2012"]`

**Interfaces:**
- Consumes: `load_map`, `FamilyMap.entry`, `FamilyMap.objects`, `resolve_value` (lesson 2's example); the commands of Tasks 10–11.

**Content.** Each lesson starts with an at-a-glance box (*Goal · Case · Time · Compute · You need*), and every command block says where it runs from. Lessons work in a copy: `cp -r cases/learn/01-cai-collisionless cases/learn/runs/01` (`runs/` is gitignored). They teach this repository, not DSMC, and define a repository term (*NON-PRODUCTION*, *transit*, *parcel occupancy*) where it first appears, linking the glossary.

- **`index.md`:**
  - the path as a table: lesson, what it teaches, time, compute;
  - prerequisites: `./Allsetup`, OpenFOAM v2512 under WSL or Linux for lessons 1, 3 and 4, and lesson 2 needing no solver;
  - the compute budget statement from the Global Constraints.
- **Lesson 0** (`explanation/how-a-case-runs.md`):
  - a plain-text lifecycle diagram (`case.yaml → load → derive → assess → write dictionaries → mesh → initialise → solve → post`);
  - for each stage, what it reads and writes;
  - a table with one row per active family: the module or script per stage, linking the trace page where one exists;
  - where physics enters (tags in `case.yaml`) and where it is checked (`checks.run_checks`, the feasibility table).
- **Lesson 1:**
  1. Read `case.yaml` (every `[TEACHING]` value and why).
  2. `./Allmesh --dict-only`, reading its exit state, run settings and checks.
  3. Set `mesh.max_core_cell_over_D: 0.1`, rerun, and read the `FAIL nozzle_area` refusal (−6.89 % against 5 %), then restore 0.08.
  4. `./Allmesh`, `./Allrun --serial`, `./Allpost`.
  5. Judge `results/centerline.csv` against `plumetools.cai2012.analytical`: what agreement to expect at Kn = 100 on 7,938 cells with 2 sampling transits, where noise dominates, and that it is not a validation.
- **Lesson 2:**
  1. Open `docs/trace/cai2012.md` at `nEquivalentParticles`.
  2. Follow `producer` into `inflow.derive_run_settings`.
  3. Find the `derive   dsmcProperties:nEquivalentParticles` line in lesson 1's `log.trace` and the entry in `constant/dsmcProperties`.
  4. Predict the new weight when `resolution.target_particles_per_cell` goes from 10 to 20 (it halves).
  5. Check the prediction with `./Allmesh --dict-only`.

  It includes the page example `learn-trace-one-number`, marked with the repository's docs-example comment and registered in `EXAMPLES`. The example loads `load_map("cai2012")`, takes `m.entry("constant/dsmcProperties", "nEquivalentParticles")`, and prints its `producer` and `inputs`. It then resolves its value on `objs = m.objects(Path("cases/learn/01-cai-collisionless"))` and asserts it equals `objs["exit_state"].number_density_per_m3 * objs["plan"].exit_cell_volume_m3 / 10.0` to a relative 1e-9, so the assertion shows the formula.
- **Lesson 3:**
  1. Walk through the case's sections.
  2. `python -m plumetools.nozzleplume describe cases/learn/runs/03`.
  3. Exercise: change `plate.standoff_m` to 0.30, predict and see the cost change, then restore 0.25.
  4. `./Allmesh` (`snappyHexMesh`), `./Allrun`, `./Allpost`.
  5. Read `results/plate_summary.yaml` and the centreline.
  6. A closing section, *Making it yours*: which keys are physics, which are cost, and the checks that must pass before a production run.
- **Lesson 4:**
  1. Set `dsmc.run_time_multiplier: 0.75`, then `./Allmesh --dict-only` and `./Allrun -np 2`.
  2. Read the `decide   checkpoints plan` line.
  3. Restore `1.0`, then `./Allmesh --dict-only` and `./Allrun -np 2`, which resumes.
  4. Read the second `Allrun` block's plan (`DECOMPOSE=none`) and confirm 4 writes.
  5. Explain why 0.5 is refused (end time = averaging start).
- **Guide (`guides/add-a-diagnostic.md`):** a numbered how-to that follows the existing "Adding a diagnostic" section of `rarefaction-diagnostics.md` and the derived-quantity row of `development/extending-workflows.md`, ending with: add or update the trace-map entry, run `python -m plumetools.trace check`, then `python -m plumetools.trace pages`.

Until Task 15 measures them, the solver times in the lessons read "about 1–2 minutes on a laptop (estimate)".

- [ ] **Step 1: Write the pages and register the example.**
- [ ] **Step 2: Verify.**
  - `python docs/_tools/run_examples.py --only learn-trace-one-number` should PASS.
  - `python docs/_tools/check_links.py --markdown-only` should exit 0.
  - With the docs environment: `mkdocs build --strict` should exit 0.
- [ ] **Step 3: Commit:** `git commit -m "docs: the Learn section, lessons 0-4, and the add-a-diagnostic guide"`.

### Task 14: Reference and development documentation

**Files:**
- Modify:
  - `docs/reference/cli.md`: a `## Trace` section with `python -m plumetools.trace pages [--check]`, `check [FAMILY ...]` and `record CASE FAMILY`, with exit statuses
  - `docs/reference/outputs.md`: under `## Run bookkeeping`, `log.trace`, its line grammar verbatim, who writes each line kind, and that `Allclean --mesh` removes it
  - `docs/reference/configuration.md`: `[TEACHING]` added to the tag list
  - `docs/getting-started/repository-map.md`: `plumetools/trace/`, `cases/learn/`, `docs/trace/` and `docs/learn/` in the tree; `[TEACHING]` under conventions; `log.trace` in the generated-files table as not tracked
  - `docs/reference/case-catalog.md`: a `## Teaching cases` table with both cases, *Directory · Family · Purpose · Evidence · Runnable today · Size (cells, steps, parcels) · Lesson*
  - `docs/development/extending-workflows.md`: three rows in the "Adding …" table:
    - "a derived quantity / dictionary entry → also a trace-map entry (or boilerplate if a literal), then `python -m plumetools.trace pages`"
    - "a runner script → embed `plumetools/trace/runner_block.sh` verbatim; add it to `TRACED_SCRIPTS`"
    - "a teaching case → `cases/learn/`, scripts byte-identical to the family's `baseCase`, `[TEACHING]` tags, cost guard"
  - `docs/development/documentation.md`: the trace pages are generated and committed, the drift check, and the regenerate command
  - `docs/development/known-gaps.md`: a `## Proposed follow-on work` section with the 1-D spike from spec §9, its question and its expected output, worded as there

- [ ] **Step 1: Make the edits.**
- [ ] **Step 2: Verify.** `check_links.py --markdown-only` and `check_api.py` should exit 0. With the docs environment, `mkdocs build --strict` should exit 0.
- [ ] **Step 3: Commit:** `git commit -m "docs: reference and development pages for the trace and the teaching cases"`.

### Task 15: Solver runs on the owner's laptop (WSL, OpenFOAM v2512)

**Files:**
- Create: `tests/openfoam/test_learn_openfoam.py`
- Modify: `docs/learn/01-collisionless-plume.md`, `03-your-own-case.md`, `04-parallel-and-resume.md` (the measured times)

**Interfaces:**
- Consumes: the `openfoam` fixture and `run_or_fail` from `tests/openfoam/conftest.py`.

- [ ] **Step 1: Write the test** (marked `needs_openfoam`).
  ```python
  def test_lesson1_runs_end_to_end_within_the_budget(openfoam, tmp_path):
      case = copy_case("cases/learn/01-cai-collisionless", tmp_path)
      run_or_fail(["bash", "./Allmesh"], case)
      t0 = time.monotonic(); run_or_fail(["bash", "./Allrun", "--serial"], case); solver_s = time.monotonic() - t0
      run_or_fail(["bash", "./Allpost"], case)
      assert solver_s <= 180, f"lesson 1 solver took {solver_s:.0f} s"
      assert (case / "results" / "centerline.csv").is_file()
      trace = (case / "log.trace").read_text()
      assert "== end Allrun  exit=0" in trace and "step     dsmcFoam" in trace
  ```
- [ ] **Step 2: Run it on the laptop.** `python -m pytest -m needs_openfoam tests/openfoam/test_learn_openfoam.py -v -s` should PASS. Note the printed solver time.
- [ ] **Step 3: Run lessons 3 and 4 by hand,** exactly as the pages say, and time each `Allrun` from `log.trace`'s `== end Allrun` line. If any solver step exceeds 180 s, lower only `[TEACHING]` values, for example `target_particles_per_*` toward 8 or `sampling_domain_transits` toward 1.5. Then rerun the cost-guard and size tests, updating the planned sizes in Tasks 10–11's tests and in the case catalog, and repeat.
- [ ] **Step 4: Record the times.** Replace "about 1–2 minutes on a laptop (estimate)" in each lesson with "measured: <n> s on <CPU>, <date>".
- [ ] **Step 5: Commit:** `git commit -m "test(learn): lesson 1 end to end on OpenFOAM; measured lesson times"`.

### Task 16: Final verification and independent review

**Files:** none new.

- [ ] **Step 1: Run the full fast suite.** `python -m pytest -m "not needs_openfoam" -q -p no:cacheprovider` should show the same skips as before plus the new tests passing, and no failures.
- [ ] **Step 2: Check the trace.** `python -m plumetools.trace check && python -m plumetools.trace pages --check` should exit 0 for both.
- [ ] **Step 3: Check the invariants.** Capture into `~/rarefiedPlume-after-phase1`, then `compare ~/rarefiedPlume-baseline-phase1 ~/rarefiedPlume-after-phase1`. Expect exit 0. The teaching-case scenarios are reported only as new notes.
- [ ] **Step 4: Check the docs.** `python docs/_tools/check_api.py`, `python docs/_tools/check_links.py --markdown-only` and `python docs/_tools/run_examples.py` should all exit 0. With the docs environment, `python -m mkdocs build --strict && python docs/_tools/check_links.py --site site` should exit 0.
- [ ] **Step 5: Independent review.** Use superpowers:requesting-code-review with a reviewer that didn't write the code. Give it the spec, this plan and the branch diff. It must confirm:
  - the four invariants hold;
  - every trace-map formula matches its producer's code;
  - no boilerplate entry hides a derived value;
  - the frozen packages and archived cases are untouched (`git diff --stat main -- plumetools/cai2012 plumetools/markelov1999 cases/1d cases/2d-planar cases/2d-wedge cases/caseFoamEx cases/wake-cylinder` is empty);
  - the lessons run as written.
- [ ] **Step 6: Address the review,** then commit the fixes.

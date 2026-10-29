# Command-line reference

Every script and `python -m` entry point, from parser and source inspection.
Options are listed as the parsers define them; nothing here was invented to
fill a gap. "Writes" is what a command creates or changes; **bold** marks a
destructive effect. The *Needs* column uses: **Py** = Python with `plumetools`
(native Windows works for Python scripts invoked as `python script.py`),
**sh** = a bash shell (WSL, Linux, Git Bash), **OF** = OpenFOAM v2512 on
`PATH`, **lib** = `libplumeDsmcBoundaryModels.so`, **thermo** =
`libhybridPlumeThermo.so`, **PV** = ParaView's `pvpython`, **ff** = ffmpeg,
**CF** = CaseFoam, **PBS** = a PBS scheduler.

Every shell runner prints its own usage with `--help` (it echoes the script's
header comment and exits 0). An unknown option exits 2, except in the hybrid
`Allrun` (exits 1) and 3d-inflow's `Allmesh`, which treats anything but
`--yes` or `--dict-only` as a request to explain itself and exits 1.

## Repository setup

| Command | Start in | Options | Writes | Needs | Exit |
|---|---|---|---|---|---|
| `./Allsetup` | repository root | `--python-only` (skip OpenFOAM), `--check` (verify only; install and build nothing), `--extras a,b` (default `test,plots,cases`), `--help` | `pip install -e ".[extras]"` into the active Python; builds `libplumeDsmcBoundaryModels.so` into `$FOAM_USER_LIBBIN` | sh, Py; OF for the library | 0 all checks pass; 1 a check failed; 2 bad option |
| `applications/dsmcBoundaryModels/Allwmake` | anywhere | wmake options | `$FOAM_USER_LIBBIN/libplumeDsmcBoundaryModels.so`; build products in the source tree | sh, OF (`wmake`, `$FOAM_USER_LIBBIN`) | non-zero on failure |
| `applications/thermophysicalModels/Allwmake` | anywhere | wmake options, e.g. `-j8` | `$FOAM_USER_LIBBIN/libhybridPlumeThermo.so` | sh, OF | non-zero on failure |
| `applications/*/Allclean` | anywhere | — | **removes the build products** | sh, OF | — |

`Allsetup` does **not** build `libhybridPlumeThermo.so`; build it separately for
the hybrid family.

## Shared tools

| Command | Purpose | Options | Writes | Needs | Exit |
|---|---|---|---|---|---|
| `python -m plumetools.foamio.checkpoints report <case>` | newest complete checkpoint in the case directory and across `processor*`, and every incomplete time | — | nothing | Py | 0; 2 no such directory |
| `python -m plumetools.foamio.checkpoints finished <case>` | is the newest checkpoint complete on every rank at `endTime` (half-step slack)? | — | nothing | Py | 0 yes; 1 no; 2 no such directory |
| `python -m plumetools.foamio.checkpoints plan <case> --np N` | what `Allrun` must do; prints `RESUME= CASE_TIME= PROC_TIME= PROC_COUNT= REACHED_END= DECOMPOSE= RECONSTRUCT_FIRST=` | `--np N` (required) | nothing | Py | 0 with the plan; **3** the case must not be continued as it is (reason on stderr); 2 no such directory |
| `python -m plumetools.casemanifest [manifest.yaml]` | print every case path of a study manifest, after checking the whole manifest | `--only PATH …` (each must be in the manifest) | nothing | Py | 0; 1 invalid manifest (nothing printed) |
| `python -m plumetools.iteration <source> <destination>` | start a new mesh iteration: copy a case's inputs, nothing from its run | `--note TEXT` (kept in `iteration.yaml`), `--include PATH` (repeatable, relative to the source) | the new directory, `iteration.yaml`, `provenance/`; reads the source only | Py | 0; 1 refused |
| `python -m plumetools.rarefaction <case> --settings FILE` | mean free path, Knudsen numbers and mesh resolution for **any** reconstructed dsmcFoam case | `--time T` (default latest), `--summary FILE` (default `<case>/results/diagnostics_summary.yaml`) | diagnostic volume fields into `<case>/<time>/`; the summary | Py | 0; 1 settings or data error |

The settings file of `plumetools.rarefaction` is in the
[configuration reference](configuration.md#rarefaction-settings).

## Imagery

| Command | Purpose | Options | Writes | Needs | Exit |
|---|---|---|---|---|---|
| `pvpython --force-offscreen-rendering plumetools/viz/render_slices.py <case> [<case> …]` | slice PNGs of solver fields, as a spec lists them; several cases in one ParaView process | `--spec FILE` (default `plumetools/viz/slices.yaml`), `--field NAME`, `--plane NAME` (repeatable), `--time latest`/`first`/`all`/value, `--output DIR`, `--preset NAME`, `--field-type CELLS`/`POINTS`, `--instantaneous`, `--stop-on-error`, `--allow-constant`, `--list`, `--dry-run`, `-v` | `<case>/results/viz/` (or `--output`): PNGs, `manifest.yaml`, logs | **PV** (never plain `python`) | 0; 1 a case failed; 2 set-up or import failure. It reports its own status, so a `GLXBadContext` at exit does not turn success into failure |
| `python plumetools/viz/video.py <viz-dir> …` | encode each frame folder into a video | `--spec FILE` (its `video:` section), `--framerate F`, `--codec C`, `--container C`, `--quality N` (`-crf`, default 18), `--ffmpeg PATH` | a video beside each frame folder | Py, ff (not ParaView) | 0; 1 failure; 2 set-up error |
| `python plumetools/viz/dissect.py sheet --matrix FILE` | a matrix contact sheet | `--output` (default `contact-sheet.png`), `--cell-width` 300, `--cell-height` 780, `--trim-bottom N`, `--font F`; global `--ffmpeg` | the sheet PNG | Py, ff | 0; 1 |
| `python plumetools/viz/dissect.py converge <viz-dir> …` | image-difference curves per series | `--field` (default `dsmcRhoN`), `--metric psnr`/`ssim`, `--output` (default `convergence.yaml`) | the YAML | Py, ff | 0; 1 |

## 3d-inflow

Start in `cases/3d-inflow/`.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `./Allmesh --dict-only` | — | `system/blockMeshDict`, `snappyHexMeshDict`, `meshQualityDict` | sh, Py |
| `./Allmesh --yes` | with no argument it prints a warning and exits 1 | **overwrites `constant/polyMesh`; deletes `constant/polyMesh/sets/` and snappyHexMesh bookkeeping**; logs | sh, Py, OF |
| `./Allrun` | `--serial`, `-np N`, `--no-solve` | `0/`, logs, time directories; `system/decomposeParDict.Allrun` for `-np` | sh, Py; OF unless `--no-solve`; always a mesh |
| `./Allclean` | `--dry-run`, `--mesh` (also `system/blockMeshDict`) | **removes `0/`, numeric time directories, `processor*`, `log.*`, `postProcessing`, `boundaries`, `fieldMeasurements`** | sh |
| `python runInflow.py [case_dir]` | — | `0/boundaryU`, `boundaryT`, `boundaryNumberDensity_<sp>`; nine measurement fields (standard dialect) | Py, a mesh |
| `python -m plumetools.verify_mesh <case>` | — | nothing; exit 1 if the mesh is not the declared geometry | Py, a mesh |

## cai2012

Start in `cases/cai2012/`.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `./generate_cases.py` | `--dry-run` (list and stop), `--knudsen KN …`, `--particles LEVEL …` (`2` or `np2x`), `--overwrite`, `--study PATH` | `Cases/<Kn>_<level>/`, `manifest.yaml`; **`--overwrite` replaces cases with results** | Py |
| `./validate_cases.py` | `-v`/`--verbose`, `--manifest PATH` | nothing; exit 0 the tree is the study, 1 not | Py, generated dictionaries |
| `./AllmeshCases` | `--dict-only`, `--keep-going`, `--list`, `--index N`, `--case NAME` | per case, as `Allmesh` | sh, Py; OF unless `--dict-only` |
| `./AllrunCases` | `--serial`, `-np N`, `--no-solve`, `--keep-going`, `--list`, `--index N`, `--case NAME` | per case, as `Allrun` | sh, Py, OF, lib |
| `./AllpostCases` | `--keep-going`, `--no-viz`, `--viz-only`, `--viz-spec F` | per-case `results/`; `results/study-table.csv`; imagery | sh, Py, OF (`postProcess`); PV optional |
| `./AllcleanCases` | `--dry-run`, `--mesh`, **`--purge` (deletes `Cases/` and `manifest.yaml`)** | **removes solver output per case** | sh, Py |
| `Cases/<c>/Allmesh` | `--dict-only`, `--force` (re-mesh a case with results) | dictionaries; mesh | sh, Py; OF unless `--dict-only` |
| `Cases/<c>/Allrun` | `--serial`, `-np N`, `--no-solve` | `0/`, time directories, logs; resumes | sh, Py, OF, lib |
| `Cases/<c>/Allpost` | `--time T` | `results/` | sh, Py, OF |
| `Cases/<c>/Allclean` | `--dry-run`, `--mesh` | **removes solver output and the generated `0/`** (the mesh too with `--mesh`) | sh |
| `Cases/<c>/postProcess.py` | `--time T`, `--mesh-time T`, `--results DIR` | `results/` | Py |

## cai2012-health

Start in `cases/cai2012-health/`.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `./generate_cases.py` | `--dry-run` (cost and destinations; writes nothing; no CaseFoam), `--weight NAME …`, `--sampling NAME …`, `--overwrite`, `--study PATH` | `Cases/<weight>/<sampling>/`, `manifest.yaml` | Py; CF unless `--dry-run` |
| `./AllmeshCases` | `--dict-only`, `--keep-going` | as cai2012 | sh, Py, OF |
| `./AllrunCases` | `--only PATH` (repeatable), `--serial`, `-np N`, `--no-solve`, `--keep-clouds`, `--keep-going` | runs; `results/cost-model.yaml`; **prunes reconstructed `lagrangian/` clouds from all but the newest time** unless `--keep-clouds` | sh, Py, OF, lib |
| `./AllpostCases` | `--no-viz`, `--viz-only`, `--analyse-only`, `--no-audit`, `--no-analyse`, `--stride N`, `--viz-spec F`, `--keep-going` | per-case results; `results/` sheets and tables | sh, Py, OF; ff, PV |
| `./AllcleanCases` | `--dry-run`, `--mesh`, **`--purge`** | as cai2012 | sh, Py |
| `python analyse.py` | `--field` (default `dsmcRhoN`), `--no-sheets`, `--no-overlap`, `--ffmpeg PATH` | `results/contact-*.png`, `convergence.yaml`, `sweep-table.csv`, `sweep-error.png`, `overlap.yaml`; exit non-zero if the overlap check fails | Py; ff for sheets |
| `Cases/…/Allaudit.py` | `--latest-only`, `--stride N` | `results/occupancy.yaml`, `results/convergence.yaml` | Py |

## cai2012-kn-extension

Start in `cases/cai2012-kn-extension/`.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `python generate_cases.py` | `--dry-run`, `--refresh`, `--regenerate`, **`--overwrite`**, `--study PATH` — semantics in the [tutorial](../tutorials/cai2012-kn-extension.md#generation-exactly) | the active cases tree and manifest; a temporary `.staging-*` beside `study.yaml` (removed on exit, also under `--dry-run`) | Py |
| `python validate_cases.py` | `-v`, `--inputs-only` (skip checks that read generated dictionaries), `--manifest PATH` | nothing | Py |
| `python smoke_cases.py` | `--study PATH`, `--report PATH` | **generates** the tree, writes every case's dictionaries and `results/full-suite-preparation/smoke.json`; refuses a tree holding a mesh or solver output | Py |
| `python -m plumetools.cai2012ext.policy [manifest] --study PATH [--all]` | list schedulable case paths; `--all` includes non-production ones | nothing | Py; exit 2 if the manifest breaks the suite policy |
| `python generate_interfaces.py` | `--dry-run`, `--report`, `--refresh`, `--regenerate`, **`--overwrite`**, `--study`, `--interfaces` | `Interface/`, `interface-manifest.yaml`; per-case `results/interface-report.yaml` with `--report` | Py |
| `python -m plumetools.cai2012ext.interface_case STEP <case>` | steps `report` (`--decision-table`), `mesh-dicts`, `verify-mesh`, `inflow`, `check-decomposed`, `collect` (`--output FILE`) | the step's report, dictionaries, identity or inflow | Py (a mesh for `verify-mesh`, `inflow`, `check-decomposed`); `report` exits 3 if not meshable, `mesh-dicts` 3 if refused |
| `python -m plumetools.cai2012ext.overlap <case> [--baseline PATH]` (= `Alloverlap.py`) | override `overlap.baseline_case` (read-only) | `results/overlap.yaml` | Py; exit 1 on error |
| `python -m plumetools.cai2012ext.hpc STEP <case>` | steps `guard`, `label`, `finished`, `record <seconds>` | `label`: a `NON-PRODUCTION` marker; `record`: `results/cost-model.yaml` under an `fcntl` lock (Unix only) | Py; `guard` exit 2 refused; `label` exit 3 refused (interface cases); `finished` 0/1 |
| `python Allscale.py MODE` | modes `inventory`, `prepare`, `fit`, `report`, `evaluate`, `smoke`; `--config`, `--output`, `--full` (required by `fit`), `--no-plots`, `--predictor`, `--points`, `--kn`, `--sensitivities` | under `results/scaling/full-suite/` by default; never runs DSMC | Py; exit 2 on error |
| `./AllmeshCases` | `--dict-only`, `--keep-going`, `--list`, `--index N`, `--case NAME` | as cai2012 | sh, Py, OF |
| `./AllpostCases` | `--case NAME` (repeatable), `--no-viz`, `--viz-only`, `--no-analyse`, `--analyse-only`, `--no-diagnostics`, `--stride N`, `--keep-going` | per-case results, imagery, `results/full-suite/` | sh, Py, OF; PV |
| `./AllcleanCases` | `--dry-run`, `--mesh`, **`--purge` (the active tree and its manifest)** | | sh, Py |
| `python analyse.py` | — | `results/full-suite/study-table.csv`, `study-table.md`, `knudsen-comparison.png` | Py |
| `<case>/Alldiagnose.py` | `--time T`, `--no-fields` | diagnostic fields and summaries | Py |
| `<case>/Allcompare.py` | `--time`, `--mesh-time`, `--radius` (r/D, default 10), `--shell` (0.5), `--bin` (2.5°), `--no-plots` | `results/comparison.yaml`, `centerline_extended.csv`, `meridional_plane.csv` | Py |
| `<case>/Allconverge.py` | `--stride N`, `--tolerance F`, `--no-plots` | `results/stationarity.yaml` | Py |
| `qsub runMeshAll.sh`, `qsub runAnalyse.sh`, `qsub runScaling.sh` | study level | as the steps they run | PBS + a site `hpc_env.sh` |
| `cd <case> && qsub runCase.sh` (also `runMesh.sh`, `runPost.sh`, `runViz.sh`) | per case | `output.log`, `error.log`, as the steps they run | PBS + a site `hpc_env.sh` |

## cai2012-virtual-collision

Start in the repository root.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `python cases/cai2012-virtual-collision/assess.py MODE` | modes `inventory`, `smoke`, `full`; `--root`, `--config`, `--output` (default the case's `results/` — **`results/smoke/` is tracked**), `--no-plots` | inventory, report, model arrays and assessment files into `--output` | Py |
| `python cases/cai2012-virtual-collision/compare_exports.py` | `--study`, `--samples` (2048), `--hazard-steps` (24), `--grid-refinement 0`/`1`/`2`, `--kn …`, `--no-plots` | into the Kn-extension cases' `results/virtual_collision/` | Py |
| `python cases/cai2012-virtual-collision/summarize_exports.py` | — | refreshes the study index | Py |
| `python cases/cai2012-virtual-collision/audit_hazard.py` | `--samples`, `--seed`, `--steps …`, `--kn …`, `--output` | a hazard audit | Py |
| `python cases/cai2012-virtual-collision/review_results.py` | — | `…/results/virtual_collision/review/*_review.{json,csv}` | Py |
| `qsub runAssess.sh` (from the case directory) | `-v VC_HPC_ENV=…` | `results-hpc/` | PBS |

## cai2012-hybrid

Start in the repository root.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `cases/cai2012-hybrid/Allrun <case> <stage>` | stages `setup` (needs `--profile`), `cfd`, `post-cfd`, `dsmc`, `post-dsmc`, `crossings`, `all`; `--profile NAME`, `--config FILE`, `-np N`, `--subdir DIR`, `--allow-unsuitable` | the case directory; resumes; never deletes | sh, Py; OF + thermo (`cfd`), OF + lib (`dsmc`) |
| `python -m plumetools.hybridplume setup <case> --config FILE` | `--profile`, `--np`, `--background-ratio R` | `hybrid.resolved.yaml`, `results/setup.json`, `cfd/` (or `cfd_bg…/`) | Py |
| `python -m plumetools.hybridplume STEP <case>` | `check-mesh` (`--subdir`), `cfd-report` and `breakdown` (`--time`), `interface` (`--time`, `--allow-unsuitable`), `dsmc-setup` (`--np`), `dsmc-fields`, `dsmc-post` (`--time`), `crossings`, `overlap` (`--reference`, required), `reference`, `report`, `hpc` (`--scheduler slurm`/`pbs`/`none`, default `slurm`; `--np`) | `results/<step>.json` and the step's files | Py; exit 2 when a step cannot proceed |

## markelov1999

Start in `cases/markelov1999/`.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `./generate_cases.py` | `--dry-run`, `--pressures PSI …`, `--overwrite`, `--study PATH` | `Cases/gap06in/<case>/` (CaseFoam's tree, including `Cases/baseCase`), `manifest.yaml`; CaseFoam also writes `Allrun`, `Allclean`, `rmCases` here | Py, **CF (also for `--dry-run`)** |
| `./AllmeshCases` | `--dict-only`, `--keep-going` | per case | sh, Py, OF |
| `./AllrunCases` | `--serial`, `-np N`, `--no-solve`, `--keep-going` | per case, sequentially | sh, Py, OF, lib |
| `./AllpostCases` | `--keep-going`, `--no-viz`, `--viz-only`, `--viz-spec F` | `results/study-table.csv`, pressure plots | sh, Py |
| `./AllcleanCases` | `--dry-run`, `--mesh`, **`--purge`** | | sh, Py |
| per case `Allmesh`, `Allrun`, `Allpost`, `Allclean` | as cai2012 (`Allpost` has no options) | | |
| `python -m plumetools.markelov1999.verify <case>` | — | nothing; exit 1 on a geometry mismatch | Py, a mesh |

## r-4d

Start in `cases/r-4d/<variant>/`.

| Command | Options | Writes | Needs |
|---|---|---|---|
| `./Allmesh` | `--describe` (= `python -m plumetools.nozzleplume describe .`; writes nothing), `--dict-only`, `--force` | dictionaries; mesh | sh, Py; OF unless `--describe`/`--dict-only` |
| `./Allrun` | `--serial`, `-np N`, `--no-solve` | `0/`, time directories; resumes | sh, Py, OF, lib |
| `./Allpost` | `--time T`, `--no-viz`, `--viz` (a missing `pvpython` is an error), `--no-diagnostics` | `results/`, diagnostic fields in `<time>/`, `results/viz/` | sh, Py, OF; PV optional |
| `./Allviz` | `--field`, `--plane` (repeatable), `--time T`/`all`, `--dry-run`, `--list`, `--output DIR`, `--spec FILE` | `results/viz/` only | sh, PV |
| `./Allclean` | `--dry-run`, `--mesh` | **removes solver output** (and the mesh) | sh |
| `python -m plumetools.nozzleplume COMMAND [case]` | `describe`, `setup`, `verify`, `fields` (`--no-checks`), `post` (`--time`, `--results`), `diagnostics` (`--time`, `--results`) | `describe`: nothing; the others as named | Py (a mesh for `verify`, `fields`; a run for `post`, `diagnostics`) |

## Legacy cases

`cases/iss_solar_panels*` and `cases/solar_panel_particle_resolution_study/Fnum_*`
have no supported commands; their tutorials say why. **Do not run their
`clean.sh`.**

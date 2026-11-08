# Run the adaptive radial sizing study

Use WSL/Linux with Python and the installed OpenFOAM v2512 hyStrath port. The
owner's setup below is an example:

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
export PATH="$HOME/OpenFOAM/andy-v2512/platforms/linux64GccDPInt32Opt/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/OpenFOAM/andy-v2512/platforms/linux64GccDPInt32Opt/lib:${LD_LIBRARY_PATH:-}"
cd /path/to/rarefiedPlume
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
export PYTHON=/path/to/python-with-numpy-scipy-pyyaml
"$PYTHON" cases/cai2012-radial-resolution/generate_cases.py --dry-run
```

Inspect [study.yaml](../../cases/cai2012-radial-resolution/study.yaml). Defaults
use a 2.5-degree half-angle, outlet plane at 10 D_e, and 100 uniform radial cells.
Generation reports predicted requirements and costs. The initial mesh is
exploratory and may miss the final requirements. It writes dictionaries only,
into the fresh `Cases-adaptive` tree:

```bash
"$PYTHON" cases/cai2012-radial-resolution/generate_cases.py
cd cases/cai2012-radial-resolution
bash ./Allrun --dry-run
bash ./Allrun
```

`Allrun` queues Kn 100 down to 0.0001 with bounded case-level parallelism. Each
solver uses one core. The default is min(7, available CPUs minus1), at least1;
use `bash Allrun --jobs 4` to choose a limit, or `--jobs 1` for serial execution.
Allrestart accepts the same option. Failure stops further launches while active
cases finish safely. Shared summary writes are protected. It builds/verifies meshes, then
runs short native segments ending at complete writes. After startup, four
post-startup windows supply two blocks for a per-cell 5% n/U/T change check.
Stable measurements allow assessment of edge/MFP, occupancy and timestep
requirements. Missing or unstable evidence continues the current run.

Passing cases finish early. Trusted failures prepare `r1_ppc20`, then later
revisions if needed, record a compact summary and correction plan, clean the
superseded attempt's native outputs, and initialize the corrected case afresh.
Population corrections aim for 25 parcels/cell; acceptance stays at 20. Five
retries are allowed by default. Inlet, angle and outlet remain fixed.

Watch the dedicated summary, which is never truncated by cleanup or retries:

```bash
tail -f Cases-adaptive/adaptive-study.log
```

Every assessment records Kn, attempt, simulated time, wall time, charged steps,
cells, particle weight, timestep, minimum mean parcels, maximum edge/MFP,
failed/unresolved counts, and action/reason. Retry records include old/new
settings. Each Kn directory has `controller-state.json` for recovery; attempts
retain dictionaries, settings, mesh metadata and compact assessment/plan records.
Accepted outputs remain available for postprocessing.

To continue after an interruption:

```bash
bash ./Allrestart --dry-run
bash ./Allrestart
```

`Allrestart` skips verified accepted cases, resumes complete checkpoints on an
unchanged attempt, recovers interrupted retry transitions, and starts untouched
cases. It refuses incomplete latest writes or unrelated input/mesh edits with
recovery instructions. Do not remove controller state or edit running inputs.
Known temporary endTime overrides are restored exactly; arbitrary edits are
refused. Solver failures retain their outputs and logs.

Cell, parcel, cumulative-step, sampling and retry limits can prevent a pass.
Both scripts stop on such an outcome or a solver failure. Inspect the summary
before explicitly using `--allow-budget-exceedance` for step/parcel limits;
this flag never overrides retries or the mesh cell cap. A wider fixed angle can
still require many cells at low Kn. Unresolved/budget-limited results never pass.

Both scripts accept `--cases /path/to/OTHER_TREE`. A dry run needs Python only;
it creates no files and launches no solver. Earlier `Cases` runs are preserved.
To continue their legacy fixed-duration workflow:

```bash
bash ./Allrestart --cases "$PWD/Cases"
```

To disable automatic control in a fresh study, set `controller.enabled: false`.
Set `mesh.initial_cells: null` for the earlier predicted-target initial mesh.
Use a new tree whenever changing geometry or sampling policy. Matching
generation is idempotent; incompatible existing inputs are refused.

An individual generated case's `Allrun` follows its embedded controller policy;
`bash ./Allrun --resume` continues it. Use the accepted attempt's own `Allpost`
for plots/diagnostics, and `python analyse.py` at family level for a study report.
Numerical sizing uses recent stable windows. Reports label this basis and do
not certify statistical independence or flow convergence; legacy full-run
analysis retains its original block checks.

Read the [method](../cai2012-radial-resolution.md) and the
[portable unrelated-case workflow](../mesh-refinement.md) before interpreting
streamtube totals as guidance for another geometry.

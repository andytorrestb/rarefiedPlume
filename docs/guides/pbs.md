# PBS jobs

The repository's batch scripts, what they assume, and what to change before
using them on your cluster. **No scheduler job was submitted for this
documentation**; every statement below comes from reading the scripts.

## The scripts

| Script | Family | `#PBS` request as committed | Does |
|---|---|---|---|
| `cases/cai2012-kn-extension/runMeshAll.sh` | Kn extension | `shortq`, 4 h, `nodes=1:ppn=1` | `generate_cases.py`, `validate_cases.py --inputs-only`, `AllmeshCases`, `validate_cases.py` |
| `<case>/runMesh.sh` | Kn extension | `shortq`, 2 h, 1 core | `hpc guard`, then `./Allmesh` |
| `<case>/runCase.sh` | Kn extension | `longq`, 120 h, `nodes=1:ppn=192` | `hpc guard`, `hpc label`, `hpc finished`, `./Allrun`, `hpc record` |
| `<case>/runPost.sh` | Kn extension | `longq`, 120 h, 192 cores | `AllpostCases --case <case> --no-viz --no-analyse` (interface cases: `Allpost`, `Alldiagnose.py`, `Allcompare.py`, `Allconverge.py`, `Alloverlap.py`, `collect`) |
| `<case>/runViz.sh` | Kn extension | `longq`, 120 h, 192 cores (interface: refuses) | `AllpostCases --case <case> --viz-only` |
| `cases/cai2012-kn-extension/runAnalyse.sh` | Kn extension | `shortq`, 1 h, 1 core | `AllpostCases --analyse-only` |
| `cases/cai2012-kn-extension/runScaling.sh` | Kn extension | `shortq`, 1 h, 8 cores | `Allscale.py inventory`, `fit --full`, `report` |
| `cases/cai2012-virtual-collision/runAssess.sh` | virtual collisions | `shortq`, 2 h, 8 cores, 24 GB | `assess.py full --output results-hpc` |
| `<case>/hpc/run_<stage>.pbs` (generated) | hybrid | `longq`, `nodes=1:ppn=192`, walltime 3 × estimate (1–120 h) | one `Allrun` stage |

No other family ships batch scripts. For `cai2012`, `cai2012-health`,
`markelov1999` and `r-4d`, run `Allrun` (or the study runners) inside a job
script of your own; the study runners' `--list` / `--index N` / `--case NAME`
options exist for job arrays.

## What is site-specific

Every job script sources `cases/cai2012-kn-extension/hpc_env.sh` (found in
`./`, `../` or `../../` of `$PBS_O_WORKDIR`; the virtual-collision job takes
`VC_HPC_ENV`, the hybrid scripts `HYBRID_HPC_ENV`). That file is the owner's
cluster: a Spack set-up path, Lmod, specific module names and versions for
Python, Open MPI and ffmpeg, an OpenFOAM v2512 build under a home directory, a
Python virtual environment, and a ParaView path. The `#PBS` headers name the
owner's queues (`shortq`, `longq`), a 192-core node, and a notification
address (`-M`, `-m abe`). None of these will work elsewhere unchanged.

A minimal replacement keeps the file's contract — it is **sourced, not
executed**, must not be sourced under `set -u`, and leaves OpenFOAM,
`plumetools` and (optionally) `pvpython` usable:

```bash
# hpc_env.sh -- placeholders: replace every <...>
module purge
module load <your-python-module> <your-mpi-module>          # whatever your site calls them
set --                                                      # hand OpenFOAM's bashrc no arguments
source <path-to-OpenFOAM-v2512>/etc/bashrc
source <path-to-your-checkout>/.venv/bin/activate           # with pip install -e ".[test,plots]"
export PVPYTHON_PATH=<path-to-pvpython>                     # optional: imagery
export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow    # headless rendering, as the owner's file sets
```

and the job headers need your queue, node shape, walltime and mail settings:

```bash
#PBS -q <queue>
#PBS -l walltime=<hh:mm:ss>
#PBS -l nodes=1:ppn=<cores-per-node>
#PBS -M <your-address>
```

Keep `dsmc.n_subdomains` in `case.yaml` equal to the cores the job gets:
`runCase.sh` does not override it, and `Allrun` reads it from the generated
`decomposeParDict`.

## Order of submission

The scripts do not chain jobs and nothing resubmits itself. Submit each stage
after the previous one has finished:

1. **Build** the libraries on the cluster (`applications/*/Allwmake` after
   sourcing your `hpc_env.sh`), interactively or in a short job.
2. **Mesh**: `qsub runMeshAll.sh` (study) or `qsub runMesh.sh` (per case).
3. **Validate**: `python validate_cases.py` and read the feasibility warnings.
4. **Run**: one `runCase.sh` per case, from the case directory
   (`for case in $(python -m plumetools.cai2012ext.policy --study study.yaml --all); do (cd "$case" && qsub runCase.sh); done`).
5. **Post-process**: `runPost.sh` per case, then `runAnalyse.sh`.

PBS can express dependencies (`qsub -W depend=afterok:<jobid> …`), which would
let steps 2–5 be queued at once; the repository's scripts do not use it and it
has not been tested here.

## When a job ends early

A run killed by its walltime is **resubmitted by hand with the same `qsub`**:
`./Allrun` continues from the newest checkpoint complete on every rank, and a
case already at `endTime` is detected (`hpc finished`) and not re-timed. If
the kill happened while a time was being written, `Allrun` refuses to continue
and names the incomplete time; see
[Restart and recovery](restart-and-recovery.md#when-a-run-is-refused).

`runCase.sh` writes `output.log` and `error.log` in the case directory. The
first lines of `output.log` are `hpc label`'s verdict, including the
`NON-PRODUCTION` text for an unsuitable case. `hpc record` appends the wall
clock to the family's `results/cost-model.yaml` under a file lock, so two jobs
finishing together do not lose a measurement.

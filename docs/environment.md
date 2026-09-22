# Environment

## The honest summary

**The exact software stack that produced the historical results cannot be
reconstructed.** No dependency file, lockfile, or environment specification was
ever committed. What follows is everything recoverable from committed artifacts,
plus what to do going forward.

Bit-level reproduction of archived output is therefore **not achievable**.
Behavioural reproduction may be.

## What is recoverable

| Component | Version | Evidence |
|---|---|---|
| OpenFOAM | **v1706** | `cases/wake-cylinder/5psi/sample_out.txt`, line 8: `Build : v1706` |
| Solver | `dsmcFoam+` / `dsmcInitialise+` — the micro/nano-flow (MNF) fork | invoked in `run.sh` and `runCases.py`; **version not recorded anywhere** |
| Meshing | Pointwise **V18.5R2** | `.pw` file headers (`PWI0`, "Pointwise V18.5R2") |
| Post-processing | ParaView **5.10.0** | `util/paraView/trace_lineplot.py`, line 1 |
| MATLAB | unknown | `cases/2d-wedge/mesh/generateBlockMeshDict.m` |
| gnuplot | unknown; needs an X11 terminal | `monitor` scripts (`set term x11`) |

Note the `controlDict` headers claim **both** v1706 and v2106 across different
cases, so they are not reliable evidence.

## What is not recoverable

Python, numpy, scipy, pandas, matplotlib, `casefoam` and `fluidfoam` versions.
None was ever recorded. The only constraint that can be inferred is
`scipy >= 1.6`, because `integrate.trapezoid` was added there (earlier versions
spell it `trapz`).

## Current requirements

> **The MNF fork is not required.** The active cases target OpenFOAM's own
> `dsmcFoam`, verified against **v2512**. `dsmcFoam+` appears only in the
> archived cases and in the provenance table above. If something asks you for a
> library, that is the custom inflow model in step 3 below, not the fork.

### Setting up a new machine

```bash
./Allsetup                  # all three steps below, then verifies the result
./Allsetup --python-only    # skip everything needing OpenFOAM
./Allsetup --check          # verify only; install nothing, build nothing
```

Three things have to be in place, and `pip` does only the first:

| | | |
|---|---|---|
| 1 | the `plumetools` package | `pip install -e ".[test]"` |
| 2 | standard OpenFOAM on PATH | `dsmcFoam`, `dsmcInitialise`, `blockMesh`, … |
| 3 | `libplumeDsmcBoundaryModels.so` | `applications/dsmcBoundaryModels/Allwmake` |

**Step 3 is the one that gets missed.** Nothing pip does touches it, and a case
that needs it fails only at `./Allrun` — several steps into a study, with a
message that reads like a missing solver rather than a missing build step. Every
`cases/markelov1999` case selects `plumeFieldInflow` in `constant/dsmcProperties`
and refuses to run without the library; there is deliberately no fallback to a
uniform inflow, because that would silently discard the plume's angular structure.

### The Python library — no OpenFOAM needed

```bash
pip install -e ".[test]"
pytest -m "not needs_openfoam"
```

Pinned in `pyproject.toml`: `numpy >= 1.20`, `scipy >= 1.6`, `pyyaml >= 5.4`;
`pytest >= 7` for the test extra. Python ≥ 3.9.

The full Tier 0 + regression suite runs in about two seconds on Windows and
Linux, with no OpenFOAM installed. This is deliberate: mesh parsing, the
source-flow model, and the field writers are all testable without a solver.

### Running a case — needs OpenFOAM

`dsmcFoam` and `dsmcInitialise` come from the standard distribution. Other
standard utilities used: `blockMesh`, `snappyHexMesh`, `topoSet`, `checkMesh`,
`decomposePar`, `reconstructPar`, `postProcess`, and `mpirun` for parallel runs.

Plus `libplumeDsmcBoundaryModels.so` from step 3 — see
[`plume-field-inflow.md`](plume-field-inflow.md) for what it does and
[`solver-compatibility.md`](solver-compatibility.md) for why standard `dsmcFoam`
needs it at all.

The archived cases (`1d`, `2d-planar`, `2d-wedge`, `caseFoamEx`, `wake-cylinder`)
still name `dsmcFoam+` in their `controlDict`s, as do `util/runCases.py` and
`util/caseFoam/runCases.py`. Those are frozen historical record and would need
the fork; see `cases/ARCHIVE.md`.

Tests requiring the toolchain are marked `needs_openfoam` and are **not run by
default**:

```bash
pytest -m needs_openfoam     # only where dsmcFoam and the library are available
```

### Post-processing — needs ParaView

`paraview.simple` ships with ParaView and is **not pip-installable**. The
one-offs under `util/` import it directly and are not part of the library.

> **`*.csv` is gitignored**, and no CSV file exists anywhere in the repository —
> yet six plotting scripts require one as input (finding RP-06). The Fnum
> resolution study cannot currently be reproduced end to end.

#### Slice imagery — `plumetools.viz` and VifPara

`plumetools/viz` renders slice PNGs of any solver field through
[VifPara](https://github.com/virtual-vehicle/VifPara), which automates ParaView.
Three things have to line up, and `pip` does only the third:

| | | |
|---|---|---|
| 1 | **ParaView 5.11–5.13** installed and on `PATH` | `pvpython`, `paraview` |
| 2 | ParaView's *internal* Python minor version **==** the venv's | see below |
| 3 | the `viz` extra | `pip install -e ".[viz]"` |
| 4 | **ffmpeg**, for the time-series videos only | `apt install ffmpeg` |

Step 4 is optional and is *not* pip-installable. Without it the frames are still
written and the step prints one line saying no videos were made; encode them
later with `python plumetools/viz/video.py <viz-dir>`, which needs ffmpeg but not
ParaView. Verified against ffmpeg 6.1.1 with `libx264`.

**Step 2 is the one that bites.** VifPara layers this virtual environment onto
ParaView's interpreter, so the two Pythons must match. Check both:

```bash
pvpython -c "import sys; print(sys.version)"     # ParaView's
python    -c "import sys; print(sys.version)"    # the venv's
```

Verified working here: ParaView **5.11.2**, both Pythons **3.12.3**. ParaView 6
is not supported by VifPara.

If `pvpython` is not on `PATH`, point at it explicitly:

```bash
export PVPYTHON_PATH=/path/to/pvpython
```

Then render, **through the `vifpara` launcher, never `python`**:

```bash
vifpara plumetools/viz/render_slices.py cases/cai2012/Cases/Kn100_np1x
```

The launcher re-runs the script inside `pvpython` with `PV_VENV` set to this
environment. Running it with `python` fails at `import paraview.simple`, because
that module exists only inside ParaView's own interpreter.

> **The `viz` extra downgrades numpy.** `vifpara` 1.3.4 pins `numpy==1.26.4`
> exactly, which also pulls `scipy` back (1.18 → 1.17 here). The full suite —
> 956 tests — passes on both, and this was checked rather than assumed, but it is
> why `viz` is an optional extra and not a core dependency. Install it into the
> environment you actually render from.

On a headless or WSLg display, `pvpython` dies at exit with `GLXBadContext`
tearing down its GL contexts, *after* every image is written. `render_slices.py`
therefore reports its own exit status rather than letting that turn a successful
render into a failure; see [`viz-slices.md`](viz-slices.md).

Full guide: [`viz-slices.md`](viz-slices.md).

## Windows notes

The library and test suite are fully supported on Windows. Two rough edges:

1. **Line endings.** The legacy code wrote field files with `print()` to a
   text-mode file, so it emitted CRLF on Windows and LF on Linux for identical
   inputs. `plumetools` writes `newline="\n"` explicitly, so its output no longer
   depends on the platform. The regression test normalises line endings before
   comparing.

2. **A symlink in an archived case.** `cases/caseFoamEx/2d-wedge/system/blockMeshDict`
   is the repository's only git symlink. Without symlink privileges it
   materialises as a text file containing the path, leaving the working tree
   permanently dirty and making `git diff` on it fail with
   `Function not implemented`. Fix before cloning:

   ```bash
   git config --global core.symlinks true
   ```

   Not repaired in-place, because that would mean modifying an archived case
   (finding AR-14).

The case shell scripts (`Allrun`, `Allmesh`, `Allclean`) are `#!/bin/bash` and
assume a POSIX shell; on Windows use WSL or Git Bash.

## Capturing provenance going forward

Since the original environment is unrecoverable, the fix is to stop the same
thing happening again. Each run should record, alongside its output:

- the resolved `case.yaml`
- `git rev-parse HEAD`
- `foamVersion` (or `$WM_PROJECT_VERSION`)
- `pip freeze`
- a timestamp

`plumetools.run` is the intended home for writing that as `<case>/provenance.txt`.
It is not implemented yet — the run wrappers are still a later stage.

## Version-pinning recommendation

For a reproducible environment, pin exactly rather than by lower bound:

- a `requirements.lock` from `pip freeze` in a known-good environment, committed
  alongside `pyproject.toml`'s looser ranges;
- a container image pinning OpenFOAM v2512 and building
  `applications/dsmcBoundaryModels` into it, which pins the whole active stack —
  the standard distribution is publicly available, so nothing here is blocked.

Reproducing the **archived** results is a separate and harder problem: it would
need OpenFOAM v1706 plus the MNF `dsmcFoam+` build, and no source for that build
has been located. That remains the single largest reproducibility gap, and it is
not solvable from inside this repository. It does not affect the active cases,
which no longer depend on the fork.

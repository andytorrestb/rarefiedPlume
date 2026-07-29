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

`dsmcFoam+` and `dsmcInitialise+` come from the MNF fork of OpenFOAM, not from
the standard distribution. Standard utilities used: `blockMesh`, `topoSet`,
`checkMesh`, `decomposePar`, `reconstructPar`, `postProcess`, and `mpirun` for
parallel runs.

Tests requiring them are marked `needs_openfoam` and are **not run by default**:

```bash
pytest -m needs_openfoam     # only where dsmcFoam+ is on PATH
```

### Post-processing — needs ParaView

`paraview.simple` ships with ParaView and is not pip-installable. Scripts that
import it are one-offs under `util/` and are not part of the library.

> **`*.csv` is gitignored**, and no CSV file exists anywhere in the repository —
> yet six plotting scripts require one as input (finding RP-06). The Fnum
> resolution study cannot currently be reproduced end to end.

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
- a container image pinning OpenFOAM v1706 plus the MNF `dsmcFoam+` build, if a
  source for that build can be obtained — this is the single largest remaining
  reproducibility gap, and it is not solvable from inside this repository.

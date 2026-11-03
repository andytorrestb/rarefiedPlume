# Installation

Four independent pieces, each needed only by some steps. Install what the
work in front of you needs; `./Allsetup --check` reports which pieces are
present.

| Piece | Needed for | How |
|---|---|---|
| the `plumetools` package (Python ≥ 3.9, numpy, scipy ≥ 1.6, PyYAML) | everything | `pip install -e .` (extras below) |
| OpenFOAM **v2512** (standard distribution) | meshing, solving, `postProcess` | your platform's OpenFOAM package or build; **not** the MNF `dsmcFoam+` fork |
| the custom libraries | every DSMC run except 3d-inflow (`libplumeDsmcBoundaryModels.so`); the hybrid CFD stage (`libhybridPlumeThermo.so`) | `applications/*/Allwmake`, once per OpenFOAM installation |
| ParaView's `pvpython` (and optionally ffmpeg) | slice imagery (and videos, contact sheets) | a ParaView install; **not pip-installable** |

## 1. The Python package

From the repository root:

```bash
pip install -e ".[test,plots,cases]"    # what ./Allsetup installs by default
```

| Extra | Brings | Needed by |
|---|---|---|
| `test` | pytest | the test suite |
| `plots` | matplotlib, pandas | figures written by post-processing (tables are written without it) |
| `cases` | CaseFoam (and PyFoam) | generating `cai2012-health` and `markelov1999` studies |

If CaseFoam will not install, use `--extras test,plots` with `Allsetup`, or
`pip install -e ".[test,plots]"`; everything except those two generators works.

**On Windows.** Native Python runs every Python-only step and the test suite.
The repository's `.venv/` is a Linux environment for WSL; create a separate one
for Windows. The shell runners need bash (WSL; Git Bash for the pure-Python
ones). Before cloning, `git config --global core.symlinks true` avoids a
permanently modified file in an archived case.

## 2. OpenFOAM and the libraries

Source OpenFOAM v2512 in the shell (for example
`. /usr/lib/openfoam/openfoam2512/etc/bashrc`), then:

```bash
./Allsetup                                   # package + libplumeDsmcBoundaryModels.so + checks
applications/thermophysicalModels/Allwmake   # libhybridPlumeThermo.so (hybrid only; Allsetup does not build it)
./Allsetup --check                           # verify only
```

`Allsetup` options: `--python-only` (skip OpenFOAM), `--check`, `--extras a,b`.
It checks `blockMesh`, `snappyHexMesh`, `checkMesh`, `dsmcInitialise`,
`dsmcFoam` (failures) and `decomposePar`, `reconstructPar`, `mpirun`,
`foamDictionary` (reported, since `--serial` runs need none of them). Both
libraries install into `$FOAM_USER_LIBBIN`; build them in **every** OpenFOAM
installation that will run a case — a library is not portable between builds.
What `plumeFieldInflow` does and why stock `dsmcFoam` needs it:
[plumeFieldInflow](../plume-field-inflow.md),
[Solver compatibility](../solver-compatibility.md).

## 3. ParaView, for imagery

Install ParaView so that `pvpython` is on `PATH`, or set
`PVPYTHON_PATH=/path/to/pvpython`. Verified with ParaView 5.11.2 (Ubuntu's
package, Python 3.12.3) and 6.2.0 (the paraview.org binary, Python 3.12.14).
`pvpython` must be able to import PyYAML and numpy; check with

```bash
pvpython -c "import sys, yaml, numpy; print(sys.version)"
```

and see [Post-processing](../guides/postprocessing.md#imagery-needs-paraviews-own-interpreter)
if it cannot. Do not try to `pip install paraview` into the package's
environment: the renderer is run *by* `pvpython`, not imported. ffmpeg
(`apt install ffmpeg`) is optional; without it frames are written and videos
are skipped.

## 4. The documentation environment (optional)

Separate from all of the above, Python ≥ 3.11:

```bash
python3 -m venv .venv-docs
.venv-docs/bin/python -m pip install -r docs/requirements.txt
.venv-docs/bin/python -m mkdocs serve
```

Details: [Writing the documentation](../development/documentation.md).

## Check the result

```bash
./Allsetup --check                         # what is present
python -m pytest -m "not needs_openfoam"   # the fast suite, no OpenFOAM
python -m pytest -m needs_openfoam         # only with OpenFOAM v2512 and the library built
```

The fast suite needs no OpenFOAM; the regression tests inside it skip unless
the historical 3d-inflow mesh is restored. The last measured run is recorded in
[What the tests establish](../explanation/evidence.md#measured-for-this-documentation-2026-10-28).
Historical provenance (what produced the archived results, and what cannot be
recovered): [Environment and provenance](../environment.md).

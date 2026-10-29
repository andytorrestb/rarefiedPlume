# Troubleshooting

Failures you are likely to meet, the message that identifies each, and what to
do. Messages are quoted from the scripts and loaders.

## Setup

**`plumetools is not importable`** (any runner). Install the package into the
Python the script finds first: `pip install -e ".[test]"` from the repository
root, or `./Allsetup`. In WSL check you are not using a Windows interpreter,
and on Windows that you are not using the repository's Linux `.venv`.

**`libplumeDsmcBoundaryModels.so not found`, or a case fails at `Allrun` with a
missing boundary model.** Build it in *this* OpenFOAM environment:
`applications/dsmcBoundaryModels/Allwmake` (needs `wmake` and
`$FOAM_USER_LIBBIN`, i.e. a sourced OpenFOAM). `./Allsetup --check` reports
which piece is missing. There is deliberately no fallback to a uniform inflow.

**`libhybridPlumeThermo.so is not built in this environment`** (hybrid `Allrun`).
`applications/thermophysicalModels/Allwmake`; `Allsetup` does not build it.

**CaseFoam will not install, or `casefoam is not importable`.** Needed only to
generate `cai2012-health` and `markelov1999` studies (Markelov's `--dry-run`
too). `./Allsetup --extras test,plots` installs everything else.

**`no $WM_PROJECT_VERSION -- the OpenFOAM environment is not sourced`.** Source
OpenFOAM's `etc/bashrc` in the shell (sourcing it from inside another script
under `set -u` fails: OpenFOAM reads unset variables).

**A cloned working tree on Windows is permanently "modified".** The archived
`cases/caseFoamEx/2d-wedge/system/blockMeshDict` is the repository's only git
symlink; without symlink privileges it becomes a text file. Set
`git config --global core.symlinks true` before cloning (finding AR-14; not
repaired because the case is archived).

## 3d-inflow and the legacy cases

**`no constant/polyMesh. It is gitignored for this case`.** Restore the
historical mesh with `git checkout f83233f^ -- cases/3d-inflow/constant/polyMesh`
(the tag named in the message, `pointwise-mesh-3d-inflow`, does not exist), or
generate one with `./Allmesh --yes` — which overwrites the directory.

**`Zero boundary temperature detected`**, or `Allrun` refusing because patches
"have geometric type 'patch'". Stock `FreeStream` injects on every
`patch`-type boundary. For 3d-inflow set `mesh.outer_patch_type: wall` and
re-mesh (the outer boundary then reflects — a physics change); the other
families avoid it by using `plumeFieldInflow` with an explicit patch list. See
[Solver compatibility](../solver-compatibility.md).

**`pytest tests/regression` fails after meshing.** The golden is bound to the
Pointwise mesh; a generated mesh cannot reproduce it. Restore the historical
mesh, or expect the failure.

**ISS and `Fnum_*` cases.** They need the MNF `dsmcFoam+` and have no run
script; see their tutorials. Do not run their `clean.sh`.

## Generation

**`these cases conflict with study.yaml / baseCase and were NOT changed`**
(Kn extension, interface track). An input differs from what the study would
generate. Nothing was written. Re-run with `--regenerate` if the case has no
solver output; `--overwrite` replaces a case **with** results and destroys
them.

**`Existing cases that have results are refused unless --overwrite`**
(cai2012, health, markelov1999). Regenerating rewrites `case.yaml` under
finished output. To change scripts only, copy them in (or, in the Kn
extension, `--refresh`).

**`baseCase/case.yaml no longer inherits cases/cai2012-health's settings`**
(Kn extension). A template value changed without being declared in
`policy.DELIBERATE_OVERRIDES`. Restore it, or declare it with its reason —
that is a code change.

**`smoke_cases.py` refuses the tree.** It will not run where a mesh or solver
output exists; use `validate_cases.py` there.

**`model is '…', expected '…'`.** The file belongs to another family's loader;
the message names it.

## Meshing

**`Allmesh` refuses to re-mesh a case that has results.** Correct: the results
would address cells that no longer exist. Use `./Allmesh --dict-only` to
rewrite dictionaries only, or start a [fresh iteration](restart-and-recovery.md#changing-the-mesh).
`--force` overrides.

**Nozzle face count or area check fails.** The built staircase differs from the
Python prediction, or the area error exceeds `checks.max_nozzle_area_error`
(5 %). Change the core cell size; the checks refuse a disk resolved by fewer
than `checks.min_nozzle_faces` faces.

**A plan refuses to split the exit across refinement levels** (nozzle plume).
snappyHexMesh refines a cell when its centre is inside a region; widen the exit
cone so the first cell centres are inside it. `--describe` prints the cones.

## Running

**`the saved state of this case cannot be continued as it is`** (exit 3 from
`checkpoints plan`). An incomplete time exists on some rank; follow
[When a run is refused](restart-and-recovery.md#when-a-run-is-refused).

**`foamDictionary is not on PATH, so -np cannot override`.** Set
`dsmc.n_subdomains` in `case.yaml` and regenerate the dictionaries, or source
OpenFOAM so `foamDictionary` exists.

**MPI starts fewer ranks than requested, or refuses.** Open MPI counts physical
cores, not hardware threads; ask for no more ranks than cores.

### A different `decomposePar` is first on `PATH`

On the owner's laptop a user-bin `decomposePar` from another OpenFOAM build
(hyStrath) aborts on every case. Put the v2512 binaries first after sourcing
OpenFOAM: `export PATH="$FOAM_APPBIN:$PATH"`. The hybrid `Allrun` does this
itself.

**The hybrid `dsmc` stage is refused.** `results/interface.json` holds an
*unsuitable* verdict; that is a result. To build the DSMC region anyway, as a
machinery smoke test, run `post-cfd` again with `--allow-unsuitable`.

**`hpc guard: REFUSED`** (Kn extension). The case's Kn is outside its suite
policy, or an interface case is infeasible or was built on another mesh.
**`hpc label: REFUSED`** stops a production interface case without passing
overlap evidence for the same coupling.

## Post-processing and imagery

**`import paraview.simple` fails.** The renderer was run with `python`; use
`pvpython`. If `pvpython` itself cannot import `yaml` or `numpy`, set `PV_VENV`
to a virtual environment of the same Python `major.minor` (see
[Post-processing](postprocessing.md#imagery-needs-paraviews-own-interpreter)).

**`GLXBadContext` at the end of a render.** Harmless when the images were
written: `render_slices.py` reports its own status.

**Plots show a black rectangle or nine different-looking cases.** A pinned
colour range far from the data, or ranges not pinned across cases. Measure the
range off a control case and pin it in the spec.

**`tests/unit/test_viz_dissect.py::test_a_missing_cell_image_is_refused_before_ffmpeg_runs`
fails.** On a machine without ffmpeg the "ffmpeg is not on PATH" check fires
before the missing-image check the test expects (observed on native Windows,
2026-10-28). An environment-dependent test, not a rendering failure; see
[Known gaps](../development/known-gaps.md).

## WSL

**`wsl` fails with `Wsl/Service/E_UNEXPECTED`, or long commands die.** On the
owner's machine this followed an out-of-memory kill of the 15 GB WSL VM. Bound
array sizes, run long jobs detached, and restart WSL (`wsl --shutdown`) only
when nothing else is running in it.

**`pvpython` hangs at start with an unkillable `orted`.** Observed on the
owner's WSL as an MPI start-up stall with kernel soft lockups; only
`wsl --shutdown` cleared it.

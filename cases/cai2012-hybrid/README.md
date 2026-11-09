# cases/cai2012-hybrid — one-way CFD → DSMC free plume

Cai & Wang (2012)'s circular argon source (D = 0.2 m, S₀ = 2), computed as a
`rhoCentralFoam` precursor, a fixed continuum-breakdown surface, and standard
`dsmcFoam` + `plumeFieldInflow` beyond it. Original profiles are 3-D Cartesian;
an opt-in planar ellipse prototype is described below. Smoke
runs are intentionally under-resolved; verification and production profiles are
prepared, not run; nothing is validated.

| file | what |
|---|---|
| `hybrid.yaml` | every input, labelled [PAPER] / [REPO] / [BIRD] / [ASSUMPTION] / [NUMERICAL]; profiles overlay `base` |
| `Allrun` | the OpenFOAM stages, resumable: `setup`, `cfd`, `post-cfd`, `prepare-dsmc`, `dsmc`, `post-dsmc`, `crossings`, `all` |

```bash
# OpenFOAM v2512 sourced, BOTH libraries built in this environment:
#   applications/thermophysicalModels/Allwmake
#   applications/dsmcBoundaryModels/Allwmake
C=runs/kn0p01c                     # runs/ is gitignored
./Allrun $C setup --profile smoke-kn0p01-coarse -np 4   # Python only
./Allrun $C cfd -np 4
./Allrun $C post-cfd               # verdict in results/interface.json
./Allrun $C dsmc -np 4             # refuses an unsuitable interface unless --allow-unsuitable
./Allrun $C post-dsmc              # results/summary.json, *.png
```

Every stage resumes after an interruption; nothing is deleted (to start over,
remove the case directory). Batch scripts:
`python -m plumetools.hybridplume hpc <case> --scheduler pbs|slurm|none --np N`.

To run the default 8M machinery-test profile with the single-case PBS wrapper,
submit `qsub runCase.sh`. Select another profile with Torque's `-v` option; the
case directory defaults to `runs/<profile>`:

```bash
qsub -v PROFILE=verification-kn0p01 runCase.sh
```

The wrapper requests one node with 192 processors so it can accommodate the
largest profiles; the selected profile controls the actual rank count.

- Tutorial (profiles, stages, suitability, resume, outputs, HPC):
  [`docs/tutorials/cai2012-hybrid.md`](../../docs/tutorials/cai2012-hybrid.md)
- Method, decisions and measured smoke results: [`docs/hybridplume.md`](../../docs/hybridplume.md)

## Laptop planar ellipse prototype

Use `planar-laptop-kn0p01` for the complete planar CFD to fitted interface to
DSMC workflow. The source is a slit, with the same exit state and D-based Kn as
the circular source; it is not an axisymmetric reproduction. The initial CFD
mesh is 40 x 40 x 1 (1,600 cells), D/20 spacing and D/20 depth, with front/back
`empty` patches. Serial execution is the default; `-np 2` or `-np 4` is optional.

The fit maximizes clipped elliptical area while remaining inside the original
source-connected continuum-valid cell footprints, including interior exclusions.
It keeps the existing `Kn_GLL` and `abs(B_exp)` tests and uses a default inward
clearance of 0.25 core-cell widths. Outer-domain separation is always at least
one quarter of the smallest CFD cell, including when excluded-cell clearance
is set to zero. Invalid or disconnected pockets stay outside
the ellipse. The deterministic optimizer supplies a local fit, not a proof of
global optimality; infeasibility is reported instead of relaxing the criteria.

The separate DSMC mesh has curved edges projected onto planar sections of a
native `searchableSphere` ellipsoid. The two sections are identical and the mesh
is uniformly extruded. Real interface faces receive bounded, zero-slope states
from their containing accepted CFD cells. Uncovered slit faces inject the
original exit state directly.

Run these commands **from the repository root**, with OpenFOAM v2512 sourced and
the two existing libraries built in that environment:

```bash
C=cases/cai2012-hybrid/runs/planar-laptop-ellipse
R=cases/cai2012-hybrid/Allrun

bash "$R" "$C" setup --profile planar-laptop-kn0p01
bash "$R" "$C" cfd
bash "$R" "$C" post-cfd
bash "$R" "$C" prepare-dsmc       # build mesh and assess/map faces; no solver

# Inspect before continuing:
#   $C/results/breakdown_ellipse.png
#   $C/results/ellipse.json
#   $C/results/interface.json and interface.vtk
#   $C/cfd/planar.foam and $C/dsmc/planar.foam

bash "$R" "$C" dsmc
bash "$R" "$C" post-dsmc
```

`post-cfd` saves the raw-boundary/ellipse overlay. Final normal-Mach, return-flow,
overlap and closure diagnostics are evaluated by `prepare-dsmc` on the actual
curved mesh. An unsuitable interface stops before DSMC initialization. To
explicitly continue a machinery smoke test after reviewing that result:

```bash
bash "$R" "$C" prepare-dsmc --allow-unsuitable
bash "$R" "$C" dsmc
```

The verdict remains unsuitable and reports retain `smoke_only: true`. If no
feasible ellipse exists, this flag cannot manufacture one. Changing the fit,
source, CFD state or mesh configuration after a DSMC mesh exists requires a fresh
case directory. Re-running the same DSMC stage resumes as before.

To adjust the fit, edit `interface.fit.clearance_cells` in the planar profile and
set up a **new** case. Extents, spacing and depth live beside it in `hybrid.yaml`.
The pre-CFD DSMC resource estimate uses exit density throughout as a conservative
sizing bound; final weight, timestep and estimates use the actual curved mesh and
CFD states, including the in-plane cell widths.

After sampling, `post-dsmc` writes `planar_dsmc.png`, `planar_fields.npz` and
`centreline.csv` on the actual DSMC mesh. The centreline file covers the DSMC
portion; the CFD case remains independently viewable. Fluxes use the stated
extrusion depth. Circular-disk collisionless and azimuthal comparisons are marked
not applicable. Matplotlib is optional for plots; the `.foam` cases remain usable.

This intentionally under-resolved prototype is for assessing the interface and
workflow. Development checks covered numerical fitting, cell mapping and native
mesh generation; neither CFD nor DSMC was run for this prototype.

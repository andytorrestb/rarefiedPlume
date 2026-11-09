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

## 3D fitted ellipsoid on the HPC

`hpc-kn1em3-ellipsoid` uses the full circular argon source at **exit Kn = 0.001**,
D = 0.2 m, T = 300 K and speed ratio 2. The domain is x = 0…8 D and
y,z = −5…5 D (1.6 m downstream, 2 m across). There are no wedge, symmetry or
planar reductions. VHS collisions, Minmod and the 0.25 viscous timestep factor
are retained. This is an initial resolution for manual assessment, without a
convergence or cluster performance claim.

On the cluster, submit from the **repository root**:

```bash
qsub cases/cai2012-hybrid/run3DCase.sh
# Alternatively, from cases/cai2012-hybrid:
qsub run3DCase.sh
```

The dedicated wrapper requests one node, **192 cores**, `mem=1200gb`, `longq`
and a 120-hour ceiling. It sources the existing cluster environment and builds
the two plume libraries if absent. Its default stage is `all`: CFD precursor,
breakdown/fit, parallel DSMC meshing, face transfer and suitability checks,
`dsmcInitialise`, parallel `dsmcFoam` startup **and sampling**, then post-processing.
Nothing is submitted remotely by the development workflow.

The default output directory is `cases/cai2012-hybrid/runs/hpc-kn1em3-ellipsoid`.
Use an absolute `CASE_DIR` for a different output location. Resubmit the same
command after interruption to resume complete checkpoints. To resume only DSMC,
or regenerate reports from a completed run:

```bash
qsub -v STAGE=dsmc cases/cai2012-hybrid/run3DCase.sh
qsub -v STAGE=post-dsmc cases/cai2012-hybrid/run3DCase.sh
# Optional different output directory (same value on every restart):
qsub -v CASE_DIR=/absolute/scratch/my-plume cases/cai2012-hybrid/run3DCase.sh
```

The initial Cartesian CFD lattice is 160 × 240 × 240 (9,216,000 cells), with
D/40 core spacing and graded outer cells. The fit maximizes clipped spheroid
volume `(xc,a,b,b)` within the source-connected accepted **3D voxel footprints**,
including off-axis exclusions and interior pockets. Clearance is 0.25 CFD core
cells; an independent all-voxel verification follows optimization. No feasible
fit means a stop, with diagnostics, rather than a relaxed breakdown threshold.
The result is a local constrained optimum; no global optimum is claimed.

Native `snappyHexMesh` uses a `searchableSphere` with vector radius `(a b b)`.
The DSMC background spacing is 0.2 D; interface refinement level is 3. Hollow
refinement regions around the source rim target 0.8 exit mean free paths so the
uncovered source remains resolved without filling the entire box at that size.
Initial caps are 12 million global / 1 million local cells. These controls live
in `dsmc.mesh`; they do not establish a converged collision-cell resolution.
Meshing uses `mesh-build/processor*`, separate from DSMC restart ownership.

Actual snapped facets must stay within a certified continuum-valid enclosure
and within 2% of the fitted radius. Source area coverage must agree with the
intended disk to 1%; density is never rescaled to hide geometry error. The
existing Mach, return-flow, overlap and flux criteria still apply. A larger
domain may help, but does not guarantee a suitable coupling surface.

Resource reports use **GiB = 2³⁰ bytes**, including the historically named
`resources.max_memory_GB`. This profile stops at an estimated **1100 GiB**,
below 1200 decimal GB as well as 1200 GiB, to accommodate site-specific PBS
`gb` units and reserve headroom on the 1.5 TB node. The estimate includes both
meshes, analysis temporaries (25,000 bytes/cell), DSMC parcels/cells, 256 MiB
per MPI rank and a 1.2 safety factor. Controls are under `resources`; the
actual-cell estimate is written before DSMC initialization. It is an estimate,
not measured peak RAM. Particle weight and physical density are not adjusted
merely to pass the memory guard.

Inspect `results/ellipsoid.json`, `interface.json`, `dsmc_fields.json` and
`memory_estimate.json` for the fit, actual-facet certification, source coverage,
suitability and resolution/resource statistics. `breakdown_ellipsoid.png` shows
the raw continuum boundary and smooth spheroid at z = 0. `spatial_dsmc.png`,
`spatial_fields.npz` and `centreline.csv` use the actual independent DSMC mesh;
the centreline covers its downstream portion using local DSMC cell widths for
the finite axis tube and the same cell centres for reference averaging. Plots select thin meridional
bands rather than projecting all azimuths. Open `cfd/spatial.foam`,
`dsmc/spatial.foam` and `results/interface.vtk` in ParaView. Circular-source
reference comparisons are labelled **collisionless-limit comparisons**, not
accuracy errors for this collisional case. Empty cells stay flagged.

An unsuitable verdict stops before DSMC initialization. After examining it,
an explicitly authorized machinery smoke run can use
`qsub -v ALLOW_UNSUITABLE=1 cases/cai2012-hybrid/run3DCase.sh`; the verdict and
`smoke_only` label remain recorded. Infeasible geometry, stale CFD inputs, invalid
mesh geometry and excess memory cannot be overridden. Changes to the fit,
configuration or frozen CFD state require a fresh case directory.

Development ran only focused checks and a synthetic, very coarse native 3D
fixture: 7,134 DSMC cells, two-rank meshing, initialization and five serial DSMC
steps. The production mesh and HPC job were not run on the laptop.

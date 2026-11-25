# Plate-impingement Knudsen sweep on uniGasFoam

Status: approved by the owner, 2026-11-24.

## Intent and why

`cases/cai2012-unigas` runs Cai & Wang's circular argon plume (D = 0.2 m, S0 = 2,
T0 = 300 K) at seven exit Knudsen numbers on uniGasFoam's hybrid USP-DSMC solver,
into vacuum. This family fires **the same thruster, at the same seven Kn_e**,
normal to a flat plate the size of an ISS solar array wing, at three stand-off
distances: z/D = 25, 50 and 100 (5, 10 and 20 m). 7 × 3 = 21 cases.

What the owner wants out of it:

* the flow field, visualised;
* the loads on the plate, **mainly pressure and shear**, sampled face by face;
* those plate loads compared against **Cai's collisionless solution** for a jet
  striking a diffuse plate.

**Timing.** The owner wants to submit the full 21-case sweep to the cluster
about two hours after this design is approved. So the work is split:

* **Before submission:** the generator, the 21 cases, the cluster job scripts,
  and an under-resolved laptop smoke run that proves mesh → run → plate fields
  written → restart.
* **After submission, while the cluster runs:** visualisation, plate-load
  sampling and the Cai comparison, prototyped on under-resolved laptop runs.
  None of them runs inside the cluster job; they read its output.

## Decisions made by the owner

1. 21 cases: Kn_e ∈ {100, 10, 1, 0.1, 0.01, 0.001, 0.0001} × z/D ∈ {25, 50, 100}.
2. Thruster, gas, exit state and uniGasFoam solver settings exactly as deployed
   in `cases/cai2012-unigas` (its `baseCase/unigas.yaml`, fork branch
   `dsmc-time-step-limit`, commit pin `bd8ba5c` or later).
3. Thruster normal to the plate.
4. **No cell budget.** "I expect to have huge runs."
5. **Cells sized by the local mean free path** (rule in §3).
6. Plate loads: pressure and shear first; Cai comparison in scope.
7. Under-resolved laptop runs prototype the workflow; the owner submits the full
   sweep on the cluster.
8. Plate 34.14 × **11.88** × 0.25 m (r-4d's panel).
9. **All 21 cases are submitted without a dense-Kn smoke run.** The owner is
   confident from previous tests that uniGasFoam's USP/DSMC switching works.
   What must be made sure of is that **the comparison workflow is set up**:
   the cluster output carries everything the plate-load sampling and the Cai
   comparison read, and that chain runs end to end on the laptop runs.
10. **Disk is not a constraint** on the cluster ("we rich"): nothing is pruned.

## Calls made by me (the owner can overturn any)

* Thruster aimed at the plate centre.
* Plate fully diffuse at 300 K (`uniGasDiffuseWallPatch`): r-4d's wall and the
  exit temperature.
* **Quarter domain** (symmetry planes y = 0 and z = 0). The plate is centred,
  so a quarter of the cells and particles gives the same answer; totals are
  reported for the whole plate. Fallback: full domain, if the smoke run shows
  uniGasFoam mishandles the symmetry planes.
* A **10 mm floor** on the cell size (§3), so the two densest cases stay
  buildable; cells there larger than λ are what USP is for.
* **One particle weight per case: 40 parcels in the 10 mm exit cell**, at every
  Kn (the first family's dense-Kn rule; §5).
* The Cai comparison covers plate loads only. Once gas reflects off the plate,
  Cai's free-plume solution no longer describes the flow field.
* Names: family `cases/plate-unigas/`, package `plumetools/plateugf/`, cases
  `Kn100_zD25` … `Kn0p0001_zD100`.

## Goals

* G1. `cases/plate-unigas/` generates the 21 cases from one `study.yaml`,
  idempotently, with a `--dry-run` that prints each case's cells, particles,
  memory and steps.
* G2. Each case meshes (blockMesh → snappyHexMesh → nozzle patch) and runs on
  the patched uniGasFoam, writing time-averaged plate loads.
* G3. Per-case PBS jobs, resumable, as in `cases/cai2012-unigas`.
* G4. A laptop smoke run (Kn_e = 100, z/D = 25, coarsened) passes before the
  owner submits.
* G5. After submission: visualisation, plate-load sampling (maps, profiles
  along both plate centrelines, totals), and the Cai plate-load comparison.

## Non-goals

* Comparing the flow field against Cai.
* Any physics change from the first family: gas model, exit state, solver
  settings.
* Making the dense cases (Kn_e ≤ 0.001) resolved. Their limits carry over from
  the first family and are labelled, not fixed.
* Off-centre or inclined plates.

## Anti-goals

* Editing `plumetools.nozzleplume`, `plumetools.cai2012ugf` or the frozen
  `cai2012`/`cai2012ext` packages, or any existing case family. They are
  imported, never modified; the r-4d cases and the first family must generate
  byte-identical files before and after.
* A cell budget sneaking back in: the generator reports size, it does not
  refuse it. A guard against a mistyped level (orders of magnitude) is fine.
* A comparison that hides a difference behind a tolerance. Differences are
  reported, as elsewhere in the repository.
* Any claim of physical validation.

## Constraints

* OpenFOAM v2512; the patched uniGasFoam fork (`dsmc-time-step-limit`,
  ≥ `bd8ba5c`, which carries the empty-rank fix every near-vacuum start needs).
* Cluster: PBS, one node of 192 cores, `longq`, 120 h walltime, `hpc_env.sh`,
  the owner's `#PBS` headers. The cluster build is the first family's
  `runBuild.sh`. Node memory is not recorded in the repository.
* Laptop: WSL, 12 cores, 26 GB to WSL, ~107 GB free on C:. Laptop runs are
  preliminary only (the owner's standing rule).
* Monatomic argon, cell weighting off (first family decisions 4, 5).

## Approaches

1. **New family on r-4d's plate machinery plus the first family's uniGasFoam
   writers (selected).** `nozzleplume` already has the rectangular plate, the
   snappyHexMesh octree that carves it, symmetry planes, a nozzle patch cut from
   the exit plane, and plate-load post-processing; it accepts a Cai exit state
   (Kn, S). `cai2012ugf` already writes every uniGasFoam dictionary and
   converts uniGasFoam output to dsmcFoam's names. New code: mean-free-path
   sizing turned into `nozzleplume` refinement regions, a plate entry in
   `boundariesDict`, the generator, the job scripts, then post and Cai.

   **Builder check correction:** `cai2012ugf`'s case-level functions
   (`render_boundaries_dict`, the other dictionary writers, `convert_case`,
   `hpc.resume_plan`/`finished`, the family `Allrun`) read a cai2012-schema
   `case.yaml` (`geometry.plume_axis`, `exit_state.T0_K`, `cai2012ext.direct_case`).
   `nozzleplume`'s config and `ExitState` differ (`temperature_K`,
   `speed_m_per_s`, no `plume_axis`). Since neither package may be edited,
   `plateugf` either passes the writers a small adapter object exposing the
   attributes they read, or writes its own thin wrappers. Schema-free pieces are
   reused directly: `convert_time`, `lost_clouds`, `mark_restart`, `leg_for`,
   `aligned_interval`, `foamio.checkpoints.plan_restart`,
   `cai2012ugf.settings` (reads `unigas.yaml`). Adapter vs wrapper is the
   builder's call; a test must show the rendered dictionaries match what
   `cai2012ugf` writes for the same exit state, apart from the plate entry and
   the run window.
2. **The first family's graded blockMesh with a plate hole.** A tensor-product
   lattice cannot follow λ: fine cells run in whole slabs across a 20 m domain,
   and cutting a plate from it is new code. Rejected.
3. **Hand-written from uniGasFoam's `plumeImpingement` tutorial.** Fast to a
   first run, but it rewrites the exit state, mesh and post code and drifts from
   both families. Rejected, as in the first family.

## Design

### 1. Geometry and domain

`nozzleplume` conventions: exit disk r ≤ D/2 at x = 0, plume along +x, plate's
upstream face at x = z (5, 10, 20 m), thickness along x. Domain, quarter model:
x ∈ [0, z + 0.25 + 1.75], y ∈ [0, 17.07 + 2], z ∈ [0, 5.94 + 2] (r-4d's margins
behind and beside the plate), extents rounded up to whole background cells by
`nozzleplume`. Patches: `nozzle` (inflow), the exit-plane remainder and every
outer face deleting, `plate` a wall, two `symmetry` planes.

### 2. Exit state, gas, solver

As in the first family: Kn_e → λ0 = Kn_e D → n0 (VHS, ω 0.74, d 4.17e-10 m,
Tref 273 K), U0 from S0 = 2 at 300 K (706 m/s). `uniGasProperties`,
`hybridDecompositionDict`, `fieldPropertiesDict`, `uniGasInitialisationDict`
from `cai2012ugf.dictionaries` with this family's run window. Before
submission, a unit test must show the exit state equals the first family's
manifest values for the same Kn_e.

### 3. Mesh: cells sized by the mean free path

Target cell size at a point, from Cai's collisionless density (no reflected
gas):

    λ_local = λ0 · n0 / n                (VHS, temperature dependence ignored)
    s = clamp( min(λ_local, r/20),  floor = 10 mm,  ceiling = background cell )

r is the distance from the exit centre. The r/20 term resolves the plume's
expansion where λ is long (it is what sizes every cell at Kn_e ≥ 0.1); the
floor keeps Kn_e ≤ 0.001 buildable; USP handles cells larger than λ.

Background cell **0.64 m**, so octree levels give 0.32, 0.16, 0.08, 0.04, 0.02
and **0.01 m** (level 6) exactly. For each level the generator finds the region
where s ≤ that level's size and writes the smallest coaxial cone(s)
`nozzleplume` accepts that contain it (`refinement_regions`, `shape: cone`,
possibly hollow), checked by sampling: every sample point inside the true
region must be inside the cone. Plate surface level: 3 (0.08 m, three cells
through the plate). The exit disk sits wholly at level 6 (20 cells across D).

Idealised counts for this rule (full domain, before octree buffers and the
plate's level-3 skin; quarter domain is ¼ of these; scratch estimate,
2026-11-24, made with a 0.25 m ceiling, so the 0.64 m background gives fewer
far-field cells and the plate skin adds roughly 0.5M):

| Kn_e | z/D 25 | 50 | 100 |
|---|---|---|---|
| 100 … 0.1 | 0.46M | 0.63M | 0.95M |
| 0.01 | 0.47M | 0.64M | 0.96M |
| 0.001 | 1.3M | 1.5M | 1.8M |
| 0.0001 | 24M | 28M | 29M |

`--dry-run` prints each case's estimated cells and memory (~4.3 KB per cell,
measured on the laptop runs). `mesh.max_cells` is set high (1e9) as a typo
guard only.

Meshing: serial snappyHexMesh as r-4d does for every case up to a few million
cells; the Kn_e = 0.0001 cases mesh in parallel (`decomposePar` the background,
`mpirun snappyHexMesh -parallel`, `reconstructParMesh`) in their mesh job. The
builder confirms the parallel path on one case.

### 4. Boundaries (`system/boundariesDict`)

`cai2012ugf.render_boundaries_dict` plus one entry: the plate patch with
`uniGasDiffuseWallPatch` (velocity 0, temperature 300 K). Symmetry patches have
no uniGas model: OpenFOAM's particle tracking reflects specularly at `symmetry`
patches (`uniGasParcel::hitPatch` returns false, so the base class handles
them). Everything else open is `uniGasDeletionPatch`. A diffuse wall sets
`measurePropertiesAtWall`, so `uniGasVolFields` writes, on the plate,
`p_Ar` (pressure), `surfaceShearStress_Ar` (shear magnitude), `fD_Ar` (force
per area, vector) and `surfaceHeatTransfer_Ar` (heat flux), averaged over the
same window as the volume fields; the window set writes the same names with
`_ArWindow`. fD and q have dsmcFoam's sign: momentum and energy the gas hands
the wall, per area per time (builder check, `uniGasVolFields.C`).

### 5. Particle weight

    nEquivalentParticles = n0 · V_exit / 40,   V_exit = (0.01 m)^3

Particle count is then independent of Kn_e: about 6, 12 and 21 million plume
particles in the full domain at z/D = 25, 50, 100 (collisionless estimate; gas
reflected off the plate may double it), a quarter of that in the quarter model. The
first family's 2000-per-20-mm rule for Kn ≥ 1 is not used: on 10 mm exit cells
it would mean 0.3–1 billion particles at Kn_e = 100.

### 6. Run window and writes

`nozzleplume`'s transit rule: one transit = x_max / U0 (≈ 10, 17, 31 ms).
Startup 6 transits, averaging 3 (r-4d's values), so ~0.09, 0.16 and 0.28 s of
simulated time. Writes: once per transit during startup, four per transit
while averaging (12 averaged frames for the stationarity check). Two run legs
as in the first family; uniGasFoam keeps only the newest particle cloud, and
the first family's resume guard (`Allrun`) is reused.

### 7. Cluster jobs

Per case: `runMesh.sh` (serial, or parallel for Kn_e = 0.0001), `runCase.sh`
(192 ranks, 120 h, resumable; resubmit until done), `runPost.sh` (after G5
exists). `runBuild.sh` and `hpc_env.sh` are the first family's.
Decomposition: `scotch` (the octree makes `simple` slabs badly balanced); the
builder may switch to `simple` if the smoke run says so.

**Disk** is not a constraint (decision 10): a dense case may reach a few
hundred GB (builder-check estimate) and nothing is pruned on the cluster.
`--dry-run` still prints a disk estimate per case. The laptop prototype runs
keep the laptop's existing pruning (C: is small).

Rough cost on 192 ranks, from the laptop's 4.7e-7 s per particle-step on 12
ranks, ideal scaling: ~1 s per step at 30M particles; ~25–70k steps at
Δt ≈ 4e-6 s; ~1 day per case at Kn_e ≥ 0.01. The patched DSMC time-step
limit may make Kn_e ≤ 0.001 several times longer (first family R2), so those
jobs are expected to need resubmission.

### 8. The laptop smoke run (before submission)

Kn100_zD25 staged with a coarser background (1.28 m) and levels capped so the
mesh is ~50k cells, weight raised so ~200k particles, 12 ranks, ~500 steps.
Pass when: the mesh builds with plate, nozzle and symmetry patches; uniGasFoam
runs past step 20 (the empty-rank adaptation); particles reflect off the plate
and the symmetry planes (no particle leaves through y = 0 or z = 0); `p_Ar` and
`surfaceShearStress_Ar` on the plate are written, finite and positive where
the plume hits; a restart from the newest write continues; and (decision 9)
**every field the comparison workflow reads is present** in a reconstructed
write: the plate-patch `fD_Ar`/`p_Ar`/`surfaceShearStress_Ar`/
`surfaceHeatTransfer_Ar`, the averaged volume fields, and what
`nozzleplume.post` needs for cell centres and volumes. Takes minutes, not
hours. No dense-Kn smoke run (decision 9).

**The comparison workflow, end to end, on the laptop.** The three Kn_e = 100
cases, under-resolved (coarse background, heavy weight, run one after another)
and **shortened to 2 transits of startup + 1 of averaging** (the owner approved
this: they prove the workflow, not results; ~6, 11 and 19 min on 12 ranks,
rough), go through conversion → plate loads → visualisation → Cai comparison.
The cluster cases keep 6 + 3. This is where §9 is built and proven. It starts as soon as the
smoke run passes and continues after the owner submits; the cluster jobs do not
wait for it.

### 9. After submission

* **Conversion.** `cai2012ugf.convert.convert_time` (the per-time function;
  `convert_case` is tied to the cai2012 schema), plus `fD_Ar → fDMean`,
  `surfaceHeatTransfer_Ar → qMean` on the plate patch, so
  `nozzleplume.post.plate_loads` reads the
  loads unchanged (it computes p = fD·n̂, τ = fD − p n̂, totals, and scales
  quarter-domain totals to the whole plate).
* **Plate loads.** Maps of p and |τ| (and τ_y, τ_z) on the plate face,
  profiles along both plate centrelines, total force and heat; CSV + YAML +
  figures under each case's `results/`, plus a family table over the 21 cases.
* **Visualisation.** `plumetools.viz` slices on y = 0 and z = 0 (n/n0, |U|, T,
  the USP/DSMC map) and the plate surface maps, with fixed colour ranges across
  the sweep.
* **Cai comparison.** Collisionless plate loads computed by quadrature over the
  exit disk, as `plumetools.cai2012.analytical` does for the free plume: at
  each plate point the incident number, normal-momentum, tangential-momentum
  and energy fluxes of the exit Maxwellian seen through the disk, plus diffuse
  re-emission at 300 K:

      p = Π_nn,inc + m Γ_inc sqrt(π k T_w / (2m))
      τ = Π_nt,inc
      q = E_inc − 2 k T_w Γ_inc

  Checked against a far-field point-source limit, and against Cai's published
  circular-jet plate solution where the paper can be obtained. *My call:* the
  2016 Phys. Fluids paper (DOI 10.1063/1.4948365) is the planar-jet one; the
  circular-jet reference is to be identified when this step starts (candidate:
  Cai & He's circular-jet follow-up, MTU digital commons 13986/14244). If no
  circular-jet paper can be had, the published check is dropped and the report
  says so. Output: simulated vs collisionless maps,
  profiles and totals, differences reported without tolerances. Kn_e = 100
  should agree within noise; it is the workflow's end-to-end check.
* Prototyped on under-resolved laptop runs of the three Kn_e = 100 cases.

## Risks

* R1. **uniGasFoam on a snappyHexMesh mesh** (split-hex cells at level
  changes). The first family ran on pure hexes. The builder check found that
  sub-cells use each cell's bounding box, so they work on split hexes; the
  owner is confident in the USP/DSMC switch from previous tests (decision 9).
  The dense cases' first cluster hours are where any remaining problem shows.
* R1b. **Parallel snappyHexMesh** (Kn_e = 0.0001 only) is untested here. Mesh
  jobs are separate from run jobs, so a failure there costs a short job, not a
  run.
* R2. **Dense-case memory and meshing.** Kn_e = 0.0001: ~6–7M cells per quarter
  plus octree buffers; serial snappy is too slow, hence the parallel path.
* R3. **Node memory unknown.** `--dry-run` prints the estimate; the owner checks
  it against the node before submitting the dense cases.
* R4. **Reflected gas** is not in the sizing field: at dense Kn_e it raises
  density near the plate and the cells there may be coarser than λ. Reported
  after the run by the first family's post-run verdict, not fixed beforehand.
* R5. **Run length.** 6 + 3 transits is r-4d's choice for a plate case, not
  measured here; the stationarity check after the run is what decides.

## Left to the builder

* How to fit cones to each level's region (the containment check is required).
* Exact background/levels/weight of the smoke case.
* The parallel snappy job for Kn_e = 0.0001.
* Test layout, following `tests/unit/test_cai2012ugf_*.py` and
  `test_nozzleplume_*.py`.

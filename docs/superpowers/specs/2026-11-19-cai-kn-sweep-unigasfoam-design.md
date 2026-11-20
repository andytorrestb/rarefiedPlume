# Cai Knudsen sweep on uniGasFoam (hybrid USP-DSMC)

Status: approved by the owner, 2026-11-19.

## Intent and why

`cases/cai2012-kn-extension` defines Cai & Wang's circular argon plume
(D = 0.2 m, S0 = 2, T0 = 300 K) at seven exit Knudsen numbers, 100 … 1e-4. Each
case is compared against Cai's collisionless analytical solution and fed to the
whole-field scaling fit (`Allscale.py`). Two routes to production have failed at
the dense end:

* **Standard `dsmcFoam`.** Kn = 0.01, 0.001 and 0.0001 cannot be resolved in the
  4e6-cell budget: cell/λ0 = 2.7, 26 and 270, Δt/τ0 = 0.15, 1.4 and 14.6
  (`docs/cai2012-kn-extension.md` §5). All five recorded runs
  (`results/cost-model.yaml`, up to 41 h on 192 ranks) completed, but the three
  dense runs are labelled unsuitable for production.
* **One-way CFD → DSMC hybrid** (`cases/cai2012-hybrid`). Every recorded
  interface verdict is unsuitable (`docs/development/known-gaps.md`).

uniGasFoam (`../uniGasFoam`, the owner's fork of NVasileiadis93/uniGasFoam;
Vasileiadis et al., CPC 310, 2025) runs the whole domain as particles:

* Cells whose gradient-length Knudsen number Kn_GLL is at most 0.05 are relaxed
  with the unified stochastic particle (USP) BGK scheme. USP cells may be larger
  than λ, and their time step is limited only by Courant.
* The remaining cells use DSMC.

The dsmcFoam suite failed on exactly those two limits, cell size and time step,
in the dense near-exit gas. In the paper's conical plume case uniGasFoam was 43×
faster than pure DSMC at a throat Kn of 1/560.

**The deliverable is unchanged:** seven reconstructed, time-averaged 3-D fields,
compared against Cai's solution with the existing estimators, and entered in the
scaling inventory. Only the solver changes.

## Decisions made by the owner

1. **Full 3-D Cartesian.** Use the same box and the same per-Kn meshes as the
   dsmcFoam suite. An axisymmetric wedge is a fallback only if the pilot shows
   3-D is unaffordable, and would be a new decision.
2. **All seven Kn run on uniGasFoam.** The sweep comes from one solver, and the
   completed dsmcFoam runs at Kn = 0.1 and 1 become an independent cross-check.
3. **Output conversion (Approach 1 below).** uniGasFoam's output is rewritten
   as dsmcFoam's field names, so the existing post-processing runs unchanged.
4. **Patch uniGasFoam (the owner's fork) so DSMC cells respect the
   collision-time limit, and make the sub-cell cap configurable** (§1).
5. **Cell weighting off.** Use the dsmcFoam suite's single particle weight and
   near-vacuum start (§4). This keeps the Kn = 0.1 and 1 cross-check down to the
   solver alone, keeps particle counts near the old ones, and keeps the
   weight-consistency check meaningful.
6. Proposed and not objected to:
   * a new case family and package beside the frozen ones;
   * the paper's solver settings where they are not overridden by decisions 4
     and 5;
   * a suitability verdict made after the run from per-cell fields;
   * a pilot before production.

## Goals

* G1. A patched uniGasFoam builds against OpenFOAM v2512, both locally (WSL) and
  on the cluster.
* G2. A new case family generates the seven uniGasFoam cases from one
  `study.yaml`, reusing the existing physics, mesh, weight and run-window code.
* G3. A converter makes the existing post-processing run on a uniGasFoam case:
  `Allpost`, `Allcompare.py`, `Alldiagnose.py`, `Allaudit.py`,
  `Allconverge.py` and the viz specs.
* G4. Each case gets a post-run suitability verdict computed from its DSMC cells.
* G5. A pilot at a dense Kn measures the following before production is
  committed:
  * where the USP/DSMC split lands;
  * the time step and collision spacing in DSMC cells;
  * the particle count;
  * the cost per step;
  * the disk per write.
* G6. Seven resumable PBS jobs, one per case, as in the dsmcFoam suite.
* G7. A cross-check report compares uniGasFoam with dsmcFoam at Kn = 0.1 and 1
  with identical estimators on identical meshes.
* G8. uniGasFoam cases enter the scaling inventory as a separate solver stratum,
  with suitability decided by G4.

## Non-goals

* Axisymmetric or planar geometry.
* Fixing the direct-exit inflow-surface limitation at Kn ≤ 0.001. It stays
  deferred and is carried over as a label.
* The interface track and the CFD → DSMC hybrid.
* Pure SP or pure USP runs, the SP-DSMC hybrid, and cell weighting.
* Any change to the physics: gas model, exit state, domain, startup reference
  or averaging window.
* Contributing changes upstream to uniGasFoam.

## Anti-goals (failure even if it runs)

* Editing the frozen packages (`plumetools.cai2012`, `plumetools.cai2012ext`)
  or any existing case family, manifest or results tree.
* Relaxing a criterion so that a case looks suitable, or giving a dense case a
  production label based on pre-run checks.
* A suitability check that passes on the nominal sub-cell size while collision
  partners are actually drawn from the whole cell (§7).
* A converter that changes any value. It must only re-express the same data.
* Pooling uniGasFoam and dsmcFoam fields without labels. Every result records
  its solver, and the scaling fits keep the two apart.
* A uniGasFoam patch that changes the published code's behaviour when its new
  settings are not used.
* Any claim of physical validation. The repository's standing rule applies.

## Constraints

* **OpenFOAM v2512.**
  * On the cluster it is `/home2/torresfa/scratch/OpenFOAM/OpenFOAM-v2512`,
    built from source and sourced by `hpc_env.sh`.
  * Locally it is the WSL install already used to mesh the earlier smoke case.
  * uniGasFoam's history shows fixes up to v2412 (`3cd2060`) and its
    `INSTALL.md` says v2212 or later. **v2512 is untested.**
* **Monatomic only.** uniGasFoam's USP and hybrid modes support a single
  monatomic, non-reacting species. Argon qualifies.
* **Cluster:** PBS, one node of 192 cores, `longq` with a 120 h walltime, and
  `hpc_env.sh`. The `#PBS` headers are the owner's own settings.
* **uniGasFoam changes** go on a branch of the owner's fork
  (`github.com/andytorrestb/uniGasFoam`), and every case records the commit it
  ran.
* **Disk.** The dsmcFoam suite was estimated at about 212 GB per case at the
  health write cadence.

## Approaches

1. **Convert outputs to dsmcFoam's field names (selected).**
   * New code writes the uniGasFoam input dictionaries; everything physical and
     geometric comes from `plumetools.cai2012` unchanged.
   * After `reconstructPar`, a converter writes `rhoNMean`, `momentumMean`,
     etc. from uniGasFoam's fields, and the existing analysis then runs on it.
   * Analysis is identical across solvers, so a difference in the Kn = 0.1 and
     Kn = 1 cross-check is attributable to the solver.
2. **Standalone family templated from uniGasFoam's `plumeImpingement` tutorial,
   with new analysis reading uniGasFoam names.**
   * Fastest to a first run.
   * It duplicates the exit-state, mesh and analytical-comparison code, and its
     estimators would drift from the dsmcFoam suite's. That undermines the
     cross-check. Rejected.
3. **Teach `plumetools` uniGasFoam's field names everywhere.** About 20 modules
   read dsmcFoam names directly, several of them frozen. Rejected under the
   frozen-package rule.

## Design

### 1. Building and patching uniGasFoam

**What the published code does.** Two behaviours matter
(`src/lagrangian/uniGas/dynamicAdaptation/uniGasDynamicAdapter.C`).

*The time step.* `calculateTimeStep` (around line 347) applies
`maxTimeStepMCTRatio` only when `bgkCollisionModelName()` is not
`unifiedStochasticParticleSBGK`. With USP selected, **no cell is
collision-time-limited, DSMC cells included**, and Δt is set by Courant alone.

* That contradicts the paper's eqs. 12–14, which say the hybrid Δt is the minimum
  of a collision-time-limited DSMC step and a Courant-limited USP step.
* Left as published, DSMC cells near the nozzle lip would run at Δt/τ of about
  0.8 at Kn = 0.01 and about 7 at Kn = 0.001 (builder-check estimate). That is
  the dsmcFoam failure again.

*The sub-cell cap.* `maxSubCellLevels_` is hard-coded to 10 per direction
(line 49).

**The patch.** It goes on a branch of the fork and is opt-in, so behaviour is
unchanged unless enabled:

* A new `adaptiveProperties` key, e.g. `dsmcCellsTimeStepLimited`, default
  `false`. When `true`, the MCT-ratio limit applies to every cell whose collision
  model is DSMC (`cellCollModelId()[cell] != bgkCollModelId()`), whatever the BGK
  model. These cases set it `true`.
* A new key for the sub-cell cap, e.g. `maxSubCellLevels`, default 10.
* **Expectation:** with cell weighting off, a sub-cell holding one particle
  falls back to a partner from the whole cell (`noTimeCounter.C`, around
  line 224). A cell with N particles therefore gains nothing from more than
  about N/2 sub-cells. At 40 particles per exit cell the cap of 10 per direction
  is not what limits these runs. The key is added because it is two lines and
  harmless; the pilot reports whether any cell is ever at the cap.
* **Test:** on a small case with USP selected, Δt is
  min(Courant, 0.2 τ_min over DSMC cells) when the key is on, and matches the
  published behaviour when it is off.

**Build sequence.**

1. Locally, in WSL with v2512 sourced, run `./Allwmake` in `../uniGasFoam`.
2. Run the `plumeImpingement` tutorial for a few hundred steps as a smoke test.
   This shows the solver runs, not that its results are right.
3. On the cluster, the family provides a build step in the style of the hybrid
   wrapper (build if absent, into `$FOAM_USER_LIBBIN` / `$FOAM_USER_APPBIN`).
   It checks that the built commit is the one the cases expect.
4. If v2512 fails to compile, fix it on the same fork branch and record what
   changed.

### 2. Case family layout

`cases/cai2012-unigas/` mirrors `cases/cai2012-kn-extension/`:

* `study.yaml`, `baseCase/case.yaml`, `generate_cases.py`, `validate_cases.py`;
* per-case `runMesh.sh`, `runCase.sh`, `runPost.sh`, `runViz.sh`;
* `AllpostCases`, `AllreconstructCases`, and a `results/` tree with this
  family's own cost model.

New code lives in a new package beside `cai2012ext`, e.g.
`plumetools/cai2012ugf/`. The seven case names keep the suite's Kn labels
(`Kn100` … `Kn0p0001`), so the cross-check can pair them.

### 3. Reused unchanged

* Kn → λ0 → n0 and the exit state (`plumetools.cai2012.gas`, `inflow`).
* The particle weight: `nEquivalentParticles = n0 V_exit / target`, with the
  kn-extension's per-Kn targets (2000 at Kn ≥ 1, 40 at Kn ≤ 0.1). It must equal
  the dsmcFoam suite's manifest value for the same Kn.
* The mesh:
  * the kn-extension `study.yaml`'s per-Kn mesh overrides: Kn ≥ 1 on the
    uniform 20 mm lattice (4.0M cells), Kn ≤ 0.1 from the graded planner
    within 4e6 cells;
  * `blockMeshDict`, `topoSetDict`, `createPatchDict`, `meshQualityDict` and
    `decomposeParDict` (192 ranks, `simple (4 6 8)`) from
    `plumetools.cai2012.dictionaries`.

  The generated mesh dictionaries must be byte-identical to the dsmcFoam
  suite's for the same Kn. `validate_cases.py` checks this.
* The initial state: 1e10 m⁻³, 300 K, at rest
  (`dsmc.initial_number_density_per_m3`).
* The run window:
  * `average_start_s` = 0.0501435 s (Cai's 10⁴ t0 at the Kn = 0.01 exit state,
    17.7 transits);
  * 4.5 transits of averaging (0.01273 s);
  * the same end time.
* The post-processing scripts are copied from the kn-extension `baseCase`
  byte-identical, except as stated in §7.

### 4. New: the uniGasFoam inputs

| File | Content |
|---|---|
| `constant/uniGasProperties` | `collisionModel hybrid`; `bgkCollisionModel unifiedStochasticParticleSBGK`; `dsmcCollisionPartnerModel noTimeCounter`; `dsmcCollisionModel variableHardSphere`; `collisionProperties { Tref 273; macroInterpolation true; theta 1e-3; }`. Argon from the repository's gas model: mass 6.63e-26 kg, d 4.17e-10 m, **ω 0.74** (the tutorial's 0.81 is not used), 0 rotational and vibrational modes. `nEquivalentParticles` as in §3. **`cellWeightedSimulation false`**. `axisymmetricSimulation false`. `adaptiveSimulation true` with `cellWeightAdaptation false`, `subCellAdaptation true`, `timeStepAdaptation true`, `adaptationInterval 20`, `maxSubCellSizeMFPRatio 0.5`, `maxTimeStepMCTRatio 0.2`, `maxCourantNumber 0.5`, and the two patch keys from §1 (DSMC time-step limit on; sub-cell cap as the pilot sets it, default 10). |
| `system/hybridDecompositionDict` | `localKnudsen`, `breakdownMax 0.05`, `theta 0.2`, `smoothingPasses 5`, `decompositionInterval 100`, `resetAtDecomposition on`, `resetAtDecompositionUntilTime` = `average_start_s`. The values are the paper's text, where the tutorial differs (50 passes, interval 1000). *My call.* Every cell is USP until the first decomposition at step 100. That is harmless from a near-vacuum start, but the pilot confirms it. |
| `system/boundariesDict` | `nozzle`: `uniGasDeletionPatch` plus `uniGasFreeStreamInflowPatch` (n0, T0 = 300 K, U = (U0, 0, 0), Ar). `upstreamVacuum` and `vacuum`: `uniGasDeletionPatch`. Patch names come from `mesh.patch_names`. |
| `system/fieldPropertiesDict` | Two `uniGasVolFields` entries, both with `measureMeanFreePath true` (it is required for `subCellSizeMFPRatio_*` and `timeStepMCTRatio_*`) and `averagingAcrossManyRuns true` (resumable). (a) Field `Ar`: `resetAtOutput on` with `resetAtOutputUntilTime` = `average_start_s`, giving a running mean over the sampling window like dsmcFoam's `fieldAverage`. (b) A second field, e.g. `ArWindow`: reset at every write through the end time, giving per-write-interval averages that stand in for dsmcFoam's single-step `rhoN` (§5). |
| `system/uniGasInitialisationDict` | `uniGasMeshFill` at the §3 initial state. |
| `system/controlDict` | `application uniGasFoam`, `startFrom latestTime`, `writeControl adjustableRunTime`, two-cadence writes (§6). |
| `system/fvSchemes` | Whatever interpolation schemes uniGasFoam's cloud requires. The dsmcFoam suite's are empty, and the builder confirms against v2512 (commit `de0ed27` removed some of these requirements). |
| `0/` | Whatever uniGasFoam requires for the above. The builder confirms this; the tutorial ships `numberDensity_Ar`, `transT`, `rotT`, `vibT`, `elecT`, `U` and `uspRhoNMean_Ar`. |

### 5. New: the converter

The converter runs after `reconstructPar`, on every time directory that holds
the `Ar` and `ArWindow` fields. It writes:

| dsmcFoam field read by plumetools | Written from uniGasFoam |
|---|---|
| `rhoNMean` | `rhoN_Ar` |
| `rhoMMean` | `rhoM_Ar` |
| `momentumMean` | `rhoM_Ar · UMean_Ar` |
| `linearKEMean` | `½ rhoM_Ar |UMean_Ar|² + (3/2) rhoN_Ar k translationalT_Ar` |
| `dsmcRhoNMean` | `uniGasRhoNMean_Ar`, the time-averaged computational-particle count per cell, the same meaning as dsmcFoam's (`audit.py`). uniGasFoam writes **0.001** for empty cells; the converter writes 0 wherever `rhoN_Ar` is 0. |
| `rhoN` (instantaneous) | `rhoN_ArWindow` |

**Exactness.** The builder check confirmed from the source that `UMean` is
mass-weighted, and that `translationalT` comes from the same Δt-weighted
accumulators and window as `rhoN`. `plumetools.cai2012.post.read_case` recovers
U = momentum/ρ and T from ½⟨c²⟩, so the round trip is exact.

* A unit test runs synthetic uniGasFoam fields through the converter and then
  `read_case`, and requires agreement to 1e-12 relative.
* Empty cells (`rhoM` = 0) must stay empty.

**Safety.** Converted files carry a header note naming the converter and its
source fields. The converter never overwrites a uniGasFoam field, and it refuses
a directory where the target files already exist and differ.

**Consequences of the conversion.** These must be documented:

* **Weight consistency.** With cell weighting off, rhoN V / dsmcRhoN is the
  global weight, so `Allaudit.py`'s check keeps its meaning and must pass.
* **Convergence series.** `Allconverge.py`'s instantaneous series becomes a
  per-write-interval average. It is less noisy than a single step, so the trend
  test's power differs; document this in the family's README.
* **Diagnostics validity.** The `Alldiagnose.py` validity rule (> 20 particles
  per step) still reads a particle count and keeps its meaning.

### 6. Write cadence and restarts

This generator is not frozen, so the dsmcFoam suite's open problem 2 can be
fixed here:

* **Startup** (to `average_start_s`): write every transit.
* **Sampling:** the health cadence, 3.536e-4 s (1/8 transit).

How to switch is the builder's call: either the job script edits
`writeInterval` at `average_start_s` (`runTimeModifiable yes`), or the run is
two solver invocations.

**Particle clouds.** By default uniGasFoam deletes the particle cloud from every
write except the latest (`cleanLagrangian`; `-keep-lagrangian` disables it).
*My call:* keep the default. Nothing in post-processing reads parcel clouds, and
it removes most of the disk problem.

The cost is that a run killed while writing can leave no complete cloud to
restart from. The job script must therefore:

1. check, before resuming, that the newest time directory has a complete cloud
   on every rank;
2. otherwise fall back to the previous complete write, keeping `purgeWrite` at
   2 or more. If uniGasFoam has already deleted the previous cloud, this falls
   back to a manual restart, which must be documented.

The pilot measures the disk used per write.

### 7. Suitability, after the run

**Pre-run (generation)** keeps only the checks that still apply:

* nozzle faces and nozzle area error;
* mesh budget;
* the particle-count estimate.

cell/λ0, Δt/τ0 and the fixed-Δt Courant number are reported but are not gates.
USP cells are exempt from the first two by design, and Δt is adaptive.

**Post-run**, each case gets a verdict over the sampling window. The DSMC cells
are those with `uniGasCollisionModelId` = 1 (USP = 0). Over those cells:

* **Collision spacing.** Use the *effective* collision spacing over λ, not
  `subCellSizeMFPRatio_Ar` alone. Effective spacing is the sub-cell size where
  sub-cells typically hold two or more particles, and the cell size where they
  don't, because of the whole-cell fallback in §1. The gate is
  `checks.max_cell_over_mfp` (1.0, Cai's Δx/λ0 = 1), with the distribution
  reported. *My call* on the definition; the builder documents the estimator.
* **Time step.** `timeStepMCTRatio_Ar` ≤ 0.2 (enforced by the §1 patch,
  verified here).
* **Occupancy.** Particles per DSMC cell ≥ `checks.min_particles_per_cell` (5).

Over all cells:

* the USP fraction and where it sits;
* the `KnGLL` map;
* the `Allconverge.py` stationarity verdict.

Labels are `production-candidate` or `non-production`, with reasons. Kn ≤ 0.001
also carries the inflow-surface limitation whatever the verdict.

**Where the verdict is shown.** The kn-extension `Allcompare.py` and
`Allconverge.py` call `interface_case.post_context`, which calls
`feasibility.assess_case`. That writes the dsmcFoam pre-run label into
`comparison.yaml`, `stationarity.yaml` and every figure banner.

*My call:* this family's copies differ from the originals in one place only,
which verdict they load. They use a family-level `post_context` that supplies
the post-run verdict. A test diffs each copy against its original and allows
only that change.

### 8. Pilot

1. **Local chain smoke.** Use a coarsened copy of a case, about 58k cells like
   the earlier `Kn10` smoke, run with the patched solver for a few hundred
   steps. Then run generate → mesh → run → reconstruct → convert → `Allpost` /
   `Allcompare.py` / `Alldiagnose.py` / `Allaudit.py`. This proves the plumbing
   only.
2. **Cluster pilot at Kn = 0.001** on its production mesh with 192 ranks, for
   about 2 transits (not to steady state). The production mesh is needed
   because Kn_GLL depends on Δx.

   *Why Kn = 0.001:* near the lip, gradients resolve only at cell scale, so
   Kn_GLL ≈ λ0/Δx:

   | Kn | λ0/Δx at the lip | Expected cell type |
   |---|---|---|
   | 0.0001 | 0.004 | USP |
   | 0.001 | 0.04 | borderline (threshold 0.05) |
   | 0.01 | 0.4 | DSMC |

   With the patch, a lip that goes DSMC at Kn = 0.001 gives Δt ≤ τ0/5 ≈ 1.0e-7
   s, against the dsmcFoam run's 7.1e-7 s: about 7× more steps.
3. **The pilot reports:**
   * the Δt history, and which cells set it;
   * DSMC/USP cell fractions and their locations;
   * the effective collision spacing / λ distribution in DSMC cells (§7), and
     whether any cell reaches the sub-cell cap;
   * the particle count;
   * wall time per step and rank load imbalance;
   * disk per write;
   * a projected wall time for the full 22.2 transits.
4. **Go/no-go is the owner's decision**, made on that report. The levers then
   available:
   * more particles per cell for the dense cases, against the one-weight policy
     (this improves collision spacing in dense DSMC cells);
   * the axisymmetric wedge;
   * accepting non-production labels.

### 9. Production and cross-check

* Seven jobs, one per case, on 192 ranks, resumable (§6). Wall clock is
  recorded in this family's own cost model.
* **Order:** pilot, then Kn = 0.1 and 1 (the cross-check), then the rest.
* **Cross-check report.** Put `Allcompare.py` metrics and centreline, angular
  and plane profiles for dsmcFoam and uniGasFoam side by side at Kn = 0.1 and 1.
  * Same meshes and same particle weight.
  * Differences are reported, with no tolerance, as the project does elsewhere.
  * At Kn ≥ 0.1, Kn_GLL far exceeds 0.05 almost everywhere, so uniGasFoam is in
    effect DSMC throughout.
  * The remaining numerical difference is the time step: adaptive, at Co 0.5 on
    max(u, v0) and 0.2 τ, against a fixed step at Co 0.2 on U0 + 3σ. The report
    names it.

### 10. Scaling inventory

The frozen `cai2012ext/scaling_inventory.py` cannot take these cases:

* it registers only paths inside its own study root;
* it has no solver field;
* it keys strata on physics only, so two solvers would be pooled;
* its `reference_suitability` requires the dsmcFoam pre-run feasibility to pass,
  and that fails by design here.

*My call:* the new package gets its own inventory step.

* It reads this family's manifest.
* It tags each entry with its solver.
* It decides suitability from the §7 post-run verdict.
* It hands the entries to the frozen fitting and reporting functions unchanged,
  with outputs under `results/scaling/` in this family.
* A combined report, if wanted, shows the solvers as separate strata. It never
  fits across them.

The builder confirms that the frozen fit functions accept externally prepared
entries. If they don't, the builder reports back rather than copying them.

### 11. Provenance per case

Each case records:

* solver (uniGasFoam hybrid USP-DSMC);
* the uniGasFoam fork commit and its patch keys' values;
* OpenFOAM version;
* converter version;
* every `uniGasProperties` and `hybridDecompositionDict` value;
* the pre- and post-run verdicts.

### 12. Risks

* **R1. Build.** v2512 compatibility of uniGasFoam is unknown.
* **R2. Time step.** With the patch, DSMC cells at the nozzle lip may make the
  global Δt about 7× smaller at Kn = 0.001 (§8). The pilot measures it.
* **R3. Collision spacing in dense DSMC cells.** With one weight and about 40
  particles per exit cell, collision spacing in dense DSMC cells cannot fall far
  below the cell size. Cases where such cells exist may fail §7. The pilot
  measures it, and the lever (more particles) is the owner's (§8).
* **R4. USP statistics.** USP relaxation uses each cell's sampled stress and
  heat flux, time-averaged with θ = 1e-3, and 40 particles per exit cell is
  thin. Watch the USP cells' noise in the pilot.
* **R5. Load balance.** `simple (4 6 8)` with a USP/DSMC mix.
  `uniGasDynamicLoadBalancing` only reports imbalance and never rebalances. The
  pilot measures it.
* **R6. Restarts after a cut-off write** (§6).
* **R7. Inflow surface.** The direct-exit inflow-surface limitation at
  Kn ≤ 0.001 is unchanged.

## Left to the builder

* Package and directory names, and the exact names of the patch keys.
* How to switch write cadence.
* The exact `0/` files and `fvSchemes` entries uniGasFoam needs.
* The coarsened smoke mesh size.
* The effective-collision-spacing estimator (§7), documented.
* Test layout, following `tests/unit/test_cai2012ext_*.py`.

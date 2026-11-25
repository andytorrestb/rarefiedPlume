# cases/plate-unigas — Cai's thruster firing onto a plate, on uniGasFoam

The first family's thruster (`cases/cai2012-unigas`: argon, D = 0.2 m, S0 = 2,
T0 = 300 K, the same uniGasFoam hybrid USP-DSMC settings) at all seven exit
Knudsen numbers, fired normal to the centre of a 34.14 × 11.88 × 0.25 m plate
(r-4d's panel) at z/D = 25, 50 and 100 (5, 10, 20 m): **21 cases**. One quarter
of the domain is modelled (symmetry planes y = 0, z = 0). Cells are sized by the
local collisionless mean free path, clamped to 10 mm … 0.64 m, with no cell budget.

Design: [`docs/superpowers/specs/2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md`](../../docs/superpowers/specs/2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md).
Plan: [`docs/superpowers/plans/2026-11-24-plate-impingement-kn-sweep-unigasfoam.md`](../../docs/superpowers/plans/2026-11-24-plate-impingement-kn-sweep-unigasfoam.md).

> **Nothing here is validated.** The Cai comparison reports differences; it
> does not judge them.

## Files

| file | what |
|---|---|
| `study.yaml` | the 7 × 3 matrix, the sizing rule, the smoke-run overrides |
| `baseCase/case.yaml` | the template (nozzleplume schema); `baseCase/unigas.yaml` is the first family's, verbatim |
| `generate_cases.py` | writes `Cases/<case>/` and `manifest-plate.yaml` (cells, parcels, memory, steps per case) |
| `Cases/<case>/runMesh.sh`, `runCase.sh` | per-case PBS jobs: mesh (serial snappyHexMesh), run (192 ranks, resumable) |
| `submit_all.sh` | queues every case's mesh job, then its run job after it (`afterok`) |
| `Cases/<case>/Allpost`, `runPost.sh` | post-processing: plate loads, the Cai comparison, slices (below) |
| `Smoke/Kn100_zD25`, `Prototype/Kn100_zD*` | coarse laptop cases: the smoke run and the phase-B prototypes |

## Pre-run estimates (quarter domain, `manifest-plate.yaml`)

| Kn_e | cells | memory | note |
|---|---|---|---|
| 100 … 0.01 | 0.31–0.33M | 2.5–3.7 GB | the r/20 geometric rule sizes every cell |
| 0.001 | 1.1M | 6–7 GB | |
| 0.0001 | 24–34M | 106–147 GB | **check node memory before submitting** |

Parcels ≈ 4.9, 6.5 and 9.7 million at z/D = 25, 50, 100 (40 parcels in the 10 mm
exit cell, one weight per case). Steps are nozzleplume's fixed-Δt count;
uniGasFoam adapts Δt (Courant 0.5) and takes roughly 2.5× fewer.

## On the cluster

```bash
git fetch && git checkout feat/plate-unigas        # or main once merged
cd cases/cai2012-unigas && qsub runBuild.sh         # only if the fork is not built at bd8ba5c
cd ../plate-unigas
bash submit_all.sh                                 # all 21: mesh job, then run job after it
bash submit_all.sh Kn100_zD25 Kn100_zD50           # or a subset
```

A run job that hits the 120 h walltime stops cleanly at its last write:
`cd Cases/<case> && qsub runCase.sh` again resumes it. The two densest Kn_e are
expected to need resubmission.

## Settings and why

* **Solver:** the first family's `unigas.yaml` (hybrid USP-DSMC, DSMC cells
  collision-time limited, cell weighting off), fork commit `bd8ba5c`.
* **Plate:** `uniGasDiffuseWallPatch` at 300 K; uniGasFoam writes `p_Ar`,
  `surfaceShearStress_Ar`, `fD_Ar` and `surfaceHeatTransfer_Ar` on it.
* **Run window:** 6 domain transits of startup, 3 of averaging (r-4d's), one
  transit = x_max / U0. Startup writes once per transit, sampling four per transit.
* **Mesh:** octree levels of a 0.64 m background, regions fitted as cone
  segments to where the target cell is smaller than each level
  (`plumetools/plateugf/sizing.py`); the plate at level 3 (0.08 m).

## Post-processing (phase B)

After a run: `cd Cases/<case> && qsub runPost.sh` (or `./Allpost` by hand). It
reconstructs the newest writes, writes cell centres and volumes, then fills
`results/`:

| output | what |
|---|---|
| `plate_loads.csv`, `plate_summary.yaml` | per-face pressure, shear, heat flux, force; whole-plate totals (r-4d's estimators) |
| `plate_profiles.csv`, `plate-{pressure,shear,heat-flux}.png` | the upstream face: maps and the two centreline profiles |
| `cai_comparison.{csv,yaml}`, `cai-profiles.png`, `cai-ratio-maps.png` | simulation vs Cai's collisionless plate solution (`plumetools/plateugf/caiplate.py`: the exit Maxwellian integrated over the nozzle disk, diffuse 300 K re-emission) |
| `slice-{n,U,T}-{y0,z0}.png` | density/exit, speed, temperature on the two symmetry planes, fixed colour ranges |

`python -m plumetools.plateugf family-table Cases` tabulates the comparison
over every post-processed case. Maps are normalised by the exit static pressure
p0 = n0 k T0, so all 21 cases share one colour scale.

At Kn_e = 100 the simulation should approach the collisionless solution. On
the laptop smoke run (Kn100_zD25, coarse and short) the whole-plate force was
1.008 × Cai's and the centre pressure 1.007 ×; this checks the chain, it
validates nothing.

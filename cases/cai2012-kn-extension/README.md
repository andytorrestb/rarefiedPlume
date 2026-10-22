# cases/cai2012-kn-extension — the Cai 2012 Knudsen extension

Cai & Wang's 2012 circular argon free plume (`D = 0.2 m`, `S0 = 2`,
`T0 = 300 K`, uniform drifting Maxwellian into vacuum, collisions on) at **four
Knudsen numbers Cai did not report**, at **one baseline parcel weight** —
40 parcels in the exit cell, `cases/cai2012-health`'s highest-resolution
setting:

| case | Kn | verdict |
|---|---|---|
| `Kn10_np1x` | 10 | production candidate |
| `Kn1_np1x` | 1 | production candidate |
| `Kn0p001_np1x` | 1e-3 | **unsuitable for production** |
| `Kn0p0001_np1x` | 1e-4 | **unsuitable for production** |

Four configurations and no others; there is no parcel-population sweep. The
`_np1x` suffix is the parent study's naming — its loader names a case by its
particle level — and here it is always the baseline. **Kn = 0.01, 0.1 and 100,
Cai's own conditions, are excluded**: generation, meshing, running and
post-processing all refuse them, even as a disabled entry in `study.yaml` or a
hand-edited row in `manifest.yaml`. They belong to
[`cases/cai2012`](../cai2012/README.md).

Every case carries the **nozzle-plume family's resolution diagnostics**: its
pre-run checks decide the verdicts above, and after a run `Alldiagnose.py`
writes the same mean-free-path, Knudsen-number and cell-resolution fields R-4D's
`./Allpost` does ([`docs/rarefaction-diagnostics.md`](../../docs/rarefaction-diagnostics.md)).

The science — what is compared, how, against what — is in
[`docs/cai2012-kn-extension.md`](../../docs/cai2012-kn-extension.md). This file
is how to run it.

> **Nothing here is a validated result yet.** No case has been run. The two
> dense Knudsen numbers cannot be resolved within the inherited cell budget and
> are marked unsuitable; the two rarefied ones meet every pre-run criterion,
> which makes them candidates, not results. Every comparison is against Cai's
> **analytical** collisionless solution — there is no published DSMC at any
> Knudsen number in this study.

---

## What is here

| | |
|---|---|
| `study.yaml` | the matrix, in `cases/cai2012`'s own study schema |
| `baseCase/case.yaml` | the template: `cases/cai2012-health`'s, with five declared overrides |
| `baseCase/diagnostics.yaml` | the nozzle-plume `diagnostics:` thresholds, per case |
| `generate_cases.py` | idempotent generation, provenance, feasibility → `manifest.yaml` |
| `validate_cases.py` | independent recomputation + the parent family's dictionary check |
| `AllmeshCases` | the parent's, with the Kn policy gate |
| `AllpostCases` | per case: post, diagnose, compare, audit, stationarity; imagery; `analyse.py` |
| `AllcleanCases` | the parent's, byte-identical |
| `analyse.py` | the study table and the four-Kn comparison figure |
| `viz.yaml` | flow-field imagery; per-Kn specs are written from the manifest |
| `viz-diagnostics.yaml` | images of the diagnostic fields, R-4D's list on this plane |
| `hpc_env.sh` | the cluster environment every job sources: modules, OpenFOAM, the venv, ParaView 6.2.0 |
| `runMeshAll.sh` | PBS job: generate, validate, mesh all four, validate again |
| `runAnalyse.sh` | PBS job: the study table and figure, once every case is post-processed |

Per case (`baseCase/`, copied into every case):

| | |
|---|---|
| `Allmesh`, `Allrun`, `Allpost`, `Allclean`, `postProcess.py`, `runInflow.py` | `cases/cai2012/baseCase`'s, **byte-identical** (a test pins it) |
| `Allaudit.py` | `cases/cai2012-health/baseCase`'s, **byte-identical**: occupancy and per-frame centreline error |
| `Alldiagnose.py` | the nozzle-plume rarefaction and mesh-resolution diagnostics |
| `Allcompare.py` | r/D = 10 angular density and mass flux, off-axis velocity, temperature, pressure |
| `Allconverge.py` | far-field stationarity across the sampled frames |
| `runMesh.sh`, `runCase.sh`, `runPost.sh`, `runViz.sh` | PBS jobs, submitted from the case: re-mesh it; run it; post-process it; draw it |

The Python behind the new parts is `plumetools/cai2012ext/`. It imports
`plumetools.cai2012`, `plumetools.nozzleplume` and `plumetools.rarefaction`
and edits none of them.

---

## The resolution diagnostics

**Before meshing** — `generate_cases.py`, recorded in `manifest.yaml` under
`feasibility`. The nozzle-plume family's pre-run checks, by its names and with
its default thresholds where the cai2012 schema has no key, plus this family's
mesh-budget check:

| check | threshold | Kn 10 | Kn 1 | Kn 1e-3 | Kn 1e-4 |
|---|---|---|---|---|---|
| `cell_over_mean_free_path` | ≤ 1 | 0.005 | 0.05 | **26.2** | **270** |
| `mesh_budget` (cells for `dx = λ0` ≤ 4e6) | | met | met | **1.9e10** | **1.8e13** |
| `time_step_collision_time` (Δt/τ0) | ≤ 1 | 2.7e-4 | 2.7e-3 | **1.42** | **14.6** |
| `time_step_courant` | ≤ 0.5 | 0.20 | 0.20 | 0.20 | 0.20 |
| `transient_length` (start / one jet flight) | ≥ 1 | 17.7 | 17.7 | 17.7 | 17.7 |
| `internal_energy_exchange` | | argon, 0 dof | ← | ← | ← |
| `particles_per_cell` (exit) | ≥ 5 | 40 | 40 | 40 | 40 |
| `parcel_count` | ≤ 3e7 | 2.3e6 | 2.3e6 | 1.6e7 | 1.5e7 |

**After a run** — `Alldiagnose.py` (step 2 of `./AllpostCases`), the same
pipeline as R-4D's `./Allpost` step 4, fed from this `case.yaml` and
`diagnostics.yaml`:

* `<time>/mfp`, `KnD`, `KnCellMax`, `KnCellV`, `cellOverMfpMax`, `cellOverMfpV`,
  `cellSizeMax`, `cellSizeV`, `cellSizeTarget`, `KnGradN`, `KnGradT`,
  `mfpResolved`, `parcelsPerStep`, `sampledParcels`, `diagValid`,
  `diagValidGradN`, `diagValidGradT`, `diagStatus` — volume fields beside
  `rhoNMean`, open in ParaView;
* `results/diagnostics_summary.yaml` — validity, Δ/λ statistics, the worst
  cells, a λ/3 mesh estimate, and `exit_plane`: the measured `KnD` and Δ/λ in
  the exit cells beside the target Kn and the planned cell/λ0;
* `results/centreline_diagnostics.csv` / `.png` — along the axis, on
  `centerline.csv`'s stations;
* `results/viz/diagnostics/*.png` — every field on the centre plane
  (`viz-diagnostics.yaml`, needs pvpython);
* columns of `results/study-table.csv` — valid volume fraction, exit KnD and
  Δ/λ, the volume fraction coarser than λ/3, the cells a λ/3 mesh needs.

**Read `diagValid` first.** A cell's diagnostics are trusted only with *more
than 20* parcels per step (`diagnostics.yaml`, the nozzle-plume default). At
40 parcels in the exit cell that holds near the exit and not much further: on
the synthetic smoke case (docs §10), 1% of the cells and 11% of the volume. Elsewhere the
values are still written and drawn, flagged `diagValid = 0`. Lowering
`min_sampled_particles` or switching `particle_count_basis` to `accumulated` is
a choice `docs/rarefaction-diagnostics.md` explains the cost of; it is left at
the default.

---

## On the laptop: generate, write dictionaries, validate (no OpenFOAM)

From WSL, with the repository's environment active:

```bash
cd cases/cai2012-kn-extension
./generate_cases.py --dry-run      # the plan, the 4 verdicts, the cost; writes nothing
./generate_cases.py                # create Cases/ and manifest.yaml
./AllmeshCases --dict-only         # every case's OpenFOAM dictionaries, no OpenFOAM
./validate_cases.py                # exit 0 = the tree IS the study
```

`generate_cases.py` is safe to run any number of times:

| on disk | outcome |
|---|---|
| nothing | created |
| identical to what would be generated | **left alone**, mesh and results included |
| a runner script, `diagnostics.yaml` or `meta` differs (`stale`) | left alone and reported; `--refresh` copies the new files in |
| an input in `case.yaml` differs (`conflict`) | **nothing touched**, exit 1; `--regenerate` replaces it if it holds no solver output, `--overwrite` even if it does (destroying those results) |

A stray directory under `Cases/` is reported and never scheduled; one holding
an excluded Kn also fails `validate_cases.py`.

---

## On the cluster: PBS, one job per case

Every job script `cd`s to `$PBS_O_WORKDIR` (and fails if it is unset or the
wrong kind of directory), so **submit each from the directory it names**. Each
sources [`hpc_env.sh`](hpc_env.sh), found by looking in `./`, `../` and
`../../`: the Spack/Lmod modules, OpenFOAM v2512's `bashrc`, the repository's
`.venv`, and ParaView 6.2.0 for the images (`PVPYTHON_PATH`, kept if already
set: `qsub -v PVPYTHON_PATH=…`). The site's paths, modules and mail address
live there and in the `#PBS` lines.

1. **Commit and push** the repository, including `plumetools/cai2012ext/` and
   this directory. `Cases/` and `manifest.yaml` are gitignored: they are
   generated on the cluster, so the manifest records the cluster's commit.

2. **Once, on the cluster**, from a clone at the path `hpc_env.sh` names: build
   the custom inflow library (`applications/dsmcBoundaryModels/Allwmake`, or
   `./Allsetup`) and install the repository into the `.venv`:

   ```bash
   pip install -e ".[test,plots]"
   ```

3. **Generate, validate and mesh** all four, one job (`shortq`, 1 core, 4 h):

   ```bash
   cd cases/cai2012-kn-extension
   qsub runMeshAll.sh
   ```

   It runs `generate_cases.py`, `validate_cases.py --inputs-only`,
   `AllmeshCases` and `validate_cases.py` (with the dictionaries now written),
   and any failure fails the job. Read its output before going on. To re-mesh
   one case: `cd Cases/<case> && qsub runMesh.sh` (`shortq`, 1 core, 2 h).

4. **Run** each case as its own job (`longq`, one node, 192 cores, 120 h), all
   four, the dense two included:

   ```bash
   for case in Kn10_np1x Kn1_np1x Kn0p001_np1x Kn0p0001_np1x; do
       (cd Cases/$case && qsub -N knx-run-$case runCase.sh)
   done
   ```

   `runCase.sh` refuses a case at Kn = 0.01, 0.1 or 100; re-assesses the case
   from its own `case.yaml` and, if it is unsuitable for production, writes a
   `NON-PRODUCTION` marker into it (and the same text at the top of
   `output.log`); runs `./Allrun` into `output.log` / `error.log`; and records
   the wall clock in `results/cost-model.yaml`. The 192 ranks are
   `dsmc.n_subdomains` in `case.yaml`, which `./Allrun` reads from
   `system/decomposeParDict`; no job passes a rank count.

   **A job that runs out of walltime is resubmitted by hand**, with the same
   `qsub`: `./Allrun` continues from the newest time directory, `processor*/`
   included. A case that had already reached `endTime` runs nothing and records
   nothing. Nothing resubmits itself.

5. **Post-process and draw** each case once its run has finished (`longq`, one
   node, 192 cores, 120 h each):

   ```bash
   cd Cases/<case> && qsub runPost.sh    # then, when it is done:
   cd Cases/<case> && qsub runViz.sh
   ```

   `runPost.sh` is `../../AllpostCases --case <case> --no-viz --no-analyse` —
   `Allpost`, `Alldiagnose.py`, `Allcompare.py`, `Allaudit.py`,
   `Allconverge.py` — into `output_post.log` / `error_post.log`. `runViz.sh`
   is `--viz-only`: the flow-field series over the sampling window and the
   diagnostic fields `runPost.sh` wrote, into `output_viz.log` /
   `error_viz.log`.

6. **When every case is post-processed**, the study table and figure (`shortq`,
   1 core, 1 h):

   ```bash
   cd cases/cai2012-kn-extension && qsub runAnalyse.sh
   ```

The `#PBS -N` names are generic — `knx-meshall`, `knx-mesh`, `knx-run`,
`knx-post`, `knx-viz`, `knx-analyse` — so pass `qsub -N <name>` to tell one
case's jobs from another's. Repository scripts are committed without the
executable bit on Windows, so the jobs call every one of them with `bash` or
`python`.

**Disk.** Nothing is pruned: every time directory and every parcel cloud is
kept, decomposed under `processor*/` and again reconstructed. About **2.7 TB**
for all four (below); check the scratch quota before submitting.

The same driver, interactively (OpenFOAM only for `writeCellCentres` inside
`Allpost`; ParaView optional):

```bash
bash ./AllpostCases                  # every case with output, then analyse.py
bash ./AllpostCases --case Kn1_np1x
bash ./AllpostCases --no-diagnostics # skip Alldiagnose.py and its images
bash ./AllpostCases --no-analyse     # the per-case chain only
bash ./AllpostCases --analyse-only   # the study table and figure only
```

The rarefaction diagnostics read the whole mesh into memory: 450 MB for
R-4D's 0.46-million-cell mesh, so roughly 1.3 GB for the 1.3-million-cell
Kn ≥ 1 meshes and ~4 GB for the 3.9-million-cell dense ones —
comfortable on a node, tight on a laptop's WSL.

**The dense cases run, labelled.** Kn = 1e-3 and 1e-4 are unsuitable for
production. Their `runCase.sh` runs them anyway and writes the
`NON-PRODUCTION` marker, listing the unresolved criteria; every artefact from
them says so, and nothing from them is a validated result.

---

## What it will cost

Estimated by `generate_cases.py` from `cases/cai2012-health`'s nine measured
runs (12 MPI ranks on the development machine). Those were collisionless
(Kn = 100), so for the collisional cases these are **lower bounds** at
12 ranks. Each case's `runCase.sh` records this family's own wall clock in
`results/cost-model.yaml` (under a file lock: jobs can finish together), and
the estimates switch to it as soon as one case has run.

| | parcels (est.) | steps | ≥ wall clock @12 ranks | disk, every write and cloud kept |
|---|---|---|---|---|
| `Kn10_np1x`, `Kn1_np1x` | 2.3e6 each | 45,924 | 26 h each | ~212 GB each |
| `Kn0p001_np1x`, `Kn0p0001_np1x` | 1.6e7, 1.5e7 | 88,288, 85,974 | 289 h, 261 h | ~1,160 GB, ~1,070 GB |
| **all four** | | | **~603 h** | **~2.7 TB** |

The disk is dominated by the health family's write cadence (1/8 domain
transit) running through Cai's 17.7-transit startup: 178 writes per case, of
which only the last 37 are sampled — each write's fields and parcel cloud
held twice, under `processor*/` and reconstructed. See the documentation, §6.

---

## Known limitations, in one place

* **Kn = 1e-3 and 1e-4 are unresolvable within `mesh.max_cells = 4e6`.** They
  need ~2e10 and ~2e13 cells for `dx = λ0`; the budget gives 26× and 270×
  coarser cells, and Δt is 1.4× and 15× the exit collision time. Raising
  `mesh.max_cells` per Kn in `study.yaml` is a configuration change the policy
  allows; the cell counts say why it is not a practical one.
* **One parcel weight.** Statistical adequacy is judged from each run's own
  output (`Allaudit.py`, `Alldiagnose.py`'s `diagValid`, `Allconverge.py`),
  not from runs at other weights — and at this weight the diagnostics' own
  validity rule holds near the exit only.
* The **cost model** has no collision work in it.
* The startup is **Cai's 10⁴ reference collision times at Kn = 0.01** for every
  case: 17.7 domain transits.
* The **Cai figure numbers** for the angular and off-axis comparisons could not
  be verified — the paper is not in the repository. See the documentation, §8.

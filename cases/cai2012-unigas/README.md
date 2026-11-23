# cases/cai2012-unigas — the Cai Knudsen sweep on uniGasFoam

The seven cases of [`cases/cai2012-kn-extension`](../cai2012-kn-extension/README.md)
(Kn = 100 … 1e-4, Cai's circular argon plume, D = 0.2 m, S₀ = 2, T = 300 K),
run on a patched uniGasFoam's **hybrid USP-DSMC** solver instead of standard
`dsmcFoam`. The cases are the same case for case: the same `case.yaml`
outside `meta`, the same meshes, the same particle weight and the same run
window. `validate_cases.py` checks all of it against the kn-extension. Only the
solver changes.

Design: [`docs/superpowers/specs/2026-11-19-cai-kn-sweep-unigasfoam-design.md`](../../docs/superpowers/specs/2026-11-19-cai-kn-sweep-unigasfoam-design.md).
Tutorial: [`docs/tutorials/cai2012-unigas.md`](../../docs/tutorials/cai2012-unigas.md).

> **Nothing here is validated.** No production run of this family is recorded.
> The pilot (`runPilot.sh`) measures whether 3-D production is affordable;
> going ahead is the owner's decision on its report.

## Files

| file | what |
|---|---|
| `study.yaml` | the kn-extension's matrix, verbatim below `kn_cases`; `meta.unigas_commit` pins the uniGasFoam fork commit |
| `baseCase/case.yaml` | the kn-extension's template, unchanged outside `meta.study` |
| `baseCase/unigas.yaml` | the uniGasFoam settings (spec §4), every key explicit |
| `generate_cases.py`, `validate_cases.py` | idempotent generation; the equivalence check |
| `runBuild.sh`, `runMeshAll.sh`, `runPilot.sh` | PBS jobs: build the fork, mesh every case, run the pilot |
| `Cases/<case>/runCase.sh`, `runPost.sh` | per-case PBS jobs |

## Commands

```bash
python generate_cases.py --dry-run     # the plan; writes nothing
python generate_cases.py               # Cases/ and manifest-unigas.yaml
python validate_cases.py               # equivalence with the kn-extension
qsub runBuild.sh                       # cluster: build the pinned fork commit
qsub runPilot.sh                       # Kn = 0.001, 2 transits: cost and resolution
qsub runMeshAll.sh                     # after go: mesh all seven
cd Cases/Kn1_np1x && qsub runCase.sh   # Kn = 0.1 and 1 first (the cross-check)
```

## What differs from the kn-extension, and why

* **The verdict is made after the run** (spec §7). USP cells are larger than
  the mean free path and are not collision-time limited *by design*, so the
  kn-extension's pre-run cell/λ₀ and Δt/τ₀ gates do not apply. Every case is
  `post-run-verdict-pending` until `Allverdict.py` has judged its DSMC cells.
  The pre-run assessment is kept in the manifest as
  `dsmcfoam_prerun_feasibility`, for reference only.
* **Analysis output before the verdict says so.** `Allconverge.py` runs before
  `Allverdict.py` (the verdict needs its stationarity result), so
  `results/stationarity.yaml` carries `post-run-verdict-pending` as its
  production status. `Allcompare.py` runs after the verdict and shows it.
* **The convergence series is a per-write-interval average.** dsmcFoam's
  single-step `rhoN` has no uniGasFoam equivalent. The `ArWindow` field averages
  over each write interval and is converted to `rhoN`. It is less noisy than a
  single step, so `Allconverge.py`'s trend test has more power than in the
  dsmcFoam suite.
* **Only the newest particle cloud is kept** (uniGasFoam's default). A job killed
  while writing can leave a rank whose older cloud is already deleted.
  `Allrun` refuses to resume such a case, because uniGasFoam would silently re-fill
  that subdomain, and prints the manual restart. If no write with complete clouds
  remains, the case starts again.
* **Two run legs.** The startup leg writes once per domain transit, and the
  sampling leg writes at the health cadence. This fixes the kn-extension's write
  cadence open problem (§11 item 2) without touching the frozen code.

## Settings fixed by the owner

* Full 3-D Cartesian (no wedge), all seven Kn on uniGasFoam.
* uniGasFoam fork branch `dsmc-time-step-limit`: DSMC cells respect the
  collision-time limit (`dsmcCellsTimeStepLimited true`) and the sub-cell cap
  is configurable (`maxSubCellLevels`, default 10).
* Cell weighting **off**: one particle weight, as in the dsmcFoam suite, and the
  same near-vacuum start.

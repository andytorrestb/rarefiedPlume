# Tutorial: the Cai Knudsen sweep on uniGasFoam

The seven cases of the [Kn extension](cai2012-kn-extension.md) (Kn = 100 … 1e-4),
run on a patched **uniGasFoam** with its hybrid USP-DSMC solver instead of
standard `dsmcFoam`. The cases are the same case for case: physics, meshes,
particle weight and run window. Only the solver changes.

* Family: [`cases/cai2012-unigas`](https://github.com/andytorrestb/rarefiedPlume/tree/main/cases/cai2012-unigas)
* Design: `docs/superpowers/specs/2026-11-19-cai-kn-sweep-unigasfoam-design.md`
* Code: `plumetools/cai2012ugf/`

!!! warning "Not validated"
    No production run of this family is recorded. The pilot decides whether
    3-D production is affordable, and going ahead is the owner's decision.

## Why uniGasFoam

The dsmcFoam suite cannot resolve Kn ≤ 0.01 inside its 4e6-cell budget:
cell/λ₀ reaches 270 and Δt/τ₀ reaches 14.6. uniGasFoam (Vasileiadis et al.,
*Comput. Phys. Commun.* 310, 2025) handles cells that are near-continuum
(gradient-length Knudsen number Kn_GLL ≤ 0.05) with the unified stochastic
particle (USP) BGK scheme. USP cells may be larger than λ, and their time step is
limited only by Courant. The rest stay DSMC.

## The uniGasFoam fork

The family runs the owner's fork, branch `dsmc-time-step-limit`, at the commit
pinned in `study.yaml` (`meta.unigas_commit`). The branch adds two opt-in
`adaptiveProperties` keys:

| key | default | why |
|---|---|---|
| `dsmcCellsTimeStepLimited` | `false` | With the USP model, published uniGasFoam limits **no** cell by its collision time, DSMC cells included. `true` limits every DSMC cell, as the paper describes. This family sets it. |
| `maxSubCellLevels` | 10 | the sub-cell cap per direction, previously hard-coded |

On the plume tutorial (OpenFOAM v2512) the keys at their defaults reproduce
the published behaviour, within run-to-run noise: the random generator is seeded
from the clock. With the limit on, Δt/τ in DSMC cells stayed at or below 0.11.
Without it, Δt/τ reached 4.

The branch also fixes a crash, with no key. At every adaptation (each
`adaptationInterval` steps) uniGasFoam fills empty cells with the smallest
occupied density and temperature. It took that minimum per processor, so a
processor holding no particle filled its cells with VGREAT (10³⁰⁰), and the
collision rate overflowed: a floating-point exception at step 20. Every case
here starts from near vacuum, so every rank away from the nozzle was empty at the
first adaptation. The minimum is now taken over all processors, as a serial
run already took it, and serial runs are unchanged.
`tests/openfoam/test_cai2012ugf_empty_rank.py` reproduces the crash on two
ranks.

Build it locally, with OpenFOAM v2512 sourced:

```bash
cd ../uniGasFoam && git switch dsmc-time-step-limit && ./Allwmake 8
git rev-parse HEAD > "$FOAM_USER_APPBIN/uniGasFoam.commit"   # what Allrun checks
```

On the cluster: `cd cases/cai2012-unigas && qsub runBuild.sh` checks out the
pinned commit, builds it and records it.

## Settings

`baseCase/unigas.yaml` states every solver setting (spec §4):

* hybrid USP-DSMC, NTC and VHS; Kn_GLL breakdown at 0.05;
* sub-cells ≤ λ/2, Δt ≤ τ/5 in DSMC cells, Courant 0.5;
* **cell weighting off**: one particle weight and the same near-vacuum start as
  the dsmcFoam suite.

The particle weight, mesh and run window come from `case.yaml`, which is the
kn-extension's.

## Commands

```bash
cd cases/cai2012-unigas
python generate_cases.py --dry-run     # the plan
python generate_cases.py               # Cases/ + manifest-unigas.yaml
python validate_cases.py               # the cases ARE the kn-extension's
```

Then, on the cluster:

1. `qsub runBuild.sh`: build the pinned fork commit.
2. `qsub runPilot.sh`: Kn = 0.001 on its production mesh, 2 transits.
   Read `Pilot/Kn0p001_np1x/results/pilot_report.md`.
3. **Go/no-go (the owner).** The report gives:
   * the time step and what sets it (DSMC collision times, or Courant);
   * the DSMC/USP split;
   * the effective collision spacing in DSMC cells;
   * the particle count and rank imbalance;
   * disk per write and a projected wall time.

   The levers if it's not affordable: more particles for the dense cases, an
   axisymmetric wedge, or accepting non-production labels.
4. `qsub runMeshAll.sh`, then `cd Cases/<case> && qsub runCase.sh`. Run Kn = 0.1
   and 1 first: they are the cross-check against the dsmcFoam suite.
5. After each run: `qsub runPost.sh` in the case.
6. `python -m plumetools.cai2012ugf.crosscheck` and
   `python Allscale.py inventory | fit --full | report`.

## How a run works

`Allrun` runs two legs:

* **startup**, to `average_start_s` (Cai's 10⁴ t₀), writing once per domain
  transit;
* **sampling**, 4.5 transits at the health cadence.

Before each leg it asks `python -m plumetools.cai2012ugf.hpc resume-plan` what
to do. That plan refuses to resume a case where a rank that held gas has lost its
particle cloud. uniGasFoam keeps only the newest cloud, and would silently
re-fill that subdomain. The refusal prints the manual restart.

After the last leg, `Allrun` reconstructs and converts uniGasFoam's fields to
the dsmcFoam names the analysis reads (`rhoNMean`, `momentumMean`, …). The
conversion is exact: through `plumetools.cai2012.post.read_case` it returns
uniGasFoam's own density, velocity and translational temperature.

## Post-processing and the verdict

`Allpostugf` runs these steps in order:

1. convert;
2. `Allpost`;
3. `Alldiagnose.py`;
4. `Allaudit.py`;
5. `Allconverge.py`;
6. **`Allverdict.py`**;
7. `Allcompare.py`.

Every analysis script except `Allverdict.py` is the kn-extension's.
`Allcompare.py` and `Allconverge.py` differ only in importing this family's post
context. `validate_cases.py` checks that.

The **verdict is made after the run** (spec §7), because USP cells fail the
dsmcFoam suite's pre-run checks by design. Over the DSMC cells that hold
particles, it requires:

* effective collision spacing / λ ≤ 1. The *effective* spacing allows for the
  fact that a sub-cell with one particle draws its partner from the whole cell;
* Δt/τ ≤ 0.2;
* at least 5 particles in every DSMC cell above n/n₀ = 1e-3.

It also requires a demonstrated stationary window, and no known source
limitation (Kn ≤ 0.001 carries the direct-exit inflow-surface one). Cells where
uniGasFoam writes its `GREAT` sentinel (no defined λ or τ) are counted as
undefined, not judged.

## Local end-to-end check

```bash
python -m pytest tests/openfoam/test_cai2012ugf_smoke.py -m needs_openfoam
python -m pytest tests/openfoam/test_cai2012ugf_empty_rank.py -m needs_openfoam
```

The first runs a ~38k-cell Kn = 0.01 case on two ranks for 0.3 transits, through
every step above. It checks the plumbing, not the physics. The second splits
the same mesh by hand so that one rank holds no particle at the first
adaptation, and checks that uniGasFoam runs on past it (the fork fix above).
Both need the project's Python first on `PATH` (`Allmesh` and `Allrun` call
`python`).

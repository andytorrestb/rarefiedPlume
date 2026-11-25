# Plate-unigas: particle-weighted decomposition and the collisionless density comparison (short design)

Approved by the owner on 2026-11-25 in conversation. Extends
[`2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md`](2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md).

## Why

Scotch balancing cells put 90% of the smoke run's particles on one rank. Cells
crowd near the nozzle (sized by the mean free path) while particles spread
through the plume and pile up in front of the plate. Run time follows the
busiest rank. The owner asked for a particle-weighted split "where deemed
necessary", before the 21 cases go to the cluster.

The owner also ruled that the collisionless model differing from the collisional
simulation is "a feature and central point of this analysis", so the same model
joins the Cai comparison.

## What was measured first (Kn_e = 100 zD100 prototype mesh, scotch, busiest rank ÷ mean)

| ranks | cells (as today) | Cai plume only | Cai + plate gas | measured particles |
|---|---|---|---|---|
| 12 | 8.7× | 2.6× | 1.31× | 1.30× |
| 48 | 24.5× | 4.4× | 1.87× | 1.60× |

At 192 ranks the plain split left 87 ranks empty. The laptop mesh is too coarse
to judge the weighted splits at 192 ranks: its largest cell holds about one
rank's share. Cai's plume alone predicts 2–3× too few particles: the plate's
re-emitted gas is missing. With that gas added, the prediction matches the
measured particles cell by cell (median ratio 1.00–1.03).

## Approaches

1. **Predicted weights, Cai plume only.** No run interruption, but imbalance
   stays at 2–4×.
2. **Predicted, then one measured re-split after 2 transits.** As even as it
   gets, at the cost of one serial reconstruct/re-decompose stop. This was my
   first pick.
3. **Predicted weights with the plate's re-emitted gas modelled (chosen by
   the owner).** No stop. As even as the measured split at Kn_e = 100. Blind to
   collisions at low Kn_e, which the owner accepts.

## Design

- `plumetools/plateugf/collisionless.py`: `Model(cfg)`.
  - **Plume:** Cai's plume, zero in the plate's shadow.
  - **Plate gas:** the upstream face re-emits Cai's incident flux diffusely
    at the wall temperature. For each 0.25 m patch, this is a closed-form
    solid-angle sum, exact for a piecewise-constant flux. It covers the whole
    plate, mirrored quarters included.
  - **Large meshes:** above 0.5M points in front of the plate, the plate gas
    is tabulated on a (y, z, log h) grid and interpolated.
- `plumetools/plateugf/weights.py`: `0/decompositionWeight` = 1 + n V / W,
  where W is `nEquivalentParticles` (my call: a cell costs about one
  particle).
  - `decomposeParDict` gets `weightField decompositionWeight;`.
  - `Allmesh` writes the weights. `Allrun` rewrites them if missing and copies
    them to any later time it re-decomposes.
- All 21 cases use it. Where cells outnumber particles (Kn_e = 1e-4), the
  weights fall back toward a cell split by themselves.
- Comparison (`plateugf viz`):
  - slices of the collisionless density, and of simulated ÷ collisionless
    density, on both symmetry planes;
  - `results/cai_density.yaml`, giving simulated ÷ collisionless molecules
    near the nozzle, mid-plume, before the plate and beside/behind it;
  - a `family-table` column, "n front".
- Testing, cut back at the owner's request:
  - unit tests for the model (wall limit, shadow, grid against the exact sum),
    the weight file, and every `plateugf` call in the case scripts;
  - the comparison run on two prototypes;
  - one short laptop run: weights, then a weighted `decomposePar`, then
    uniGasFoam, then `reconstructPar`.

Found on the way: `Allmesh` carried a literal `\n` that made `flatten-symmetry`
fail in every case's mesh job. It is fixed, and the script-call test now
catches that class of error.

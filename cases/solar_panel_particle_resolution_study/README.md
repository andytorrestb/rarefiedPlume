# cases/solar_panel_particle_resolution_study — legacy particle-weight study

Three copies of `iss_solar_panels_2`'s case (`Fnum_1e16`, `Fnum_5e15`,
`Fnum_2.5e15`) that differ only in `nEquivalentParticles`, with plotting scripts
for sampled plume pressure (commit `3b7b5e4`). They name the MNF fork's
`dsmcFoam+`, have no run script, and every plotting script reads a CSV that was
never committed (`*.csv` is gitignored), so **neither the runs nor the figures
can be reproduced** from the repository. No result here is validated.

> **Do not run any `Fnum_*/clean.sh`.** It is the unguarded `rm -r 0* …` of
> finding RP-01.

Variant table, the exact missing inputs, and a figure placeholder:
[`docs/tutorials/particle-resolution.md`](../../docs/tutorials/particle-resolution.md).

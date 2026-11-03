# Solar-panel particle-resolution study (`Fnum_*`)

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Cases** | `cases/solar_panel_particle_resolution_study/Fnum_1e16`, `Fnum_5e15`, `Fnum_2.5e15` |
    | **Purpose** | Sensitivity of sampled plume pressure to the DSMC particle weight (commit `3b7b5e4`, 2023-02-28) |
    | **Model** | Per-case copies of the legacy `processInflowData.py`, as in the [ISS cases](iss-solar-panels.md) |
    | **Evidence** | None. The committed figures cannot be regenerated: their CSV inputs were never committed ([blockers](#blockers)) |
    | **Status** | Labelled *active — under refactor*; **not runnable or re-plottable** with what is in the repository |
    | **Start in** | `cases/solar_panel_particle_resolution_study/` (plotting script) or a `Fnum_*` directory |

## Purpose

Three copies of one case that differ in `nEquivalentParticles` — the number of
real molecules each simulated particle represents. Halving it doubles the
particle count and lowers statistical noise; agreement between levels is the
usual evidence that sampled quantities are converged in particle weight.

The study samples pressure along lines and on slices at several heights
(5, 7.5, 10 and 20 m in the committed image names). Despite the study's name,
the mesh has **no panel patch**: it is `iss_solar_panels_2`'s mesh, with only
`inflow` (2044 faces, `type patch`) and `vacuum` (17 378 faces, `type patch`).
The "panel" is a sampling location in the free plume.

## Variant table

| | `Fnum_1e16` | `Fnum_5e15` | `Fnum_2.5e15` |
|---|---|---|---|
| `nEquivalentParticles` | 1e16 | 5e15 | 2.5e15 |
| `controlDict` | `dsmcFoam+`, `endTime` 1e-2 s | same | same |
| `processInflowData.py` | `To` = 3000 K, `Po` = 744634 Pa ("108 psi"); density still uses 475 psi and 300 K (SM-01, SM-02) | identical | identical |
| Tracked extras | plotting scripts | plotting scripts | plotting scripts, 40 `plume_pressure_h_*` PNG frames, `pressure_over_a_line/*.png` |

The three `processInflowData.py` files are byte-identical; `system/` differs
only in how `endTime` is written. SM-02 matters most here: with `To` = 3000 K
for velocity and temperature but 300 K for density, mass and energy flux are
inconsistent by about a factor of ten ([source-flow model](../source-flow-model.md#4-known-defects-frozen-deliberately)).

## Blockers

1. **The solver is the MNF fork** (`application dsmcFoam+;`), as for the
   [ISS cases](iss-solar-panels.md#blockers), with the same missing run script
   and the same `processInflowData.py` prerequisites (`log.checkMesh`,
   `constant/polyMesh/sets/inflow`).
2. **The plotting inputs do not exist.** Every plotting script reads a CSV:

    | Script | Reads |
    |---|---|
    | `graph_resolution_study.py` | `Fnum_{2.5e15,5e15,1e16}/line_h_7-5_offset_1-1.csv` (columns `Points_1`, `p_Ar`) |
    | `Fnum_*/process_slice.py` | `slice_h_10m.csv`; writes `pressure_distribution_h_10m.csv` |
    | `Fnum_*/pressure_over_a_line/graph_impingement_pressures.py` | `width_05m.csv`, `width_10m.csv`, `width_20m.csv` |
    | `Fnum_*/displacement_graph_data/graph_displacement.py` | `plume_pressure_distribution.csv`, `uniform_pressure_distribution_{max,avg,min}.csv` |

    `*.csv` is gitignored repository-wide and no CSV is tracked (finding RP-06 in
    [Environment](../environment.md)). The CSVs were presumably ParaView exports;
    which filter, time and field produced them is not recorded.

!!! danger "Do not run `clean.sh`"
    Each `Fnum_*/clean.sh` is the unguarded `rm -r 0* …` of RP-01; see the
    [ISS tutorial](iss-solar-panels.md#blockers).

## What you can do safely

Read-only comparison of the variants, from the repository root:

```bash
cd cases/solar_panel_particle_resolution_study
diff Fnum_1e16/constant/dsmcProperties Fnum_5e15/constant/dsmcProperties
```

To make the study reproducible the missing pieces are: a `plumetools`
conversion of the inflow step (as `3d-inflow` has), a run script for the
standard solver, and a recorded sampling step that writes the CSVs the plotting
scripts expect. None exists; see [Known gaps](../development/known-gaps.md).

## Outputs

Historical, committed: `Figure_1.png` at the study root (`graph_resolution_study.py`
only calls `plt.show()`, so the file was saved by hand from the plot window,
under matplotlib's default name), per-case
`centroids.png`/`norms.png`/`sym_theta.png`, the `Fnum_2.5e15` pressure frames.
Their provenance — solver build, sample times, sampling geometry — is not
recorded, so they are not evidence of convergence.

## Resources

Configured only: `endTime` 1e-2 s, `numberOfSubdomains` 8 in each case's
`system/decomposeParDict`. No run time or memory is recorded.

## Figure placeholder

!!! figure "Figure placeholder — pressure along a sampling line, three particle weights"
    **Subject:** `p_Ar` against arc length on the line `line_h_7-5_offset_1-1`, one curve per `Fnum`.
    **Source run:** all three `Fnum_*` cases, solved and sampled; not possible today ([blockers](#blockers)).
    **Plotting command:** `python graph_resolution_study.py` in `cases/solar_panel_particle_resolution_study/` once the three CSVs exist (read from source, not run; it opens a window with `plt.show()` and saves nothing, so save the figure by hand).
    **Destination:** `docs/assets/figures/particle-resolution-line-pressure.png`
    **Suggested caption:** "Sampled plume pressure at h = 7.5 m for three particle weights (state solver build and sampling time)."

## Related

- [ISS solar panels](iss-solar-panels.md) — the parent cases and the shared blockers
- [Rarefaction diagnostics](../rarefaction-diagnostics.md) — the supported way to judge particle and mesh resolution in the converted families
- [Case catalog](../reference/case-catalog.md)

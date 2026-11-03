# Figures to add

This directory holds the documentation's figures. **It contains none yet.**
Every tutorial carries an empty, labelled placeholder instead, so no page
shows a figure that was not produced by a real, identified run, and no page
links to an image that does not exist.

To add one: produce it from the run named in its placeholder, save it here
under the listed name, and replace the placeholder admonition in the page with
an image and a caption that states the case, the time window or sample, the
git commit and anything non-production about it (for example
`![Centreline density, Kn100_np1x …](../assets/figures/cai2012-centerline-density.png)`).
Do not add analytical or synthetic plots in their place.

| Destination | Subject | Source run | Page |
|---|---|---|---|
| `3d-inflow-boundary-density.png` | legacy source-flow number density on the inflow hemisphere | `3d-inflow` `./Allrun --no-solve` on the restored Pointwise mesh | [3d-inflow](../../tutorials/3d-inflow.md#figure-placeholder) |
| `cai2012-centerline-density.png` | centreline n/n₀, DSMC against Cai (same-cell average) | `cai2012` `Kn100_np1x`, `Kn0p1_np1x` + `Allpost` | [cai2012](../../tutorials/cai2012.md#figure-placeholder) |
| `cai2012-health-contact-final.png` | the 3 × 3 `dsmcRhoN` contact sheet at each case's final frame | all nine health cases + `AllpostCases` | [cai2012-health](../../tutorials/cai2012-health.md#figure-placeholder) |
| `kn-extension-knudsen-comparison.png` | departure from Cai across Kn | the candidate Kn-extension cases + `analyse.py` | [Kn extension](../../tutorials/cai2012-kn-extension.md#figure-placeholder) |
| `kn-interface-overlap.png` | overlap region, truncated vs full-domain density | a `val_Kn0p01_*` case that is not infeasible (none exists yet) | [Interface track](../../tutorials/cai2012-interfaces.md#figure-placeholder) |
| `kn-scaling-fit.png` | fitted coordinate map against sampled density | `Allscale.py fit --full` on real fields | [Field scaling](../../tutorials/cai2012-scaling.md#figure-placeholder) |
| `virtual-collision-contours.png` | three-way density contours | `compare_exports.py --kn 10 1` on the saved exports | [Virtual collisions](../../tutorials/cai2012-virtual-collision.md#figure-placeholder) |
| `hybrid-breakdown-surface.png` | breakdown indicators and the interface on the meridional plane | hybrid `smoke-kn0p01` after `cfd` and `post-cfd` | [Hybrid](../../tutorials/cai2012-hybrid.md#figure-placeholder) |
| `markelov1999-windward-pressure.png` | windward plate pressure, four reservoir pressures | all four Markelov cases + `AllpostCases` | [markelov1999](../../tutorials/markelov1999.md#figure-placeholder) |
| `r-4d-plate-pressure.png` | wall-normal pressure on the plate's upstream face | a full-length `r-4d/baseCase` run + `Allpost` | [r-4d](../../tutorials/r-4d.md#figure-placeholder) |
| `iss-solar-panels-3-panel-pressure.png` | pressure on the notional panel | none possible with the supported stack | [ISS](../../tutorials/iss-solar-panels.md#figure-placeholder) |
| `particle-resolution-line-pressure.png` | line pressure, three particle weights | the three `Fnum_*` cases, solved and sampled (blocked) | [Particle resolution](../../tutorials/particle-resolution.md#figure-placeholder) |

Formats: PNG for raster output as the families write it; keep the file the
family's own post-processing produced rather than redrawing it, so the figure
and the run stay traceable.

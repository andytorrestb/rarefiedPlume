# r-4d: an R-4D exit-plane plume, with or without a plate

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | The exit plane of an R-4D nozzle (D = 0.377 m, 3187 m/s, 498 K, 1.3 × 10²² m⁻³; Kn ≈ 4 × 10⁻⁴) as the DSMC inlet, expanding along +x into vacuum onto a 34.14 × 11.88 × 0.25 m plate — or, in one variant, with no plate |
    | **Model** | The general nozzle-plume family, `plumetools.nozzleplume`: a uniform exit state given directly, any single VHS gas, an optional plate, a `snappyHexMesh` mesh, `plumeFieldInflow` |
    | **Evidence** | Unit-tested; **solver-tested** in `tests/openfoam/test_nozzleplume_openfoam.py`, including a collisionless conservation check on the plate load. One full-length run recorded on an earlier mesh. **No result validated** against data or another code |
    | **Limitation** | Under-resolved by necessity at this exit density (exit cell ≈ 220 mean free paths), and it says so on every `Allmesh` |
    | **Runs on** | Python for `--describe` and `--dict-only`; OpenFOAM v2512 + `plumeFieldInflow` for mesh and run; `pvpython` for imagery |
    | **Start in** | `cases/r-4d/<variant>/` |

Everything about the case — what was given, what was assumed and why, the
mesh, why it is under-resolved, the cost and the knobs, every output field, and
what has been checked on which mesh — is in
[r-4d, explained](../r-4d-case.md). The family's design (the Cai problem with
the paper taken out, and why it does not import `cai2012`) is in
[`plumetools.nozzleplume`](../reference/api/nozzleplume/index.md).

## Variants

Three tracked directories, each a complete case with its own `case.yaml` and
`All*` scripts:

| | `baseCase` | `baseCaseHPC` | `baseCaseFreePlumeHPC` |
|---|---|---|---|
| Plate | 34.14 × 11.88 × 0.25 m at 20 m, 300 K diffuse wall | same | **none** (`plate: null`) |
| Domain (x_max, y_half, z_half) | 22, 19, 8 m | 50, 40, 16 m | 50, 40, 16 m |
| Background cell | 0.5 m | 1.0 m | 1.0 m |
| Refinement cones | 4 (levels 4/3/2/1 to 1/3/8/20 m) | 5 (levels 6/5/4/3/2 to 1/2/10/20/30 m) | same as HPC |
| `plate_refinement_level` | 2 | 3 | 3 (no plate) |
| `mesh.max_cells` | 3 × 10⁶ | 3 × 10⁷ | 3 × 10⁷ |
| `dsmc.n_equivalent_particles` | 5 × 10¹⁶ | 5 × 10¹⁶ | 1 × 10¹⁸ |
| `diagnostics:` section | present | absent (defaults) | absent (defaults) |
| `dsmc.n_subdomains` | 12 | 12 | 12 |
| Purpose | the case the README describes and the runs were made on | a finer production mesh for a cluster | the free plume at production resolution (≈ 7.35 × 10⁶ cells estimated, per its README) |

The scripts are identical between `baseCase` and `baseCaseHPC`;
`baseCaseFreePlumeHPC` has its own `Allmesh`, `Allpost` and `Allviz`, and
`Allpost` writes no plate-load files for it. None of the three ships a PBS
script. Treat all three as templates: run in a fresh iteration (below), not in
place, once a directory holds results you want to keep.

## Geometry and coordinates

Exit disk at x = 0 on the +x axis; the rest of the x = 0 plane is open vacuum
(molecules deleted); the outer box faces delete. The plate, when present, is
normal to the axis and centred on it (`plate.centre_y_m`, `centre_z_m` move it
off-axis; a tilted plate is not expressible). `domain.symmetry_planes: [y, z]`
models one quarter, with plate totals still reported for the whole plate.

## Prerequisites

| Step | Needs |
|---|---|
| `./Allmesh --describe`, `--dict-only` | Python with `plumetools` |
| `./Allmesh` | `blockMesh`, `snappyHexMesh`, `topoSet`, `createPatch`, `checkMesh` |
| `./Allrun` | `dsmcInitialise`, `dsmcFoam`, `libplumeDsmcBoundaryModels.so`; `decomposePar`, `mpirun`, `reconstructPar` for parallel |
| `./Allpost` | `postProcess` (cell centres and volumes); matplotlib for figures; `pvpython` for the images (skipped with a note if absent, unless `--viz`) |
| `./Allviz` | `pvpython` only (`PVPYTHON_PATH` to name one) |

## Inputs

One file per variant: `case.yaml`, loaded and validated by
`plumetools.nozzleplume.config`. Every dictionary, the mesh and `0/` are
generated from it and gitignored (`cases/r-4d/*/system/`, `constant/`,
`results/`).

| Section | Holds |
|---|---|
| `nozzle` | `exit_diameter_m` |
| `exit` | `temperature_K`; **exactly one** of `number_density_per_m3`, `mass_density_kg_per_m3`, `knudsen`; **exactly one** of `speed_m_per_s`, `speed_ratio` |
| `gas` | one VHS species: `molecular_mass_kg`, `molecular_diameter_m`, `omega` (0.5 = hard sphere), `internal_degrees_of_freedom` |
| `plate` | `standoff_m`, `size_y_m`, `size_z_m`, `thickness_m`, `centre_y_m`, `centre_z_m`, `temperature_K` — or `null` |
| `domain` | `x_max_m`, `y_half_m`, `z_half_m`, optional `symmetry_planes` |
| `mesh` | `background_cell_size_m`, `plate_refinement_level`, `n_cells_between_levels`, `refinement_regions` (cones, or boxes), `max_cells` |
| `dsmc` | `n_equivalent_particles` (written exactly as given), `target_particles_per_exit_cell`, collision and wall models, transient/sampling transits or pinned times, `n_subdomains`, `decomposition_method` |
| `diagnostics` | optional thresholds for the rarefaction diagnostics |
| `execution_profile` | `production` (the committed value) or `smoke` |

Full schema, units and constraints: [configuration reference](../reference/configuration.md#nozzleplume-caseyaml-r-4d).

## Cheap checks (no OpenFOAM)

```bash
cd cases/r-4d/baseCase
./Allmesh --describe     # python -m plumetools.nozzleplume describe . -- prints, writes nothing
./Allmesh --dict-only    # writes system/ and constant/ dictionaries (gitignored); no mesh
```

`--describe` prints the derived exit state (area, mass density, cross-section,
Kn on D, mean free path, speed ratio and Mach, mass flow, thrust, Isp), the
mesh plan with each refinement cone and its half-angle, the estimated cell
count, the time step, the run length, the parcels per exit cell and in the
domain, and the resolution warnings. Run it after every edit to `case.yaml`.

## Production run

```bash
cd cases/r-4d/baseCase
./Allmesh                # dictionaries, blockMesh, snappyHexMesh, nozzle patch, checkMesh, verify
./Allrun                 # 0/ fields, dsmcInitialise, decomposePar, dsmcFoam, reconstructPar
./Allpost                # plate loads, centreline, diagnostics, images -> results/, results/viz/, <time>/
./Allviz                 # slices only, any time: --field, --plane, --time T|all, --dry-run, --list, --output, --spec
```

`./Allrun` takes `--serial`, `-np N`, `--no-solve`; `./Allpost` takes
`--time T`, `--no-viz`, `--viz` (a missing `pvpython` becomes an error),
`--no-diagnostics`. `./Allmesh --force` re-meshes a case with results.

**Restart.** `./Allrun` continues from the newest checkpoint complete on every
rank, as the Cai family does. To run past the end time, raise it in
`case.yaml`, then `./Allmesh --dict-only && ./Allrun`. A run in progress cannot
change particle weight: clean first.

**Watch the parcel count.** `Number of dsmc particles` in `log.dsmcFoam` should
flatten before anything is read; an estimate of how long gas is held up at the
plate is not a substitute.

!!! danger "Resetting"
    `./Allclean` removes solver output (`--mesh` also the mesh; `--dry-run`
    previews). In a directory whose results you want to keep, do not clean and
    do not re-mesh — start a new iteration.

## Change the mesh: a fresh iteration

`dsmcFoam` cannot refine a mesh during a run, and nothing here refines it
automatically. The workflow is: run, inspect the diagnostics, edit by hand, and
start again from t = 0 **in a new directory**:

```bash
cd cases/r-4d
python -m plumetools.iteration baseCase iter01-exit5 --note "exit cone to level 5"
#   copies case.yaml, viz.yaml, open.foam and the All* scripts -- nothing from the run
#   (--include PATH copies another input; repeatable)
cd iter01-exit5
# edit case.yaml
./Allmesh --describe && ./Allmesh && ./Allrun
diff provenance/parent-case.yaml case.yaml     # exactly the edit that made this iteration
```

`plumetools.iteration` only reads the source, refuses a destination that is not
new and empty, and checks that no time or processor directory, parcel or
averaging state was carried over. It writes `iteration.yaml` (parent, note) and
`provenance/`. Keep iterations beside `baseCase` (as `cases/r-4d/<name>`) so the
ignore rules for generated files cover them. Re-meshing a result-bearing case in
place would leave results that address cells that no longer exist.

## Outputs

| Path | Contents |
|---|---|
| `results/plate_loads.csv` | pressure, shear and heat flux on every plate face |
| `results/plate_summary.yaml` | force and heat load on the plate beside the jet's thrust and power |
| `results/plate-pressure.png`, `plate-heat-flux.png` | maps of the upstream face |
| `results/centreline.csv`, `centreline.png` | n, U_x, T along the axis over the exit values, with the cell size of each sample |
| `<time>/mfp`, `KnD`, `KnCellMax`, `mfpResolved`, … | rarefaction diagnostics as volume fields, beside `rhoNMean` |
| `results/diagnostics_summary.yaml`, `centreline_diagnostics.csv`, `centreline-diagnostics.png` | the remeshing summary |
| `results/viz/*.png`, `results/viz/manifest.yaml` | slices of number and mass density, parcels per cell, speed, temperature, the plate loads and the diagnostics |

Plate pressure is `fD · n` from the time-averaged `fDMean`, not `nkT` from the
adjacent cell. The centreline steps down wherever the mesh coarsens — the
sample widening across a peaked profile, not the flow changing. Read
`mfpResolved` before `KnCellMax`: a cell with no data holds −1, which the
renderer draws in the bottom colour of the bar. Field-by-field meanings:
[r-4d, explained](../r-4d-case.md#inspect-remesh-start-fresh) and
[Rarefaction diagnostics](../rarefaction-diagnostics.md).

**Judging the output.** A valid diagnostic cell needs more than 20 sampled
parcels per step; at the shipped weight the exit region is reported invalid.
The cell/λ and Δt/τ warnings are real: near the exit the simulated gas is more
viscous than the real one. Compare a run on the shipped cone mesh with the
recorded full-length run on the earlier cylinder mesh before drawing a
conclusion from either.

## Resources

Recorded by the owner (see the explanation for the mesh each belongs to):
`baseCase` builds about 459 000 cells; the earlier cylinder mesh ran its full
45 348 steps on 8 ranks in 18.3 h and levelled off at 1.23 × 10⁶ parcels; a
written time is 90–280 MB. `--describe` gives the estimate for any edit. The
HPC variants' sizes are estimates only; neither has a recorded run.

## Figure placeholder

!!! figure "Figure placeholder — plate pressure on the upstream face"
    **Subject:** `results/plate-pressure.png` for `baseCase` at the end of sampling.
    **Source run:** a full-length `baseCase` run on the shipped cone mesh (not yet recorded).
    **Plotting command:** `./Allpost` (verified in source).
    **Destination:** `docs/assets/figures/r-4d-plate-pressure.png`
    **Suggested caption:** "Time-averaged wall-normal pressure on the R-4D plate (under-resolved near field; state the mesh and weight)."

## Related

- [r-4d, explained](../r-4d-case.md)
- [Rarefaction diagnostics](../rarefaction-diagnostics.md) · [Slice imagery](../viz-slices.md)
- [Post-processing](../guides/postprocessing.md) · [Restart and recovery](../guides/restart-and-recovery.md)
- API: [`nozzleplume.config`](../reference/api/nozzleplume/config.md), [`workflow`](../reference/api/nozzleplume/workflow.md), [`mesh`](../reference/api/nozzleplume/mesh.md), [`post`](../reference/api/nozzleplume/post.md), [`iteration`](../reference/api/shared/iteration.md)

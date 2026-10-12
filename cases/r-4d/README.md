# R-4D exit-plane plume onto a flat plate

A single DSMC case: the exit plane of an R-4D nozzle is the inlet, the plume
expands along `+x` into vacuum, and a flat plate stands in its way.

```
      y,z ^                                 vacuum (molecules deleted)
          |   +---------------------------------------------+
          |   |                                   ##        |
          |   |          plume                    ## plate  |
   -------+---#====#------------------------------##------->+  x
          |   |    ^ nozzle exit                  ##        |
          |   |      D = 0.377 m                  ##        |
          |   +---------------------------------------------+
              x = 0                             20 m      22 m
```

It is the first case of the general nozzle-plume family,
[`plumetools/nozzleplume`](../../plumetools/nozzleplume/__init__.py), which is
`plumetools/cai2012` with the paper taken out: the same circular exit into the
same kind of vacuum box, but with the exit state given directly, any single gas,
and a plate. `plumetools/cai2012` and `cases/cai2012*` are untouched.

**Status: a minimum viable case.** It meshes, runs, and post-processes end to
end, and the pieces that have a known answer have been checked against it
([below](#what-has-been-checked)). **No result from it has been validated**
against data or another code, and at this exit density it is under-resolved by
necessity ([below](#it-is-under-resolved-and-says-so)).

## Run it

```bash
cd cases/r-4d/baseCase

./Allmesh --describe   # what case.yaml implies and what it costs. Writes nothing.
./Allmesh              # dictionaries, blockMesh, snappyHexMesh, nozzle patch, verify
./Allrun               # 0/ fields, dsmcInitialise, decomposePar, dsmcFoam, reconstructPar
./Allpost              # plate loads and centreline -> results/   (--viz: slice imagery)
./Allviz               # the flow field as 2-D slices -> results/viz/   (needs pvpython)
./Allclean             # reset (--mesh also removes the mesh; --dry-run previews)
```

Needs OpenFOAM v2512 on `PATH`, `libplumeDsmcBoundaryModels.so` built
(`./Allsetup` from the repository root), and `plumetools` importable.
`./Allviz` needs ParaView's `pvpython` instead of OpenFOAM.

`./Allrun` **continues** a case from its newest time directory, so an
interrupted run is restarted by running it again. To run past the end time,
raise it in `case.yaml`, then `./Allmesh --dict-only && ./Allrun`.

## One file

Everything is in [`baseCase/case.yaml`](baseCase/case.yaml). Every OpenFOAM
dictionary, the mesh and the `0/` fields are generated from it and gitignored.
After any edit, `./Allmesh --describe` prints what the file now implies.

### What was given

| | `case.yaml` | value |
|---|---|---|
| exit diameter | `nozzle.exit_diameter_m` | 0.376936 m |
| exit temperature | `exit.temperature_K` | 498.04 K |
| exit bulk speed | `exit.speed_m_per_s` | 3187.19 m/s |
| exit number density | `exit.number_density_per_m3` | 1.306e22 m⁻³ |
| molecular mass | `gas.molecular_mass_kg` | 3.40055e-26 kg |
| molecular diameter | `gas.molecular_diameter_m` | 3.5e-10 m |
| plate | `plate.size_y_m`, `size_z_m`, `thickness_m` | 34.14 × 11.88 × 0.25 m |

The rest of the data sheet is **derived, and printed, not typed in** — so it is
a cross-check rather than a second source of truth:

| derived | value | data sheet |
|---|---|---|
| exit area | 0.111590 m² | 0.111590 m² |
| mass density | 4.4411e-4 kg/m³ | 4.442e-4 kg/m³ |
| collision cross-section | 3.84845e-19 m² | 3.84845e-19 m² |
| Knudsen number (on D) | 3.73e-4 | "10⁻⁴" |
| mean free path | 0.141 mm | |
| speed ratio / Mach | 5.01 / 5.49 | |
| mass flow, thrust, Isp | 0.158 kg/s, 513 N, 331 s | an R-4D is a ~490 N engine |

The exit density may instead be given as `mass_density_kg_per_m3` or `knudsen`,
and the speed as `speed_ratio` — exactly one of each.

### What was assumed

None of these were supplied. Each is one line of `case.yaml`, tagged
`[ASSUMPTION]` there.

| assumption | value | why, and what it changes |
|---|---|---|
| plate orientation | normal to the plume axis, centred on it | "impinge". 34.14 m along `y`, 11.88 m along `z`, 0.25 m along `x` — the footprint of the `panel` patch in `cases/iss_solar_panels_3`. `plate.centre_y_m` / `centre_z_m` move it off-axis; a tilted plate is not expressible yet. |
| stand-off | 20 m, exit plane to upstream face | where that same panel sits. `plate.standoff_m`. |
| wall | diffuse, fully accommodated, 300 K | `dsmc.wall_interaction_model`, `plate.temperature_K`. Sets the re-emitted momentum and the heat flux. |
| collision model | hard sphere (`gas.omega: 0.5`) | it makes the cross-section exactly the πd² on the data sheet at every temperature. A VHS exponent needs a reference temperature the sheet does not give. |
| internal modes | none (`gas.internal_degrees_of_freedom: 0`) | so the **simulated** gas has γ = 5/3. Real exhaust has rotational modes; set the count and `dsmc.collision_model: LarsenBorgnakkeVariableHardSphere` to model them. |
| exit profile | uniform over the disk | one `n`, `U`, `T`. No boundary layer, no divergence. |
| upstream of the exit | open vacuum | the rest of the `x = 0` plane deletes molecules. There is no nozzle body or spacecraft for back-scattered gas to hit. |
| one species | the mixture's mean molecular mass | as given. |

## The mesh

All OpenFOAM-native, from analytic primitives — no STL, no external mesher:

| step | does |
|---|---|
| `blockMesh` | a uniform **cubic** background box, 0.5 m cells |
| `snappyHexMesh` | octree-refines the plume volumes; carves the plate as a `wall` |
| `topoSet` + `createPatch` | cut the `nozzle` disk out of the `x = 0` plane |
| `checkMesh`, then `verify` | quality; then that the mesh **is** what `case.yaml` describes |

Cell size at refinement level *N* is `background_cell_size_m / 2**N`:

| where | level | cell | |
|---|---|---|---|
| the plume cone, to 1 m | 4 | 31.25 mm | 12 cells across the exit |
| to 3 m | 3 | 62.5 mm | |
| to 8 m | 2 | 0.125 m | |
| to the plate | 1 | 0.25 m | |
| plate surface | 2 | 0.125 m | 2 cells through its thickness |
| elsewhere | 0 | 0.5 m | |

That is **459,007 cells** (`Mesh OK`, max non-orthogonality 28°), built in about
a minute and a half.

### The refinement regions are cones

`mesh.refinement_regions` is a list of **cones** coaxial with the plume, each
written to `snappyHexMeshDict` as a `searchableCone`:

```yaml
- {name: plumeCore, level: 2, x_max_m: 8.0, radius1_m: 0.4, radius2_m: 2.0}
```

`radius1_m` is the radius at `x_min_m` (default 0, the exit plane) and
`radius2_m` the radius at `x_max_m` — the `radius1` and `radius2` of the
dictionary — with a straight line between. Equal radii are a cylinder.
`inner_radius1_m` and `inner_radius2_m` hollow it out; left out, it is solid.
(`shape: box`, with `y`/`z` bounds, is the other shape.) `--describe` prints
each one with the half-angle its radii imply.

The four shipped lie on **one surface, r = 0.4 m + 0.2 x** — a half-angle of
11.3° — cut off at 1, 3, 8 and 20 m. A plume widens as it goes, and a cone
follows it; a cylinder wide enough for the far end of a stage is far wider than
the plume at the near end. This case was first meshed with four cylinders
reaching the same radii at those stations, and that took 529,126 cells: the
cones give up 70,000 of them, none on the exit or on the plate, whose patches
have the same faces either way. In the run made on the cylinders, the gas on
the cone's surface is at 0.3% of the on-axis density at 1 m, 4% at 3 m and 16%
at 8 m, with a mean free path longer in proportion — so the cells given up are
the ones that needed refining least.

The slope is a choice, not a given. `1/S` is where to start for another exit:
0.2 is this one's, at a speed ratio of 5.0.

Edit the list freely — `--describe` re-estimates the cell count, to within 2%
of what snappyHexMesh builds for this case. Three things worth knowing first:

* **The exit must sit at one refinement level**, and snappyHexMesh refines a
  cell when its *centre* is inside a region — so a region smaller than the cells
  around it is never seen at all. For a cone that is decided where the centres
  are, not on the exit plane: the first background cells have theirs at
  x = 0.25 m, 0.354 m off the axis, and the exit cone has to have widened past
  that by then. The plan applies the rule and refuses a configuration that
  would split the exit, rather than letting it mesh.
* **The nozzle is a staircase.** Its area is 1.98% short of the disk at 12 cells
  across, so the injected mass flow is 1.98% low. That is not monotone in
  resolution: it is +0.64% at level 5 and −0.02% at level 6. A finer exit also
  halves the time step per level.
* **`n_cells_between_levels` is not a number of cells**, and it moves the cell
  count in steps. snappyHexMesh builds buffers out of whole cells of the level
  below, so the skin of cells kept at the plate's level is 2 deep at 1 or 2,
  4 deep at 3 or 4, and 6 at 5 — and on 811 m² of plate that skin is half the
  mesh:

  | `n_cells_between_levels` | 1 | 2 | **3** | 4 | 5 |
  |---|---|---|---|---|---|
  | cells built | 318,104 | 341,512 | **459,007** | 489,905 | 640,916 |

The plate is exact: every face of its patch is axis-aligned and its upstream
face is 405.5832 m² = 34.14 × 11.88.

## It is under-resolved, and says so

DSMC wants cells smaller than a mean free path and steps shorter than a mean
collision time. At an exit Knudsen number of 4e-4 that is not available: the
0.14 mm mean free path would need 2,700 cells across the exit.

| at the exit | this case | criterion |
|---|---|---|
| cell / mean free path | 222 | ≤ 1 |
| time step / mean collision time | 7.0 | ≤ 1 |

`./Allmesh` prints both as **warnings**, with what each costs, every time. They
ease downstream as the plume thins — roughly with the density — but near the
exit the simulated gas is more viscous than the real one and gradients finer
than a cell are smeared. Read the near field with that in mind; this is why
"expensive" is the first thing anyone says about this case.

## Cost, and the knobs

`dsmc.n_equivalent_particles` — real molecules per simulated parcel — is written
to `dsmcProperties` exactly as given. `dsmcFoam` has one weight for the whole
domain, so it is the trade between noise and run time:

| `n_equivalent_particles` | parcels per exit cell | parcels per step | parcels in the domain (rough) |
|---|---|---|---|
| 2e16 | 20 | 312 | 3.9e6 |
| **5e16** (shipped) | 8.0 | 125 | 1.5e6 |
| 1e17 | 4.0 | 62 | 7.7e5 |
| 2e17 | 2.0 | 31 | 3.9e5 |

Run `./Allmesh --describe` after changing it; nothing needs re-meshing, but a
run in progress cannot change weight — `./Allclean` first.

As shipped the run is **45,348 steps** of 1.37 µs (6 domain transits of
transient, then 3 of averaging).

The timings below were all taken on the mesh this case had first — cylinders
where the cones are now, 529,126 cells. The parcels are what a step costs, and
they are the same on either mesh; the cells are 13% fewer now.

Measured over the first 10,000 of those steps on a 12-core workstation, 12
ranks: the wall-clock time per step is **0.05 s + 0.57 µs per parcel**, plus
about 0.12 s once `fieldAverage` is accumulating. Those 10,000 steps took 74
minutes end to end and finished at 1.07 million parcels and 0.85 s a step.

The parcel count was levelling off by then — net growth had fallen from 62 to 29
thousand per 500 steps — which pointed at a plateau near 1.2–1.3 million,
against the table's rough 1.5 million, and at a whole run of roughly eleven
hours. The whole run has since been made, on 8 ranks: it levelled off at **1.23
million parcels** and took **18.3 hours**. Watch `Number of dsmc particles` in
`log.dsmcFoam` flatten rather than trusting an estimate of how long gas is held
up at a plate. A written time is 90–280 MB, growing with the parcel count.

The other levers, in rough order of value:

* **`domain.symmetry_planes: [y, z]`** models one quarter — a quarter of the
  cells *and* of the parcels for the same answer, as long as the plate stays
  centred. The totals `./Allpost` reports are still for the whole plate.
* **Run length**: `dsmc.transient_domain_transits`, `sampling_domain_transits`,
  or pin `end_time_s` / `average_start_s`.
* **The exit level** only moves one way. Each level finer halves the time step
  and divides the parcels per exit cell by eight. One level coarser is not on
  offer: a 6-cell nozzle is 12% out in area, and the checks refuse it.
* `dsmc.n_subdomains`, and `decomposition_method: simple` if `scotch` leaves one
  rank with most of the parcels (compare the ranks in `log.dsmcFoam`).

## What comes out

`./Allpost` writes `results/`:

| file | |
|---|---|
| `plate_loads.csv` | pressure, shear and heat flux on every plate face |
| `plate_summary.yaml` | force and heat load on the plate, beside the jet's thrust and power |
| `plate-pressure.png`, `plate-heat-flux.png` | maps of the upstream face |
| `centreline.csv`, `centreline.png` | `n`, `U_x`, `T` along the axis, over the exit values |
| `viz/` (`./Allviz`, or `./Allpost --viz`) | slices of number and mass density, parcel count, speed and temperature through the plume |

Plate pressure is `fD · n` from the time-averaged `fDMean` — the momentum the
gas actually delivers — not `nkT` from the adjacent cell, which would miss the
directed momentum that is nearly all of the load.

The centreline is the mean over the cells beside the axis, and it **steps down
wherever the mesh coarsens**. That is the sample getting wider across a peaked
profile, not the flow changing: at `x = 3 m` the four coarse cells read exactly
what the sixteen fine ones under them average to, and the mass flux through the
interface is continuous to 0.7%. `centreline.csv` carries the cell size each
sample came from and the figure marks every change.

### The flow field

`./Allviz` hands the case to the repository's ParaView slice renderer
([`docs/viz-slices.md`](../../docs/viz-slices.md)) with
[`baseCase/viz.yaml`](baseCase/viz.yaml) — the renderer, and the time-averaged
samples, that `cases/cai2012-health` draws its mid-plane from. It writes
`results/viz/`:

| image | field | colour scale |
|---|---|---|
| `rhoNMean_<plane>.png` | number density [1/m³] | log, 8 decades below the peak |
| `rhoMMean_<plane>.png` | mass density [kg/m³] | log, 8 decades below the peak |
| `dsmcRhoNMean_<plane>.png` | simulated parcels per cell | log, 6 decades |
| `U_<plane>.png` | speed [m/s] | 0 – 3500, pinned |
| `Ttra_<plane>.png` | translational temperature [K] | 0 – 6500, pinned |
| `fDMean_plate.png`, `qMean_plate.png` | force density [Pa] and heat flux [W/m²] on the plate's loaded face | follows the data |

`<plane>` is `xz` and `xy`, the two centre planes through the plume axis, and
`plate`, a cross-section 0.2 m ahead of the plate seen from the nozzle side: 17
images, a little over a minute on a mesh of this size. `manifest.yaml` beside them records
the time and the colour range each was drawn at.

```bash
./Allviz                          # everything, at the latest time
./Allviz --field Ttra --plane xz  # one picture
./Allviz --time 0.0517723         # another written time
./Allviz --dry-run                # name the images; draw none
```

- **Look at the parcel count first.** It is the map of where the other four are
  shot noise: one particle weight serves the whole domain, so in the near-empty
  region beside the nozzle a cell averages far less than one parcel and its
  speed and temperature mean nothing.
- **The speed and temperature scales are pinned** just above what the gas can
  reach from the exit state — the limiting speed, 3342 m/s, and the stagnation
  temperature, 5502 K. Autoscaled, they follow those empty cells to 6600 m/s and
  8700 K. The numbers are in `viz.yaml`, not derived from `case.yaml`: change
  the exit state and move them (a unit test fails until you do).
- **Mass density is number density times the molecular mass**, exactly, for a
  single species: the same picture read in kg/m³.
- The averaged fields exist only at times after sampling starts
  (`./Allmesh --describe`). An earlier time is drawn from one step's parcels,
  and the renderer says so.

`./Allpost --viz` draws the same images after the plate loads and the
centreline. Both need `pvpython` on `PATH` (or `PVPYTHON_PATH`). `open.foam`
opens the case in ParaView directly.

## What has been checked

Against OpenFOAM v2512, with the case's own scripts. The case has had two
meshes, and each check below belongs to one of them.

**On the mesh as it ships**, with the cones:

* **The mesh is the one described.** `./Allmesh` ends by measuring it: 112
  nozzle faces as predicted, plate box and face area exact, patch types right.
* **snappyHexMesh refines exactly where the cones say.** The finest level holds
  nothing but the exit cone, so its cells can be counted instead of estimated:
  eight for every 62.5 mm cell whose centre is inside the cone, 25,824 — and
  25,824 were built, at each `n_cells_between_levels` from 1 to 5. A hollow
  cone added at a level of its own came out the same way: 78,016 predicted,
  78,016 built.
* **The solver takes it.** The first 1,500 steps on 12 ranks, at the shipped
  weight: 124.74 ± 0.16 parcels injected per step against 124.744 predicted
  from the Maxwellian flux, through the `nozzle` patch only, and no solver
  warning. At 187,000 parcels a step was taking 0.14 s, which is what the
  timings [above](#cost-and-the-knobs), taken on the earlier mesh, give for
  that many.
* **The plate feels what conservation says it must.** With collisions off and a
  plate that catches the whole jet, momentum and energy balance alone fix the
  totals. A quarter-symmetry run returned the force to 0.1% and the heat load to
  0.2% (F/thrust 1.1105 against 1.1120, Q/power 0.9352 against 0.9374). That is
  inflow, wall model, `fD`/`q` accumulation, symmetry and post-processing
  agreeing at once, and it is kept as `tests/openfoam/test_nozzleplume_openfoam.py`.

**On the mesh it had first** — four cylinders where the cones are, 529,126
cells, with the same exit, plate, particle weight and run length:

* **Injection is right.** 124.745 ± 0.08 parcels per step measured over 5,770
  steps, against 124.744 predicted from the Maxwellian flux — through the
  `nozzle` patch only.
* **It survives impingement, and behaves like a gas.** That mesh and the
  shipped weight were run for 10,000 steps (13.7 ms, two domain transits) on
  12 ranks; the jet reaches the plate at about step 4,600. No solver warning,
  1.07 million parcels at the end. Three things in that transient are worth
  more than "it ran", because nothing imposed them:
  * the lateral force on the plate is −0.14 N and +0.04 N beside 529 N axial —
    a full 3-D mesh, symmetric to 0.03%;
  * the plume accelerates to 1.0483 of its exit speed and stays there, against
    the 1.0486 a γ = 5/3 gas reaches when all its thermal energy is spent;
  * a shock layer stands about 1.7 m ahead of the plate, and the translational
    temperature in it peaks at 10.8 times the exit value, against a stagnation
    temperature of 11.0 times.
* **It runs to its end.** All 45,348 steps, on 8 ranks: the parcel count is
  steady at 1.23 million from about step 25,000. Averaged over the last three
  domain transits the plate takes **553.4 N**, 1.078 times the jet's thrust,
  and 0.281 of its power as heat; the pressure on it peaks at 9.10 Pa; and the
  lateral force is 0.04 N and 0.02 N — symmetric to 0.01%.

What has **not** been done: the mesh as it ships has not been run to its end
time. **The converged plate load above is the cylinder mesh's**, and it is what
a run on the cones is to be compared with — the two differ only in where the
plume is refined, so a difference between them is a measure of how much that
matters. Nothing has been compared with measurement, with another DSMC code, or
with a continuum plume model; and no resolution or weight study has been made.

## Tests

```bash
pytest tests/unit -k nozzleplume                    # 383 tests, no OpenFOAM needed
pytest tests/openfoam/test_nozzleplume_openfoam.py -m needs_openfoam    # ~1 min
```

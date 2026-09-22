# Cai & Wang 2012 — collisionless circular plume

Cai, C., and Wang, L., *Numerical Validations for a Set of Collisionless Rocket
Plume Solutions*, Journal of Spacecraft and Rockets, Vol. 49, No. 1, 2012.

Argon expands from a **circular nozzle exit directly into vacuum**: a uniform
drifting Maxwellian on the disk `r ≤ D/2` at `x = 0`, with `D = 0.2 m` and speed
ratio `S0 = 2`. Three Knudsen numbers — 100, 0.1, 0.01 — compared against Cai's
collisionless analytical solution.

Each of those three physical cases is run at **four numerical-particle
resolutions** — 1×, 2×, 5× and 10× the baseline parcel population — so the
study is 3 × 4 = 12 cases. That axis varies `nEquivalentParticles` and nothing
else; see [Two axes](#two-axes-physical-and-numerical) below.

The physical nozzle exit **is** the DSMC inlet. There is no hemispherical source
surface and no analytical inflow model anywhere in this case family; that is the
difference from `cases/3d-inflow` and `cases/markelov1999`.

> **This is not a reproduction of Cai's DSMC results.** Cai used an axisymmetric
> DSMC code (GRASP) on a uniform grid; this is 3-D Cartesian OpenFOAM
> `dsmcFoam` v2512 at a different resolution. What is being validated is the
> *trend* across the three Knudsen numbers. See
> [`docs/cai2012-case.md`](../../docs/cai2012-case.md) §8.

## Quick start

```bash
./Allsetup                       # from the repository root: pip install, build the
                                 # custom inflow model, verify both
cd cases/cai2012

./generate_cases.py --dry-run    # what would be generated
./generate_cases.py              # 3 Knudsen cases x 4 particle levels = 12
./AllmeshCases                   # blockMesh + topoSet + createPatch
./validate_cases.py              # check the 12 before spending HPC hours on them
./AllrunCases                    # solve
./AllpostCases                   # the four validation outputs, plus the study table
```

Subsets, for a partial run or a job array:

```bash
./generate_cases.py --knudsen 100        # one physical case, all four levels
./generate_cases.py --particles 1        # all three cases, the 1x control only
./AllrunCases --list                     # the 12 cases, numbered
./AllrunCases --index 0                  # one of them, by that number
./AllrunCases --case Kn100_np2x          # or by name
```

Per case, from inside `Cases/Kn100_np1x`:

```bash
./Allmesh              # dictionaries + mesh + nozzle carving + verification
./Allmesh --dict-only  # dictionaries only; needs no OpenFOAM
./Allmesh --force      # re-mesh even though the case has results
./Allrun --serial      # one process
./Allrun -np 12        # override the subdomain count
./Allrun --no-solve    # generate 0/ and run the checks; needs no solver
./Allpost              # centreline, contour, error metrics, figures
./Allclean             # remove solver output (--mesh to remove the mesh too)
```

## Runs continue; they do not restart

Re-running `./Allrun` **continues** a case from its newest time directory.
`system/controlDict` carries `startFrom latestTime`, and the two steps that would
destroy the state being resumed from are skipped:

| Step | On a resume |
|---|---|
| `runInflow.py` | skipped — the solver reads the boundary fields from the resume time, not from `0/` |
| `dsmcInitialise` | skipped — it would rebuild the initial parcel cloud |
| `decomposePar -force` | skipped if `processor*/` already hold the latest time; otherwise narrowed to `-latestTime` |
| `reconstructPar` | `-newTimes`, so earlier writes are not reconstructed again |

The **averaging continues too**: `fieldAverage` stores its accumulators in each
time directory's `uniform/functionObjects/`, so the mean fields span the whole
sampled period rather than only the last leg.

So an interrupted run is resumed by running `./Allrun` again, and a finished run
is carried further by raising `dsmc.end_time_s` (or
`dsmc.sampling_domain_transits`) in `case.yaml` and then:

```bash
./Allmesh --dict-only    # rewrite controlDict; does NOT re-mesh
./Allrun                 # continue from the latest time
```

Two guards come with this:

* **`./Allmesh` refuses to re-mesh a case that has results.** blockMesh would
  replace the geometry those results were computed on while leaving them in
  place, addressing cells that no longer exist — and nothing downstream would
  notice. `--dict-only` is always allowed; `--force` overrides.
* **On a resume `0/` is not regenerated.** `boundaryT`, `boundaryU` and
  `boundaryNumberDensity_Ar` are `AUTO_WRITE`, so every time directory carries
  its own copy and the solver reads the resume time's. Rewriting `0/` would
  change nothing while looking as though a `case.yaml` edit had taken hold. To
  change the physics, start over with `./Allclean`.

`./Allclean` is the only thing here that deletes a time directory, and
`./generate_cases.py` already refuses to overwrite a case that has results unless
`--overwrite` is given.

## Layout

```
cases/cai2012/
  study.yaml            the case matrix: three Knudsen numbers x four particle levels
  generate_cases.py     clone baseCase per pair, derive everything, write manifest
  validate_cases.py     check the generated tree against study.yaml
  baseCase/             the template — case.yaml plus the scripts
    case.yaml           THE authoritative source of physical inputs
  Cases/                generated; gitignored
    Kn100_np1x/ Kn100_np2x/ Kn100_np5x/ Kn100_np10x/
    Kn0p1_np1x/  ...                       12 in total
      results/          centerline.csv, density_plane.csv, metrics.yaml, *.png
  manifest.yaml         generated; the auditable case → inputs map
  results/              generated; study-table.csv
```

`Cases/`, `manifest.yaml` and `results/` are gitignored: they are reproducible
from `baseCase` and `study.yaml`, and committing a generated tree would let the
template and its copies drift apart.

## Two axes: physical, and numerical

```
kn_cases         Kn = 100, 0.1, 0.01                [PAPER]
particle_levels  np1x, np2x, np5x, np10x            numerical resolution
```

The product is generated, one directory per pair. A level changes exactly one
thing:

```
nEquivalentParticles = baseline / multiplier
```

because one DSMC parcel stands for `nEquivalentParticles` real molecules, so the
numerical parcel population goes as the **reciprocal** of the weight and 10× the
parcels means a tenth of the weight.

**The baseline is never written down.** It is whatever
`plumetools.cai2012.inflow.derive_run_settings` derives for that case from its
own `n0` and exit cell volume — the same number this study has always used — so
`np1x` *is* the case as it was, and the sweep cannot drift away from the study it
extends. `study.yaml` is refused if it has no 1× level, because without the
control the other three are unreadable.

| | Kn100 | Kn0p1 | Kn0p01 |
|---|---|---|---|
| `np1x` weight | 3.310039e+09 | 3.310039e+12 | 5.066386e+12 |
| `np2x` | 1.655019e+09 | 1.655019e+12 | 2.533193e+12 |
| `np5x` | 6.620078e+08 | 6.620078e+11 | 1.013277e+12 |
| `np10x` | 3.310039e+08 | 3.310039e+11 | 5.066386e+11 |
| parcels at 10× | 1.2e+07 | 1.2e+07 | 7.5e+07 |

Everything else is **identical across a case's four variants** — mesh, geometry,
gas, species, boundary conditions, collision model, density, `deltaT`, and the
transient. That is not a convention, it is checked: `./validate_cases.py`
compares the generated dictionaries and fails if any of it moved. `study.yaml`
may not set `dsmc.numerical_particle_multiplier` through a per-case `overrides`
block, the same way it may not set `exit.knudsen`.

### Longer runs, more frequent writes

Two study-wide scales, in `study.yaml`:

```yaml
run_time_multiplier: 3.0
output_frequency_multiplier: 4.0
```

`run_time_multiplier` scales `endTime` and **not** the averaging start. The
transient is a measured statement about when the plume is established — the
parcel count is flat after 1.32 domain transits — so it stays where it is and
every second the multiplier adds lands inside the averaging window. The run is
3× longer; the *sampled* window is 5.7× longer.

`output_frequency_multiplier` divides the write interval. `writeControl` stays
`timeStep` — see below for why it is not `runTime` — so the interval is a whole
number of steps and 1443 → 361 is a 3.997× increase rather than exactly 4×. The
manifest records both the request and what was achieved.

| | before | after |
|---|---|---|
| `endTime` (Kn100) | 9.90534e-03 s | 2.97366e-02 s |
| `deltaT` | 1.37288e-06 s | **unchanged** |
| `writeControl` | `timeStep` | **unchanged** |
| `writeInterval` | 1443 steps | 361 steps |
| writes | 5 | 60 |
| writes inside the averaging window | 3 | 49 |

That is roughly 12× the output of the original study per case, and 18× the
parcel-seconds across the sweep. Both are intended: these run on HPC, and the
fields are for a statistical-convergence comparison that needs frames.

### Variance, not just means

`fieldAverage` now keeps `prime2Mean` for the two number densities:

```
rhoNMean       rhoNPrime2Mean
dsmcRhoNMean   dsmcRhoNPrime2Mean
```

so the statistical scatter in a cell can be **measured** off the solution,

```
CV = sqrt(dsmcRhoNPrime2Mean) / dsmcRhoNMean
```

rather than estimated from an assumed parcel count — which is the quantity a
numerical-particle convergence study is about. `dsmcRhoN` is the parcel tally, so
its CV is the raw sampling noise; `rhoN` is the real number density, so its CV is
the error on the quantity actually compared with Cai. `prime2Mean` stays **off**
for the vector fields: their variance is a `symmTensor` per cell per write,
nothing reads it, and this study writes 60 times per case.

## Only the Knudsen number is a per-case input

`study.yaml` sets `Kn` and nothing else physical. Everything the solver runs at
is derived:

```
Kn, D          →  λ0 = Kn·D
λ0, T0         →  n0
S0, T0         →  U0 = S0·√(2RT0)
n0, cell       →  particle weight
cell, U0       →  deltaT
```

So a case named `Kn0p1` cannot run another case's density: no density is stored
anywhere. `study.yaml` is **forbidden** from overriding `exit.knudsen` — that is
rejected with an error, not merged.

`manifest.yaml` records `Kn`, `D`, `λ0`, `T0`, `U0`, `n0`, the particle weight,
the minimum cell size and `deltaT` for every case, so a results directory can
always be traced back to its inputs. It also records, per case, the physical
case it belongs to (`cai_case`), its particle level and multiplier, the baseline
weight the level was derived from, a `mesh_id` that is equal across a group,
`startTime` / `endTime` / `writeControl` / `writeInterval` and the multipliers
that produced them, and the whole `fieldAverage` configuration including which
fields keep a `prime2Mean`.

## `case.yaml` is the only place physics lives

`constant/dsmcProperties`, `system/controlDict`, `system/dsmcInitialiseDict`,
`system/decomposeParDict`, the three mesh dictionaries and every `0/` field are
**generated** from `case.yaml`. They are gitignored and overwritten on every
`./Allmesh`. Do not edit them — edit `case.yaml` and regenerate.

## What the mesh does, and why it is not Cai's

Cai's grid is uniform at `Δx = λ0(Kn = 0.01) = 2 mm`. In this 3-D box that is

```
1000 × 2000 × 2000 = 4.0 × 10⁹ cells
```

which `./Allmesh` computes and prints **before it writes anything**. Cai can
afford his grid because his solver is axisymmetric — the mesh is a 2-D `(x, r)`
sheet.

So this builds a graded multi-block box instead: a uniform fine core around the
nozzle and the near plume, with geometric expansion outward. The expansion ratios
are solved for so the first outer cell matches the core cell; `mesh.max_cells` is
a hard budget, and a mesh that had to be coarsened to fit says so on every report
and in the manifest.

| Case | `λ0` | core cell | `Δx/λ0` | cells | parcels | steps |
|---|---|---|---|---|---|---|
| Kn = 100 | 20 m | 10 mm | 0.0005 | 1.3 M | ~1.2 M | 7 212 |
| Kn = 0.1 | 20 mm | 10 mm | 0.5 | 1.3 M | ~1.2 M | 7 212 |
| Kn = 0.01 | 2 mm | 5.4 mm | **2.69** | 3.7 M | ~7.5 M | 13 523 |

**Kn = 0.01 does not meet Cai's `Δx/λ0 = 1`.** The checks warn about it every
time; see [`docs/cai2012-case.md`](../../docs/cai2012-case.md) §4.3 for what it
costs and how to raise it.

Rough cost on 12 physical cores, at the transient lengths `study.yaml` sets:
Kn = 100 is about an hour, Kn = 0.1 somewhat more, and Kn = 0.01 several hours
with a few gigabytes for its parcels. `dsmc.n_subdomains` is 12 rather than 24
because Open MPI counts cores, not hardware threads, when deciding how many
slots exist.

## The nozzle patch

`blockMesh` cannot make a disk — it can only name block faces — so `x = 0` comes
out as one patch, `upstreamVacuum`, and the nozzle is carved out of it:

```
topoSet       faces of upstreamVacuum with r ≤ R0   →  faceSet nozzleFaces
createPatch   nozzleFaces                            →  patch nozzle
```

The result is a **staircase** approximation to the circle. `./Allmesh` predicts
the face count and area in Python first, then checks the built mesh against the
prediction and fails if they disagree — the two are supposed to be the same rule
(`topoSet`'s `cylinderToFace` and "face centre inside `r ≤ R0`"), and if they are
not, the mesh that was checked is not the mesh that was built.

| Case | faces | area vs `πR0²` |
|---|---|---|
| Kn = 100, Kn = 0.1 | 316 | +0.59% |
| Kn = 0.01 | 1092 | −0.24% |

The injected flow is proportional to the area, so this is a direct multiplier on
the whole solution. Errors above `checks.max_nozzle_area_error` (5%) are fatal.

## Boundaries

| Patch | Type | Behaviour |
|---|---|---|
| `nozzle` | `patch` | injects |
| `upstreamVacuum` | `patch` | the rest of `x = 0`; **deletes** |
| `vacuum` | `patch` | the other five faces; **deletes** |

No boundary is a `wall`, and the checks fail if one ever becomes one: a wall
reflects, and the plume would expand into a closed box.

This is only possible because the case runs **`plumeFieldInflow`**, which takes
an explicit patch list. Stock `FreeStream` injects on every `patch`-type boundary
with no selection, so all six would become inlets — and the only way to stop that
is to make them walls. The case **fails** if the library is not built; there is
no fallback. See [`docs/plume-field-inflow.md`](../../docs/plume-field-inflow.md).

## What `./Allpost` produces

Per case, in `Cases/<name>/results/`:

| File | Contents |
|---|---|
| `centerline.csv` | `n/n0`, `U₁√β0`, `T/T0` vs `X/D`, DSMC and analytical |
| `density_plane.csv` | `n/n0` at the cell centres of the `X–Z` centre plane |
| `metrics.yaml` | max / mean / RMS relative error against Cai |
| `centerline-*.png` | the three centreline figures |
| `density-contour.png` | `n/n0` = 0.1, 0.01, 0.001, DSMC over analytical |

and `results/study-table.csv` across the three cases.

Everything is read from the **time-averaged** fields. A single-timestep `rhoN` is
one step's worth of parcels and comparing it with an analytical curve measures
shot noise.

## Results so far

`Kn = 100` and `Kn = 0.1` have been run to completion; `Kn = 0.01` is generated
and meshed but has not been run (several hours on 12 cores).

| Case | density max | density mean | velocity max | temperature max | Cai's max density |
|---|---|---|---|---|---|
| Kn = 100 | 5.6% | **1.1%** | 1.5% | 4.2% | 0.07% |
| Kn = 0.1 | 17.6% | **6.8%** | 3.5% | 54.9% | 3.63% |

The expected trend: the departure from the collisionless solution grows by 6.4x
in the mean density error and 22x in the temperature. Over `X/D ≤ 3` — where the
shortened transient has no effect — the `Kn = 100` density agrees to about 1%,
and its contour reproduces all three of Cai's levels within the sampling noise.

At `Kn = 0.1` the departure has the right sign everywhere: density above the
collisionless solution near the exit and below it downstream, velocity 2-3%
higher, temperature up to 55% lower. Collisions collimate the near plume, convert
thermal energy into directed motion, and keep cooling where free-molecular flow
freezes.

The `Kn = 100` far field (`X/D > 5`) is **transient-limited**, not converged: the
slow tail of the exit distribution has not arrived by the time sampling starts.
See [`docs/cai2012-case.md`](../../docs/cai2012-case.md) §5.1 — the arithmetic
predicts the measured deficit, and `dsmc.transient_basis: cai` fixes it.

## Status: what this does and does not establish

**Does**: the case matrix generates from one study file; every case runs standard
OpenFOAM v2512 `dsmcFoam`; the physical nozzle is the DSMC inlet; only the nozzle
injects; the outer boundaries are vacuum; `Kn`, `S0`, `D`, `n0`, `U0` and `T0` are
traceable from YAML to the solver.

**Does not**: reproduce Cai's DSMC numbers. Cai's reported maximum centreline
density errors — 0.07% / 3.63% / 12.03% — are carried in the code as **reference
values, not tolerances**. Different solver, different mesh, different resolution.
The Kn = 0.01 case in particular is under-resolved against Cai's own criterion,
which biases its collision rate low.

No claim of agreement is made anywhere in this directory beyond what
`results/metrics.yaml` measures.

## Differences from `cases/markelov1999`

That family is the design this one follows: `study.yaml` → generated cases →
`case.yaml` per case → generated dictionaries → manifest. Two deliberate
departures:

* **No CaseFoam.** The Markelov study is a *hierarchy* (gap over pressure) and
  CaseFoam builds hierarchies. This one is a flat list of three cases, so
  `generate_cases.py` uses `shutil.copytree` and has no extra dependency.
* **Its own config loader.** `plumetools.config` cross-validates a source-flow
  sphere radius and a throat radius that this case does not have, and the frozen
  regression golden is bound to it. `plumetools.cai2012.config` is separate, and
  the two loaders reject each other's files by `model:`.

## Related documents

* [`docs/cai2012-case.md`](../../docs/cai2012-case.md) — paper values,
  assumptions, and every difference from Cai's setup
* [`docs/plume-field-inflow.md`](../../docs/plume-field-inflow.md) — the inflow model
* [`docs/solver-compatibility.md`](../../docs/solver-compatibility.md) — standard
  `dsmcFoam` versus the MNF fork

# Rarefaction and mesh-resolution diagnostics

After a run, how fine does the **next** mesh need to be, and where?
[`plumetools/rarefaction`](../plumetools/rarefaction/__init__.py) answers that
from the time-averaged fields the run sampled: the local mean free path, several
Knudsen numbers, and how many mean free paths each cell spans. It writes them
as volume fields you can open in ParaView, and summarises them in a YAML file.

It does **not** refine anything. There is no adaptive mesh refinement here:
`dsmcFoam` v2512 builds a static mesh, and the decision about what to change is
made by a person editing `case.yaml`. The workflow is:

| step | command | |
|---|---|---|
| 1. run | `./Allrun` | samples accumulate in `<time>/` |
| 2. inspect | `./Allpost`, then ParaView on `open.foam` | writes the diagnostics |
| 3. plan | read `results/diagnostics_summary.yaml` | where, and how fine |
| 4. new iteration | `python -m plumetools.iteration <case> <new>` | a fresh case dir |
| 5. edit, remesh | edit `<new>/case.yaml`; `./Allmesh` | refinement by hand |
| 6. start fresh | `./Allrun` | new parcels from t = 0 |

Nothing is mapped from one iteration to the next: not parcels, not fields, not
`fieldAverage`'s accumulators.

## What is computed

Per cell, from the averaged fields of **one** sampling window:

| field | definition | units | what it measures |
|---|---|---|---|
| `mfp` | λ, equation (G1) below | m | |
| `KnD` | λ / D, D the nozzle exit diameter | – | rarefaction, nozzle scale |
| `Kn<name>` | λ / L for each `diagnostics.reference_lengths_m` entry | – | rarefaction, that scale |
| `cellSizeV` | Δ_V = V^(1/3) | m | |
| `cellSizeMax` | Δ_max, the largest bounding-box dimension | m | |
| `KnCellV`, `KnCellMax` | λ / Δ_V, λ / Δ_max | – | **mesh resolution** |
| `cellOverMfpV`, `cellOverMfpMax` | Δ_V / λ, Δ_max / λ | – | **mesh resolution** |
| `KnGradN` | λ \|∇n\| / n | – | rarefaction, density gradient scale |
| `KnGradT` | λ \|∇T\| / T | – | rarefaction, temperature gradient scale |
| `cellSizeTarget` | λ / k, k = `cells_per_mean_free_path` = 3 | m | the cell that resolves λ |
| `mfpResolved` | 1 where λ/Δ_max ≥ k, 0 where not, −1 no value | – | **the mesh criterion**, as a map |
| `parcelsPerStep` | `dsmcRhoNMean`: parcels in the cell, per step | parcels | DSMC occupancy |
| `sampledParcels` | parcels per step × steps averaged | parcel samples | statistics behind the means |
| `diagStatus` | 0 valid, else a sum of flags (below) | – | |
| `diagValid` | 1 where the cell's values are valid | – | |
| `diagValidGradN`, `diagValidGradT` | 1 where the gradient could be fitted | – | |

### Two different questions

The **cell** numbers — `cellOverMfp*`, `KnCell*` — are about the *mesh*: how
many mean free paths one cell spans. The usual DSMC guidance is Δ ≲ λ/3, and a
cell above that is reported. The **reference-length** and **gradient** numbers
— `KnD`, `KnGradN`, `KnGradT` — are about the *flow*: how rarefied it is,
measured against the nozzle and against the lengths over which its own
density and temperature change. A flow can be near-continuum on the nozzle
scale (`KnD` ≈ 4e-4 at the R-4D exit) and strongly rarefied on a gradient
scale (in a shock layer), and a mesh can resolve one region's λ and miss
another's.

**No one of these numbers, and no threshold on any of them, establishes that a
DSMC solution is valid.** The time step against the mean collision time, the
parcels per cell that collision partners are drawn from, the sampling, and
the length of the transient all matter as well — and some of those are
reported by `./Allmesh --describe`, not here.

## Inputs, and what they really are

| input | read from | |
|---|---|---|
| n | `rhoNMean` | molecules per m³ |
| ρ, ρ**u**, ½ρ⟨c²⟩ | `rhoMMean`, `momentumMean`, `linearKEMean` | the moments T is built from |
| parcels per step | `dsmcRhoNMean` | a **count**, not a density |
| steps averaged | `<time>/uniform/functionObjects/functionObjectProperties` | `fieldAverage`'s `totalIter` |
| V, cell centres | `V`, `C` | `postProcess -func writeCellVolumes` / `writeCellCentres` |
| Δ_max, face neighbours | `constant/polyMesh` | |

The field names are `diagnostics.field_names` in the settings, so a run with
other names needs no code change.

**`dsmcRhoN` is a parcel count.** `DSMCCloud::calculateFields` adds one per
parcel and never divides by the cell volume, while `rhoN` is scaled by
`nParticle/V` — yet both files declare `[0 -3 0 0 0 0 0]`. Every run checks the
identity that follows, `rhoNMean · V / dsmcRhoNMean = nEquivalentParticles` in
every occupied cell, and the summary's `checks.particle_weight` says whether it
holds. If it does not, a field mapping is wrong.

**One sampling window.** `fieldAverage` records `totalIter` and `totalTime` for
every field it averages. If the fields were averaged over different windows —
a field added to `fieldAverage` on a restart, say — they are refused, because
a temperature from moments of different samples is not a temperature.

**Temperature from moments.** T is the translational temperature of the
averaged moments,

    3 R T = 2 ⟨½ρc²⟩/⟨ρ⟩ − |⟨ρu⟩/⟨ρ⟩|²,      R = k_B/m,

never an average of instantaneous temperatures. In a thinly sampled cell it can
come out at or below zero; that cell is flagged, not clamped.

## The mean free path

    λ = 1 / (√2 π d_ref² n) · (T / T_ref)^(ω − 1/2)        (G1)

which is [`plumetools.nozzleplume.gas`](../plumetools/nozzleplume/gas.py)'s
equation, applied element by element and tested against it. The species and
`gas.mean_free_path_convention` come from `case.yaml`.

- **Hard sphere** (`gas.omega: 0.5`, the R-4D case): the correction is exactly
  one, λ = 1/(√2 π d² n), and T does not enter.
- **VHS** (0.5 < ω ≤ 1): the equilibrium VHS mean free path at the local
  *translational* temperature. It assumes the molecules there form one
  Maxwellian at that temperature. In an expanding plume (parallel and
  perpendicular temperatures diverge), in a shock layer, and where gas
  re-emitted from a wall meets the jet (two streams, whose collision rate is
  set by their relative speed rather than by T), λ is a nominal local length,
  not how far a molecule actually flies.
- `mean_free_path_convention: hard_sphere` drops the temperature factor.

## Gradients

`KnGradN` and `KnGradT` use a weighted least-squares gradient over each cell's
face neighbours (weights 1/|d|², OpenFOAM's `leastSquares`). It is taken of the
**averaged** fields — the gradient of the mean, never the mean of gradients —
and it is exact for a linear field on any stencil, including the 2:1 faces of
an octree mesh, where a Gauss gradient is not.

- Only **valid** neighbours enter the fit. A cell beside an empty one gets a
  one-sided gradient; a cell left with fewer than three non-coplanar valid
  neighbours gets none, and `diagValidGradN`/`T` is 0 there.
- Across a `symmetry` or `empty` patch the cell's mirror image is a neighbour,
  which imposes the zero normal gradient the plane implies. Every other
  boundary is one-sided.
- **Noise.** A gradient amplifies sampling noise: a relative noise ε in n
  puts roughly ε·λ/Δ into `KnGradN`. That is negligible where cells span many
  mean free paths and dominates where they are smaller than one — read a
  speckled `KnGradN` map against `sampledParcels`.
- **Mesh transitions.** A plume's density is not linear across a level change,
  and the lopsided stencil there gives a first-order error: expect a seam in
  `KnGradN` along every refinement-region boundary.
- A valid zero gradient gives a valid `KnGradN` of exactly zero.

## Validity

A cell's diagnostics are valid when it has **more than**
`min_sampled_particles` (20) sampled parcels, a positive finite n and a
resolved positive T. Otherwise `diagStatus` holds the sum of:

| flag | meaning |
|---|---|
| 1 | no parcel was ever sampled in the cell |
| 2 | sampled, but not more than `min_sampled_particles` |
| 4 | n ≤ 0 or not finite |
| 8 | no resolved T > 0: T ≤ 0, or its thermal part is round-off (below) |

What is counted is set by `particle_count_basis`:

- `per_step` (the default): the sampled — time-averaged — number of parcels
  in the cell, `dsmcRhoNMean`. It is how many *distinct* parcels the cell's
  moments come from at any step, and it is what "particles per cell" means
  everywhere else in this repository (`dsmc.target_particles_per_exit_cell:
  20`). At R-4D's shipped weight an exit cell holds about 8, so the exit region
  is reported invalid — with its bounds, under `validity.invalid_by_reason` —
  until the weight comes down to about 2e16.
- `accumulated`: parcels per step × steps averaged, every parcel sample behind
  the means. It is far larger and far weaker, because successive steps see the
  same parcels: a slow parcel counts once per step it lingers. In the run made
  by `tests/openfoam/test_nozzleplume_openfoam.py`, cells sampled by **one**
  parcel reached 26–149 "samples" this way and passed a threshold of 20.

The summary reports both counts either way. The threshold is a
diagnostic-confidence flag. It does not trigger refinement and does not certify
convergence.

**A temperature has to be resolved.** T is a difference, ⟨c²⟩ − |u|², of two
moments written to 10 significant figures. Where the difference is under 1e-8
of ⟨c²⟩ it is round-off, and the cell is flagged 8 whatever its count. That is
exactly what one parcel sampled many times produces (a single molecule has no
thermal spread; measured at ~1e-12), and no flow comes near it (1e-8 is 7e-5 K
at 3000 m/s). Without this check those cells reported T ≈ 1e-9 K and
`KnGradT` ≈ 1e11.

Invalid cells come in two kinds, and are written differently:

- **No value** (flags 1, 4 or 8: no parcels, n ≤ 0, no resolved T): the cell
  holds `invalid_value` (−1) in every written value field — a number every
  reader parses, and that no diagnostic can take, since all are ≥ 0.
- **Low confidence** (flag 2 alone: sampled, but not more than the threshold):
  the cell keeps its computed values, so the fields show the whole sampled
  flow, and `diagValid` is 0 there. Its gradient comes from its computable
  neighbours; a valid cell's gradient still comes from valid neighbours only.

Either way the cell is left out of every summary statistic, and **no value is
ever replaced by the exit state or by a neighbour's.** A run whose weight
leaves every cell under the threshold — R-4D at 5e18, with 0.003 parcels per
cell per step — therefore still shows its Knudsen numbers in ParaView, all
flagged, while its summary is empty and its console report says why. On a
logarithmic colour bar −1 takes the bottom colour, so read `diagValid` beside
the others.

## Outputs

`./Allpost` (step 4; `--no-diagnostics` skips it) or
`python -m plumetools.nozzleplume diagnostics [case] [--time T]` writes:

| | |
|---|---|
| `<time>/mfp`, `KnD`, … | `volScalarField`s in the sampled time directory, on its mesh: open `open.foam` in ParaView and they are beside `rhoNMean` |
| `results/diagnostics_summary.yaml` | below |
| `results/centreline_diagnostics.csv` | each diagnostic along the plume axis |
| `results/centreline-diagnostics.png` | the same, drawn: lengths, Δ/λ against 1/3, the Knudsen numbers |

A centreline value is the volume-weighted mean of the written field over the
cells `centreline.csv` samples at that station that have a value; the
`valid_cells` column says how many of those are trusted, and the figure draws
a station with none hollow. A station with no value at all is `nan`.

`./Allpost` then draws every diagnostic, beside the flow fields, on the `xz`
plane and the plate cross-section ([`viz.yaml`](../cases/r-4d/baseCase/viz.yaml));
`--no-viz` skips that and `./Allviz` redraws it. The mesh criterion comes
first: `mfpResolved` as a three-colour map, and `KnCellMax` = λ/Δ_max on a
diverging log scale with white at 3, pinned so iterations compare. Read
`mfpResolved` first — in `KnCellMax` a no-data cell takes the bottom colour.

The summary holds:

- `sampling` — the window (start, end, steps averaged) and the field names read;
- `settings` — the gas, the convention and every threshold used;
- `checks` — the particle-weight identity;
- `validity` — valid and invalid cell and volume fractions, by reason; how many
  cells got a gradient;
- `particles` — the sampled-parcel counts on both bases;
- `resolution` — statistics of Δ_max/λ and Δ_V/λ over valid cells; the cell
  and volume fractions above 1/k; a row per cell size (per octree level) with
  where its violations are; and the `worst_cells`, spread out so they name
  places rather than neighbours;
- `rarefaction` — statistics of `KnD`, `KnGradN`, `KnGradT`;
- `planning` — below;
- `skipped` — any diagnostic that could not run, and why;
- `notes` — what the numbers are and are not.

## Planning the next mesh

An octree mesh has only the sizes B/2^L of its background cell B. A practical
minimum cell h_min (`practical_min_cell_size_m`, 10 mm) allows the levels with
B/2^L ≥ h_min:

| level | B = 0.5 m (`baseCase`) | B = 1.0 m (`baseCaseHPC`) |
|---|---|---|
| 3 | 62.5 mm | 125 mm |
| 4 | 31.25 mm | 62.5 mm |
| 5 | **15.625 mm** — finest permitted | 31.25 mm |
| 6 | 7.8125 mm — excluded | **15.625 mm** — finest permitted |
| 7 | | 7.8125 mm — excluded |

The minimum and `mesh_budget_cells` (10 million) are planning constraints. They
decide what is listed as available; they never clip the reported λ/k target or
hide a violation. The summary's `planning` block gives:

- `refinement_levels` — the table above for this case;
- `required_level_histogram` — valid cells and volume by the level that would
  meet λ/k, marked permitted or not;
- `target_below_practical_minimum` — cells whose λ/k is under h_min, with their
  bounds and the smallest target;
- `mesh_budget` — the current cell count, and lower-bound estimates of a mesh
  that met λ/k in every valid cell: stopped at the finest permitted level, and
  uncapped. Buffer layers are not included.

At the R-4D exit (λ ≈ 0.14 mm, so λ/3 ≈ 0.05 mm) no permitted level comes close;
expect the exit region to appear in `target_below_practical_minimum` however the
mesh is refined.

## Starting a fresh iteration

```bash
cd cases/r-4d
python -m plumetools.iteration baseCase iter01-exit5 --note "exit cone to level 5"
#   edit iter01-exit5/case.yaml
cd iter01-exit5 && ./Allmesh --describe && ./Allmesh && ./Allrun
```

[`plumetools/iteration.py`](../plumetools/iteration.py) copies the files at the
top of the case that make it what it is — `case.yaml`, other `*.yaml`, the
`All*` scripts, `*.foam` — and nothing else. The dictionaries, the mesh, `0/`,
and the time step and particle occupancy derived from the smallest cell are all
generated again by `./Allmesh` and `./Allrun` from the new `case.yaml`. Time
directories, `processor*`, `lagrangian/` parcels and `uniform/` averaging state
are never copied — not even by `--include`, which is for genuine extra inputs
such as an STL — and the new directory is checked for all of them before the
command returns. With no time directory, `controlDict`'s `startFrom latestTime`
resolves to 0 and `./Allrun` runs `dsmcInitialise`: new parcels, new averaging.

The source is only read. The destination must be new (or empty), and may not be
inside the source. Put iterations **beside** the base case
(`cases/r-4d/iter01-…`): the repository's ignore rules for `cases/r-4d/*/`
then keep their generated files out of git.

`iteration.yaml` records the parent, the hash of the `case.yaml` copied, the
parent's latest time, the lineage, the repository commit and your note.
`provenance/` keeps the parent's `case.yaml` and diagnostics summary:

```bash
diff provenance/parent-case.yaml case.yaml     # exactly the edit that made it
```

Each diagnostics summary carries its case's `case_yaml_sha256` and, for an
iteration, its name and parent, so two summaries can be lined up.

## Another case

A **nozzle-plume** case (`model: nozzle_exit_plume`) needs nothing: its
`./Allpost` runs the diagnostics with the gas, convention, exit diameter and
background cell taken from `case.yaml`. An optional `diagnostics:` section
changes the defaults:

```yaml
diagnostics:
  min_sampled_particles: 20
  particle_count_basis: per_step        # or accumulated
  cells_per_mean_free_path: 3.0
  practical_min_cell_size_m: 0.010
  mesh_budget_cells: 10000000
  reference_lengths_m: {Standoff: 20.0} # adds KnStandoff = lambda / 20 m
  select: null                          # or a list of output field names
  field_names: {}                       # e.g. {number_density: rhoNMean}
  invalid_value: -1.0
```

**Any other dsmcFoam case** uses the generic command with a small settings file:

```yaml
# diagnostics.yaml
gas: {molecular_mass_kg: 6.63e-26, molecular_diameter_m: 4.17e-10,
      omega: 0.81, t_ref_K: 273.0, mean_free_path_convention: vhs}
background_cell_size_m: null       # a number enables the level planning
diagnostics:
  reference_lengths_m: {D: 0.2}
```

```bash
postProcess -func writeCellCentres -latestTime
postProcess -func writeCellVolumes -latestTime
python -m plumetools.rarefaction <case> --settings diagnostics.yaml
```

or `plumetools.rarefaction.run.diagnose(...)` from Python.

## Adding a diagnostic

One decorated function in
[`plumetools/rarefaction/core.py`](../plumetools/rarefaction/core.py), naming
what it needs:

```python
@register("KnGradP", requires=("mean_free_path", "pressure", "grad_pressure"),
          dimensions=DIMENSIONLESS, units="-",
          description="lambda |grad p| / p")
def _kn_grad_p(ctx):
    ...
    return Field.masked(values, valid)
```

It is then computed, written, summarised (names starting `Kn` land in the
`rarefaction` block) and selectable with `diagnostics.select`. A requirement
that is a new input is read in `plumetools.rarefaction.foam.load`; until it is,
the diagnostic is skipped with the reason and the rest run. A length scale needs
no code at all: `reference_lengths_m`.

## Limitations

- Reconstructed cases only; `processor*` directories are not read.
- ASCII only, like the rest of this repository's readers.
- The whole mesh is read into memory. On the 459,007-cell R-4D mesh the
  diagnostics take 11 s and 450 MB; memory grows with the cell count, so a
  10-million-cell mesh needs roughly 8 GB — more than a laptop's WSL may
  spare, comfortable on the HPC node `baseCaseHPC` is for.
- One species. A mixture would need per-species λ and a mixture mean free path.
- λ is the equilibrium value (see above); the gradient is first-order at
  level transitions and noise-sensitive where cells are smaller than λ.
- The mesh estimate is a lower bound — no buffer layers — and assumes cells are
  only refined, never coarsened.
- The worst-violation list and the per-size table are descriptive. Turning them
  into `refinement_regions` is the person's decision.

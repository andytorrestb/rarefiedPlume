# Portable DSMC mesh refinement

`plumetools.refinement` assesses measured local density, translational temperature
and mean parcels per cell. Its default stopping targets are largest actual edge
`<= 0.8 lambda_measured` and mean occupancy `>= 20`. Passing cells do not request
another refinement level. This utility requires no Cai, nozzle or radial metadata.

## An unrelated OpenFOAM case

Use a reconstructed, static, single-species ASCII case. Generate actual geometry
fields in a copy if they are absent:

```bash
postProcess -case /data/my-case -func writeCellCentres -time 1
postProcess -case /data/my-case -func writeCellVolumes -time 1
```

Save standalone reader settings, for example `refinement.yaml`:

```yaml
gas:
  species_name: Ar
  molecular_mass_kg: 6.63e-26
  molecular_diameter_m: 4.17e-10
  omega: 0.74
  t_ref_K: 273
  mean_free_path_convention: vhs
fields:
  profile: custom
  number_density: myNumberDensity
  temperature: myTranslationalTemperature
  mean_parcels: myMeanParcelCount
  cell_centres: C
  cell_volumes: V
sampling:
  start_s: 0.5
  end_s: 1.0
  basis: mean_per_cell
  source: declared
targets:
  edge_over_mfp: 0.8
  min_mean_parcels: 20
```

The declared interval must describe these fields' actual averaging window. Counts
must be mean per-cell occupancy, not accumulated samples. `profile: standard`
reads stock averaged moments and fieldAverage metadata; `profile: hystrath` reads
`rhoN_Ar`, `Ttra_Ar`, and `dsmcNMean_Ar`, with an optional species `suffix`.
Custom names can override each role in either profile. The standard reader derives
temperature from moments; the hyStrath reader supplies direct temperature.

```bash
python -m plumetools.refinement assess /data/my-case --settings refinement.yaml --time 1
python -m plumetools.refinement plan /data/my-case --settings refinement.yaml --time 1 --output /data/targets-01
python -m plumetools.refinement apply /data/my-case --plan /data/targets-01 --backend hex --output /data/mesh-01 --max-cells 1000000
```

Assessment is read-only unless `--output` is supplied. `--dry-run` on assess/plan
writes nothing. A plan contains JSON provenance, a checksummed NPZ of source
geometry and per-cell bounds, and a CSV for another mesher. Changed source fields
or mesh invalidate it. Insufficient sampling or low occupancy prevents a trusted
mesh correction and identifies the missing work.

The `hex` backend requires sourced OpenFOAM v2512 and a supported 3-D hex mesh.
It splits failing cells with `refineHexMesh`, extends the selection for conformity,
regenerates `C`/`V`, and checks actual edges and mesh quality after refinement.
Its discrete splits may overshoot the target; no additional target ratio is
introduced. Budgets stop publication instead of coarsening the requested mesh.
Passing plans return no new mesh artifact. Existing output directories are refused.

The new artifact contains geometry and inherited planning fields, **not child
measurements**. Install it in a fresh case through that solver's own initializer,
choose a new parcel weight, run, and reassess. No parent parcel cloud is mapped.
Source mesh, fields and results remain intact. Failure logs remain available.

Unsupported binary, dynamic, decomposed-only, mixture, or incompatible topology
inputs receive an actionable error. Assessment and CSV target export can still
support another mesher when the provided hex backend cannot materialize a mesh.

## Python and other meshers

The pure API accepts `Geometry`, `Measurements`, gas properties and `Targets`.
See the [API](reference/api/refinement/index.md) for a runnable example. Other
solvers supply these normalized arrays without OpenFOAM dictionaries.

A mesher implements `MeshRefiner.refine(source, plan, destination)` and returns
`RefinementResult`. Preserve source identities, verify the plan's hashes, honor
all valid target edges, publish a fresh artifact, and require a new flow run.
The [radial study](cai2012-radial-resolution.md) supplies a second backend using
the same target-plan contract.

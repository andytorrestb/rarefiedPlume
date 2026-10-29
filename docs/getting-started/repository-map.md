# Repository map

```text
plumetools/                 the Python package (static analysis builds its API pages)
  sourceflow.py, inflow.py,   the LEGACY source-flow model and its orchestration (3d-inflow),
  config.py, geometry.py,     frozen by the regression golden; config.py also loads
  species.py, constants.py,   markelov1999's case.yaml
  verify_mesh.py
  cai2012/                  Cai & Wang 2012: config, gas, analytical, mesh, study, checks,
                            post, health (statistical sweep), audit
  cai2012ext/               the Kn extension: policy, generation, feasibility, comparisons,
                            hpc wrappers; interface_* (interface track); scaling_* (field scaling)
  markelov1999/             AIAA 99-3455: paper-faithful source flow, geometry, mesh, flux,
                            resolution, post (frozen package)
  nozzleplume/              general nozzle-exit plume with an optional plate (r-4d)
  hybridplume/              one-way 3-D CFD -> DSMC hybrid
  virtualcollision/         reduced virtual-collision model and its assessment
  rarefaction/              mean free path, Knudsen numbers, mesh resolution after a run
  viz/                      ParaView slice imagery; video.py and dissect.py need only ffmpeg
  foamio/                   OpenFOAM field/dictionary I/O, time directories, checkpoints
  mesh/                     polyMesh parsing and patch geometry
  casemanifest.py           the case list study runners read, checked whole
  iteration.py              a fresh mesh iteration of a case: inputs only
applications/
  dsmcBoundaryModels/       plumeFieldInflow -> libplumeDsmcBoundaryModels.so
  thermophysicalModels/     hybridPlumeThermo -> libhybridPlumeThermo.so
cases/                      OpenFOAM cases and studies (see the case catalog)
  ARCHIVE.md                the frozen cases and what is wrong in them
docs/                       this site's sources; docs/_tools/ builds and checks it
tests/unit/                 Tier 0: pure Python
tests/regression/           bit-exact freeze of the legacy 3d-inflow inflow
tests/openfoam/             Tier 2: needs OpenFOAM (marker needs_openfoam)
util/                       one-off analysis scripts from the archived lineage (not maintained)
Allsetup                    install, build the inflow library, verify
mkdocs.yml                  documentation site configuration
```

## Conventions you will meet

- **Units are in key names**: `diameter_m`, `T0_K`, `p0_pa`,
  `number_density_per_m3`, `write_interval_s`. Normalised quantities say so:
  `x_over_D`, `n_over_n0`, `U_sqrt_beta0`.
- **Values are tagged** in the committed YAML: `[PAPER]` (printed in the
  reference), `[ASSUMPTION]` (an implementation choice, documented),
  `[DERIVED]` (computed, never typed), `LEGACY` (a frozen pre-refactor value,
  with its finding ID), `[EXTENSION]`, `[HEALTH]`, `[BIRD]`, `[NUMERICAL]`,
  `[REPO]`.
- **`null` means derive it**, in the Cai and nozzle-plume schemas.
- **Generated is gitignored.** OpenFOAM dictionaries, `0/` fields, generated
  case trees, manifests and `results/` are reproducible from the templates and
  are not committed — with listed exceptions (measured `cost-model.yaml`
  files, the virtual-collision smoke result, the Pointwise meshes of the cases
  that have no generator).
- **Finding IDs** (`SM-`, `AR-`, `HA-`, `RP-`, `SF-`, `AD-`) name recorded
  defects; they appear in `case.yaml` comments and in
  [The plume source-flow model](../source-flow-model.md) and
  [`cases/ARCHIVE.md`](../../cases/ARCHIVE.md).

## Generated, tracked, and where results go

| Kind | Example | Tracked? |
|---|---|---|
| inputs | `case.yaml`, `study.yaml`, `hybrid.yaml`, `viz.yaml`, `interfaces.yaml` | yes |
| runner scripts | `Allmesh`, `Allrun`, `Allpost`, `Allclean`, `AllrunCases`, `run*.sh` | yes |
| generated dictionaries and fields | `system/controlDict`, `constant/dsmcProperties`, `0/` | no |
| generated case trees | `Cases/`, `Cases_full_suite/`, `Interface/`, `runs/` | no |
| meshes | `constant/polyMesh` | **yes** for the Pointwise-meshed legacy and archived cases (no generator exists); no for generated ones |
| solver output | time directories, `processor*`, `log.*` | no |
| results | `results/` | no, except the files listed above |

## Glossary

| Term | Here |
|---|---|
| parcel | a simulated DSMC particle, standing for `nEquivalentParticles` real molecules (the particle weight) |
| exit cell | the cell at the nozzle exit where the occupancy target is met |
| transit | the time for the flow at U₀ to cross the domain once; transients and sampling windows are counted in transits |
| t₀, τ | Cai's reference collision time; the mean collision time λ/c̄ |
| manifest | the generated record of every case of a study and its derived inputs |
| suite profile | `meta.knudsen_suite` in the Kn extension's `study.yaml`: `full_decades` or the historical default |
| NON-PRODUCTION | the label written into a case that misses a production criterion; the case still runs |
| dialect | `standard` (OpenFOAM `dsmcFoam`) or `mnf` (the fork's `dsmcFoam+`) |

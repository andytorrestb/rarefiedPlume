# markelov1999: a plume past a cylinder onto a flat plate (AIAA 99-3455)

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | The configuration of Lumpkin, Stewart & Markelov (1999): a rarefied N₂ plume from a small orifice past a finite cylinder onto a flat plate 6 in behind it, at reservoir pressures 5, 25, 100 and 475 psi |
    | **Reference** | Lumpkin, F. E., Stewart, B. D., and Markelov, G. N., *Study of 3D Rarefied Flow on a Flat Plate in the Wake of a Cylinder*, AIAA 99-3455, 1999 |
    | **Model** | A paper-faithful source-flow model on a hemispherical inflow surface (`markelov1999_axisymmetric`, deliberately separate from the frozen legacy model), injected by `plumeFieldInflow` with a per-face density; `snappyHexMesh` geometry; VHS + Larsen–Borgnakke, Z_R = 5 |
    | **Evidence** | Unit-tested model, geometry, flux verification and post-processing; **solver-tested** in `tests/openfoam` (mesh, library loading, patch selection, open boundaries, flux consistency). **No physical result is validated** against DAC or experiment |
    | **Runs on** | Python + **CaseFoam** for study generation (including `--dry-run`); OpenFOAM v2512 + `plumeFieldInflow` for meshing and solving |
    | **Start in** | `cases/markelov1999/` (study) or `cases/markelov1999/Cases/gap06in/<case>/` |

This is new work, not a fix of the archived `wake-cylinder` lineage, and it
does not inherit that lineage's defects. The model, the six places the legacy
code disagrees with the paper (SF-1 … SF-6), and every assumption with its
reasoning are in [The AIAA 99-3455 family](../markelov1999-case.md).

## Geometry and coordinates

Printed by the paper: orifice diameter 0.8255 mm; cylinder 6 in diameter ×
18 in long, centred at 11.75 in; plate 6 in wide × 15 in high; 6 in gap between
cylinder and plate. Assumptions (each tagged `[ASSUMPTION]` in
`baseCase/case.yaml`): inflow hemisphere radius 0.1524 m, domain extents,
plate thickness 0.5 in, background cell size and refinement levels, VHS
coefficients, stagnation temperature 300 K, time step, particle weighting, the
pressure-averaging windows.

The configuration is symmetric about y = 0 only. Patch roles
(`mesh.patch_names`): `inflow`, `cylinder`, `plate`, `outer`, `symmetry`, and
`upstream_vacuum` — the x = 0 plane outside the inflow cavity, which is an
**open** boundary, not a symmetry plane. Pressure windows use azimuth 180° for
windward (facing the source) and 0° for leeward.

## Prerequisites

| Step | Needs |
|---|---|
| `generate_cases.py` (and `--dry-run`) | **CaseFoam** (`pip install -e ".[cases]"`; verified against casefoam 0.2.0). It is a dependency, not an option: there is no built-in substitute |
| `Allmesh --dict-only` | Python with `plumetools` |
| `Allmesh` | `blockMesh`, `snappyHexMesh`, `checkMesh` |
| `Allrun` | `dsmcInitialise`, `dsmcFoam`, **`libplumeDsmcBoundaryModels.so`** — every case selects `plumeFieldInflow` and refuses to run without it; there is no fallback to a uniform inflow |
| `Allpost` | `writeCellCentres` (optional — without it the resolution audit is domain-wide only); matplotlib for figures |

CaseFoam is not installed in this checkout's environments, so study generation
was not exercised for this documentation.

## Inputs

| File | Role |
|---|---|
| `study.yaml` | `gap_in` (6.0) and `gap_dir` (`gap06in`); four `pressure_cases` (`name`, `pressure_psi`, `enabled`) |
| `baseCase/case.yaml` | **the authoritative physics**, loaded by `plumetools.config` with `mesh.type: snappy_markelov` and the `bodies`, `dsmc`, `resolution`, `checks` and `post` sections; cloned per case with the pressure written in structurally |
| `viz.yaml` | slice imagery |

`constant/dsmcProperties`, `system/controlDict`, `system/dsmcInitialiseDict`,
`system/decomposeParDict`, the mesh dictionaries and every `0/` field are
**generated** from `case.yaml` and gitignored. Do not edit them. This replaces
the archived string-substitution pattern (`util/caseFoam/genCases.py`), which
silently did nothing when a shipped file was reformatted.

Schema: [configuration reference](../reference/configuration.md#legacy-source-flow-caseyaml-3d-inflow)
(the same loader; Markelov-specific sections are listed there).

## Cheap checks (no OpenFOAM, no CaseFoam)

The template can be checked without generating the study — **in a copy**, so
the generated dictionaries do not land in `baseCase/` (where CaseFoam would
clone them into every case; generation removes them, but there is no reason
to create them):

```bash
cp -r cases/markelov1999/baseCase /tmp/mk-check && cd /tmp/mk-check
./Allmesh --dict-only      # model, geometry, mesh plan, flux check, DSMC checks; writes system/ and constant/ dictionaries
```

With CaseFoam installed:

```bash
cd cases/markelov1999
./generate_cases.py --dry-run          # the cases, their pressures in Pa, their destinations
./generate_cases.py --pressures 5 25   # a subset
```

## Production run

```bash
cd cases/markelov1999
./generate_cases.py      # CaseFoam clones baseCase into Cases/gap06in/<case>; values applied structurally
./AllmeshCases           # blockMesh + snappyHexMesh + checkMesh + geometry verification, per case
./AllrunCases            # sequentially, in manifest order
./AllpostCases           # per-case summaries, the study table, plots, imagery
```

!!! warning "CaseFoam's own `Allrun` launches every case at once"
    CaseFoam writes `Allrun`, `Allclean` and `rmCases` into the study directory
    (gitignored). Its `Allrun` starts every case **concurrently** in the
    background — four simultaneous multi-million-parcel runs. Use
    `./AllrunCases`, which runs them one after another.

Per case, from `Cases/gap06in/<case>/`: `./Allmesh` (or `--dict-only`,
`--force`), `./Allrun` (or `--serial`, `-np N`, `--no-solve` — which writes `0/`
and verifies the inflow flux without the solver, and needs the mesh),
`./Allpost`. Study runners: `AllmeshCases --dict-only --keep-going`,
`AllrunCases --serial -np N --no-solve --keep-going`,
`AllpostCases --no-viz --viz-only --viz-spec F --keep-going`.

`./Allrun` **continues** from the newest checkpoint complete on every rank,
exactly as in the [Cai family](cai2012.md#continuing-a-run): `dsmcInitialise`
and `decomposePar -force` are skipped on a resume, `reconstructPar -newTimes`
adds new times, and `fieldAverage` keeps accumulating — which matters here,
because the surface pressures come from those means. To run further, raise
`dsmc.end_time_s`, then `./Allmesh --dict-only && ./Allrun`.

!!! danger "Commands that destroy results"
    `./AllcleanCases --purge` (the whole `Cases/` tree and the manifest),
    `./Allclean` / `--mesh` per case, `./generate_cases.py --overwrite`, and
    `./Allmesh --force` on a case with results. Preview cleaning with
    `--dry-run`.

## Outputs

| Path | Contents |
|---|---|
| `manifest.yaml` | the auditable case → inputs map, including every `[ASSUMPTION]` |
| `Cases/gap06in/<case>/case-summary.json` | `runInflow.py` writes the inflow, the flux verification and the DSMC check results; `Allpost` keeps those and adds `mesh`, `particles`, `pressure` (one entry per window) and `resolution_audit` (measured parcels per cell, per region when cell centres exist, against the target) |
| `results/study-table.csv` | one row per case, from the summaries |
| `results/windward_pressure.png`, `leeward_pressure.png`, `pressure_ratio.png` | the plots |

**Surface pressure** is the wall-normal component of the time-averaged `fDMean`
(the momentum the gas delivers to a face, per area per time), area-averaged
over each window — not `n k T` from the cell next to the wall, which omits the
directed momentum of the plume. The derivation is in
[§3.9a](../markelov1999-case.md#39a-surface-pressure-from-fd).

**Judging the output.** Check the parcel count in `log.dsmcFoam` has levelled
off before reading any pressure: the shipped `end_time_s` is an assumption, not
a demonstrated steady state. Read the resolution audit: at the shipped settings
the 20 particles-per-cell target is met at the cylinder by construction and not
in the wake or at the plate (one global weight; see
[§3.8](../markelov1999-case.md#38-particle-weighting)). No pressure, ratio or
profile here has been compared with the paper.

What has been established by running against OpenFOAM v2512 (case README and
`tests/openfoam/test_markelov_openfoam.py`): the mesh builds and `checkMesh`
passes; the mesh is the declared geometry, measured from its own vertices
(including the 6.000 in gap); the library loads and injects on the one named
patch; it reads a per-face density; particles leave through the open
boundaries; two independent flux routes agree to round-off.

## Resources

Configured, at the shipped settings (case README, from the resolution
estimator): about 3.4 × 10⁶ particles, independent of reservoir pressure;
`resolution.sizing_region: plate` meets the target everywhere at roughly six
times that. No measured run time is recorded.

## Variants and safe customisation

| Want | Do |
|---|---|
| One pressure | `enabled: false` on the others, or `--pressures` |
| Meet the occupancy target at the plate | `resolution.sizing_region: plate` (≈ 6× the particles) |
| A different gap | `gap_in`/`gap_dir` in `study.yaml` and the body geometry in `case.yaml`; regenerate |

## Figure placeholder

!!! figure "Figure placeholder — windward plate pressure against reservoir pressure"
    **Subject:** `results/windward_pressure.png`.
    **Source run:** all four `Cases/gap06in/p*psi`, solved past steady state and post-processed.
    **Plotting command:** `./AllpostCases` (writes the three PNGs; verified in source).
    **Destination:** `docs/assets/figures/markelov1999-windward-pressure.png`
    **Suggested caption:** "Windward window pressure on the plate, four reservoir pressures (not compared with the paper)."

## Related

- [The AIAA 99-3455 family](../markelov1999-case.md) — paper values, the model, assumptions, flux verification, checks, CaseFoam
- [plumeFieldInflow](../plume-field-inflow.md) · [Solver compatibility](../solver-compatibility.md)
- [Archive](../archive/index.md) — the frozen `wake-cylinder` lineage
- API: [`markelov1999.sourceflow`](../reference/api/markelov1999/sourceflow.md), [`flux`](../reference/api/markelov1999/flux.md), [`study`](../reference/api/markelov1999/study.md), [`postprocess`](../reference/api/markelov1999/postprocess.md)

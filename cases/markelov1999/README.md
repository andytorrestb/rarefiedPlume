# AIAA 99-3455 — flat plate in a cylinder wake

Lumpkin, F. E., Stewart, B. D., and Markelov, G. N., *Study of 3D Rarefied Flow
on a Flat Plate in the Wake of a Cylinder*, AIAA 99-3455, 1999.

A rarefied nitrogen plume expands from a small orifice, passes a finite cylinder,
and impinges on a flat plate 6 inches behind it. Four reservoir pressures — 5,
25, 100 and 475 psi — at one cylinder-to-plate separation.

> **These results are not validated.** Nothing here has been compared with DAC or
> with experiment. Digitising the paper's curves is out of scope, and no claim of
> agreement is made anywhere in this directory. See
> [Status](#status-what-this-does-and-does-not-establish).

## Quick start

```bash
pip install -e ".[test]"                              # from the repository root

cd applications/dsmcBoundaryModels && ./Allwmake      # the custom inflow model
cd ../../cases/markelov1999

./generate_cases.py --dry-run                         # what would be generated
./generate_cases.py                                   # the four pressure cases
./AllmeshCases                                        # blockMesh + snappyHexMesh
./AllrunCases                                         # solve
./AllpostCases                                        # summaries, table, plots
```

Per case, from inside `Cases/gap06in/p005psi`:

```bash
./Allmesh              # dictionaries + mesh + checkMesh + geometry verification
./Allmesh --dict-only  # dictionaries only; needs no OpenFOAM
./Allrun --serial      # one process
./Allrun -np 8         # override the subdomain count
./Allrun --no-solve    # generate 0/ and verify the flux; needs no solver
./Allpost              # statistics, wall pressures, resolution audit
./Allclean             # remove solver output (--mesh to remove the mesh too)
```

## Layout

```
cases/markelov1999/
  study.yaml            the case matrix: four pressures, one gap
  generate_cases.py     clone baseCase per pressure, write manifest.yaml
  baseCase/             the template — case.yaml plus the scripts
    case.yaml           THE authoritative source of physical inputs
  Cases/gap06in/        generated; gitignored
    p005psi/ p025psi/ p100psi/ p475psi/
  manifest.yaml         generated; the auditable case → inputs map
  results/              generated; study-table.csv and the plots
```

`Cases/`, `manifest.yaml` and `results/` are gitignored: they are reproducible
from `baseCase` and `study.yaml`, and committing a generated tree would make the
template and its copies drift apart.

## `case.yaml` is the only place physics lives

`constant/dsmcProperties`, `system/controlDict`, `system/dsmcInitialiseDict`,
`system/decomposeParDict`, the two `blockMesh`/`snappyHexMesh` dictionaries and
every `0/` field are **generated** from `case.yaml`. They are gitignored and
overwritten on every `./Allmesh`. Do not edit them — edit `case.yaml`.

This replaces the pattern in `util/caseFoam/genCases.py`, which patches shipped
dictionaries by string substitution:

```python
'system/controlDict': {'#!stringManipulation':
                        {'deltaT          1.0E-05': '%s' % deltaT}}
```

That silently does nothing if the shipped file is reformatted, if the value has
already changed, or if the case was generated from an older template — and none
of those is visible in the output. `generate_cases.py` instead loads each
`case.yaml` as a data structure, sets specific keys, and re-emits it. A key that
is not there is an error.

## What CaseFoam does

[CaseFoam](https://github.com/DLR-RY/caseFOAM) clones `baseCase` into the
`Cases/gap06in/<case>` hierarchy. It is a **dependency**, not an optional
accelerant — there is no built-in substitute:

```bash
pip install -e ".[cases]"          # casefoam 0.2.0
./generate_cases.py
```

`generate_cases.py` calls `casefoam.mkCases` the way it is designed to be
called — pointed at the `baseCase` *directory*:

```python
mkCases(<baseCase dir>, [["gap06in"], ["p005psi", ...]],
        caseData, hierarchy="tree", writeDir="Cases")
```

It copies the template to `Cases/`, moves that content down into
`Cases/baseCase/`, and creates `Cases/gap06in/<case>/` from it. So
`Cases/baseCase` alongside the case directories is **CaseFoam's own layout** —
the same one `cases/caseFoamEx` has. It is not cleaned up; it belongs to
CaseFoam.

> **CaseFoam also writes `Allrun`, `Allclean` and `rmCases` into this
> directory** (`mkAllRunClean` opens `'Allrun'` relative to the working
> directory). They are kept, and gitignored. Note what its `Allrun` does:
>
> ```
> Cases/gap06in/p025psi/Allrun &
> Cases/gap06in/p005psi/Allrun &
> ...
> ```
>
> — every case launched **concurrently**. Four simultaneous 3.4 M-particle DSMC
> runs is how a study of four cases becomes four cases that all get killed. Use
> **`./AllrunCases`**, which runs them sequentially in manifest order. CaseFoam's
> is left in place for anyone who does want the parallel launch.

**Why the parameters are applied outside CaseFoam.** `caseData` is passed empty.
Its two forms both miss `case.yaml`: the dictionary-aware form goes through
PyFoam's `ParsedParameterFile`, which reads OpenFOAM dictionaries and not YAML,
and the other is `'#!stringManipulation'` — the whitespace-sensitive substitution
described above. So CaseFoam does the cloning and the hierarchy, and
`apply_case_parameters` applies the physical values structurally.

The one thing done afterwards is removing generated dictionaries from each clone.
That is template hygiene, not a second cloner: CaseFoam copies the template
faithfully, as it should, so a `constant/dsmcProperties` left in `baseCase` by
someone running `./Allmesh` there would otherwise let every case run against the
template's physics instead of its own.

## Status: what this does and does not establish

**Established, by running it against OpenFOAM v2512:**

| | |
|---|---|
| The mesh builds | 12 096 background cells → 34 346 after snapping |
| `checkMesh` | **Mesh OK** — max non-orthogonality 37.6, mean 10.2, max skewness 0.69 |
| The mesh is the declared geometry | every dimension measured from the mesh's own vertices, including the 6.000 in gap between the cylinder and plate patches |
| The custom inflow model loads | `Selecting InflowBoundaryModel plumeFieldInflow` |
| It injects on one named patch | `injecting across 1(inflow)`, not on every `patch`-type boundary |
| It reads a per-face density | the field varies across the inflow surface rather than being collapsed to a mean |
| Particles leave the domain | the count is below everything ever supplied, which reflection would forbid |
| Inflow flux conservation | two independent deterministic routes agree to 6e-15; the field round trip to 9e-12 |

**Not established:**

* **Any physical result.** No pressure, ratio or profile here has been compared
  with the paper, with DAC, or with experiment.
* **Steady state.** The shipped `end_time_s` is an assumption, not a converged
  result. Check the particle count in `log.dsmcFoam` before reading any pressure.
* **The 20 particles/cell target.** See below.

## The particle-resolution trade

Standard `dsmcFoam` has **one** particle weight for the whole domain —
`DSMCCloud::nParticle_` is a single scalar, read from `nEquivalentParticles`.
There is no radial, per-cell or adaptive weighting to configure.

So for a flow whose density falls as `1/r²`, the occupancy target can be met in
exactly **one** region; every other region follows from

```
occupancy(A) / occupancy(B) = [n(A) · V_cell(A)] / [n(B) · V_cell(B)]
```

which is a property of the flow and the mesh, not of anything configurable.

At the shipped settings, sized on the cylinder:

| region | particles/cell | note |
|---|---|---|
| cylinder | **20.0** | the target, met by construction |
| wake | 6.2 | 1 in behind the base |
| plate | 3.6 | |
| total | ~3.4 × 10⁶ particles | independent of reservoir pressure |

Setting `resolution.sizing_region: plate` meets the target everywhere at roughly
six times the particle count. Whichever is chosen, the estimator reports all
three regions, and `./Allpost` audits what actually happened against the sampled
`dsmcRhoN` field — which holds the literal parcel count per cell. **The target is
never reported as met on the strength of the estimate alone.**

## Paper values versus assumptions

Everything below is stated again, with reasoning, in
[`docs/markelov1999-case.md`](../../docs/markelov1999-case.md).

### Printed by the paper

orifice diameter 0.8255 mm · cylinder 6 in diameter × 18 in long, centred at
11.75 in · plate 6 in wide × 15 in high · 6 in gap · reservoir pressures 5, 25,
100, 475 psi · nitrogen · VHS collisions with Larsen-Borgnakke redistribution ·
Z_R = 5

### Assumptions

inflow hemisphere radius (0.1524 m) · domain extents · plate thickness
(0.5 in) · background cell size and refinement levels · the VHS coefficients
modern OpenFOAM needs but the paper does not print · stagnation temperature
(300 K) · time step · particle weighting · the pressure-averaging windows

Every one is tagged `[ASSUMPTION]` in `baseCase/case.yaml` and carried into
`manifest.yaml` and each `case-summary.json`.

## Surface pressure

There is no ready-made wall-pressure field in standard `dsmcFoam`. What it
provides is `fD`, accumulated in `DSMCParcel::hitWallPatch`:

```cpp
deltaFD = cloud.nParticle()*(preIMom - postIMom)/(deltaT*fA);
cloud.fDBF()[wppIndex][wppLocalFace] += deltaFD;
```

— the momentum the gas delivers to the face, per unit area per unit time.
Dimensions `[1 -1 -2 0 0 0 0]` = Pa. The surface pressure is its wall-normal
component, `p = fD · n̂`; the tangential part is the shear stress and is reported
separately rather than folded in.

This is **not** the gas static pressure. Substituting `n k T` from the cell next
to the wall would omit the directed momentum of the impinging plume, which at
these speed ratios is most of the load.

`fD` is reset every timestep by `DSMCCloud::resetFields`, so a written snapshot
is a one-step momentum tally. The generated `controlDict` runs `fieldAverage`
from `dsmc.average_start_s`, and the post-processing reads `fDMean`.

Pressure is averaged over an angular and axial **window**, area-weighted — never
one face, whose value would depend on where snappyHexMesh happened to put it.
Windward is azimuth 180° (facing the source); leeward is 0° (the base).

## Related documents

* [`docs/markelov1999-case.md`](../../docs/markelov1999-case.md) — the model, the
  geometry, and every assumption with its reasoning
* [`docs/plume-field-inflow.md`](../../docs/plume-field-inflow.md) — the custom
  boundary model and why standard `dsmcFoam` needs one
* [`docs/solver-compatibility.md`](../../docs/solver-compatibility.md) — what
  stock `dsmcFoam` cannot express
* [`cases/ARCHIVE.md`](../ARCHIVE.md) — the frozen `wake-cylinder` lineage, which
  attempted this configuration and whose known defects this family does not
  inherit

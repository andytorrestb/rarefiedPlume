# 3d-inflow: the legacy source-flow reference case

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | The reference case of the 3D source-flow lineage: an analytical plume model evaluated per face on a hemispherical inflow surface, DSMC downstream |
    | **Model** | The *legacy* source-flow model, `plumetools.sourceflow` (equations E1–E8), with every known defect frozen as a `LEGACY` value |
    | **Evidence** | Mechanically tested (unit tests, and a bit-exact regression golden against the pre-refactor script). **Not physically validated**; the model has no recoverable citation |
    | **Runs on** | Python for `0/` generation; OpenFOAM v2512 `dsmcFoam` for meshing and solving |
    | **Start in** | `cases/3d-inflow/` |
    | **Blocker** | A fresh checkout has **no mesh**; making one needs OpenFOAM or a `git checkout` from history ([below](#get-a-mesh)) |

## What this case is for

A rarefied plume expands from a hemispherical inflow surface into a vacuum box.
The inflow state on each face — number density, velocity, temperature — comes
from the analytical source-flow model; `dsmcFoam` solves the field from there.
Seven older cases (`iss_solar_panels{,_2,_3}`,
`solar_panel_particle_resolution_study/Fnum_*`) were copied from it; see the
[ISS tutorial](iss-solar-panels.md).

Use it to understand the legacy model and the regression machinery. Do not use
it for new physics: the model disagrees with its own stated inputs in several
places (SM-01 … SM-10), each frozen deliberately so the refactor could be proven
behaviour-preserving. The newer families replace it — `markelov1999` with a
paper-faithful source model, `cai2012` and `nozzleplume` with a nozzle exit as
the DSMC inlet.

Read [The plume source-flow model](../source-flow-model.md) before changing any
`LEGACY` value. The model's provenance is unknown; the only literature reference
in the repository covers the archived 1d/2d-wedge validation case.

## Geometry and coordinates

| | |
|---|---|
| Domain | box, x ∈ [0, 5] m, y, z ∈ [−2.5, 2.5] m |
| `inflow` | hemisphere of radius 0.5 m centred at the origin, opening toward **+x** |
| `sym` | the x = 0 plane outside the hemisphere |
| `vacuum` | the five other box faces |
| Plume axis | **+x** — inferred from the geometry and `revolutionAxis "x"`, never stated in the original code |
| Polar angle θ | measured from **+z**, not from the plume axis (SM-04, frozen) |

Generated meshes put all of this in `case.yaml` (`mesh.*`, `geometry.*`); the
mesh pipeline is described in [Meshes](../mesh.md).

## Prerequisites

| Step | Needs |
|---|---|
| Validate `case.yaml`, run the unit/regression tests | Python ≥ 3.9 with `plumetools` (native Windows works) |
| `./Allmesh --dict-only` | bash (WSL, Linux or Git Bash) + Python; no OpenFOAM |
| `./Allmesh --yes` | OpenFOAM v2512: `blockMesh`, `snappyHexMesh`, `checkMesh` |
| `./Allrun --no-solve` | a mesh in `constant/polyMesh`; `topoSet` optional |
| `./Allrun` | `dsmcFoam`, `dsmcInitialise`; `decomposePar`, `mpirun`, `reconstructPar` for parallel; `foamDictionary` for `-np` |

`./Allrun` reads the solver name from `application` in `system/controlDict`,
which is `dsmcFoam` (the standard solver), and derives the initialiser from it.
This case does **not** use the `plumeFieldInflow` library: `output.inflow_model`
is left at its default, `FreeStream`.

## Inputs

`case.yaml` is the only file to edit, and it is authoritative: the loader is
`plumetools.config.load_case_config`, which rejects unknown keys and
cross-checks the file against the mesh and `constant/dsmcProperties`. The full
schema is in the [configuration reference](../reference/configuration.md#legacy-source-flow-caseyaml-3d-inflow).

| Key | Value here | Units | Notes |
|---|---|---|---|
| `geometry.sphere_radius_m` | 0.5 | m | must equal `mesh.sphere_radius_m` (checked) |
| `geometry.plume_axis` / `polar_axis` | `[1,0,0]` / `[0,0,1]` | — | polar axis is LEGACY (SM-04) |
| `gas.gamma`, `gas.molar_mass_g_per_mol` | 1.4, 28.0134 | —, g/mol | **N₂**, while the solver simulates argon (SM-03, frozen) |
| `gas.species_name` | `Ar` | — | must appear in `dsmcProperties` `typeIdList` (checked); suffix of the `0/` field |
| `stagnation.p0_pa` | 3275009.71275 | Pa | LEGACY 475 psi; keep every digit (SM-01) |
| `stagnation.T0_K` | 300 | K | used by velocity and temperature |
| `stagnation.throat_radius_m` | 0.0041275 | m | must be smaller than the sphere radius (checked) |
| `angular.*` | 0.41, `legacy_separable`, 1.4, 500 | — | LEGACY angular model (SM-04, SM-05, SM-07) |
| `legacy.rhoN_T0_K` | 300 | K | density's own T₀ (SM-02); `null` = use `stagnation.T0_K` |
| `output.dialect` | `standard` | — | `standard` (OpenFOAM `dsmcFoam`) or `mnf` (`dsmcFoam+`); they read `0/boundaryT` as different field types |
| `mesh.type` | `snappy_hex_sphere` | — | or `block_mesh_ogrid`; see [Meshes](../mesh.md#choosing-between-the-two) |
| `mesh.outer_patch_type` | `wall` | — | **physics compromise**: the outer boundary reflects instead of absorbing ([why](../solver-compatibility.md)) |
| `mesh.background_cell_size_m`, `refinement_level` | 0.25, 4 | m, — | surface cell = 0.25 / 2⁴ = 0.015625 m |

!!! warning "Every `LEGACY` value changes results when flipped"
    The `legacy:` switches and the values annotated `LEGACY` reproduce the
    pre-refactor script bit for bit. Change one deliberately, record the
    before/after inflow statistics, and expect `pytest tests/regression` to fail
    from then on — that is what it is for.

## Cheap checks (no OpenFOAM)

All from the repository root unless stated.

=== "Native Windows or any Python"

    ```bash
    python -c "from plumetools.config import load_case_config; print(load_case_config('cases/3d-inflow').mesh)"
    ```

    Loading warns that `cylinder` and `plate` appear in `output.patches` but not
    in this case (AD-03, frozen via `legacy.emit_absent_patches`). Any other
    error names the key it rejects.

=== "WSL / Linux shell"

    ```bash
    cd cases/3d-inflow
    ./Allmesh --dict-only
    ```

    Writes `system/blockMeshDict`, `system/snappyHexMeshDict` and
    `system/meshQualityDict` (all gitignored, generated from `case.yaml`) and
    prints the background division and surface cell size. It runs no OpenFOAM
    command and touches neither `constant/` nor `0/`.

The regression test is the strongest available check, but it needs the
historical mesh (next section); without it, `pytest tests/regression` skips.
With a *generated* mesh in place it fails, because the golden is bound to the
Pointwise mesh.

## Get a mesh

`constant/polyMesh` is gitignored for this case, so a fresh checkout has none,
and `./Allrun` stops with a message saying so. Two ways out:

**Restore the historical Pointwise export** (99 240 tetrahedra, 2044 triangular
inflow faces — the mesh the regression golden is bound to). The tag the README
and `Allmesh` name, `pointwise-mesh-3d-inflow`, is **not present** in this
repository or on the remote. The mesh is still in history, in the parent of the
commit that stopped tracking it:

```bash
git checkout f83233f^ -- cases/3d-inflow/constant/polyMesh
```

That mesh declares its inflow patch `type Unspecified;`, a Pointwise export
artefact rather than an OpenFOAM patch type ([HA-07](#unresolved-settings-ha-06-and-ha-07)). Whether `dsmcFoam` v2512 accepts
it has not been recorded.

**Or generate one** from `case.yaml`:

!!! danger "`./Allmesh --yes` overwrites `constant/polyMesh`"
    It also deletes `constant/polyMesh/sets/` and the snappyHexMesh bookkeeping,
    which address the old mesh by index. A generated mesh is a *different* mesh
    (5452 polygonal inflow faces at the committed settings) and the regression
    golden will no longer match.

```bash
cd cases/3d-inflow
./Allmesh --yes       # blockMesh, snappyHexMesh -overwrite, checkMesh, verify
```

The last step is `python -m plumetools.verify_mesh .`, which checks the built
mesh against the geometry `case.yaml` declares and exits non-zero on a mismatch.

## Run

```bash
cd cases/3d-inflow
./Allrun --no-solve      # topoSet (if present), then runInflow.py -> 0/; no solver
./Allrun                 # decomposition from system/decomposeParDict
./Allrun --serial        # one process
./Allrun -np 4           # override the subdomain count (needs foamDictionary)
```

The order inside `Allrun` is load-bearing: `dsmcInitialise` **creates** `0/`,
then `runInflow.py` **overwrites** three files in it. Reversed, the inflow is
silently discarded. Every step's exit status is checked and the first FOAM
error is printed. A parallel run always ends with `reconstructPar`.

Before launching, `Allrun` refuses a standard-dialect run in which any boundary
other than the inflow has geometric type `patch`: standard `FreeStream` injects
on every such patch, and the run would abort with *Zero boundary temperature
detected*.

There is no restart logic in this case's `Allrun`; it does not use
`plumetools.foamio.checkpoints`. To rerun, clean first.

!!! danger "Resetting the case"
    `./Allclean` removes `0/`, every numeric time directory, `processor*`,
    `log.*`, `postProcessing`, `boundaries`, `fieldMeasurements`. Preview first
    with `./Allclean --dry-run`. It never touches `constant/polyMesh`;
    `--mesh` additionally removes the generated `system/blockMeshDict`.

## Outputs

| Path | Written by | Contents |
|---|---|---|
| `0/boundaryU` | `runInflow.py` | per-face inflow velocity, m/s, `volVectorField` |
| `0/boundaryT` | `runInflow.py` | per-face temperature, K; `volScalarField` in the `standard` dialect, `volVectorField (T 0 0)` in `mnf` |
| `0/boundaryNumberDensity_Ar` | `runInflow.py` | per-face number density, m⁻³ — **not read** by standard `FreeStream` |
| `0/{dsmcRhoN,fD,iDof,internalE,linearKE,momentum,q,rhoM,rhoN}` | `runInflow.py` (standard dialect) | zeroed measurement fields `dsmcFoam` needs at start |
| `log.*` | `Allrun`, `Allmesh` | one log per OpenFOAM application |
| time directories | `dsmcFoam` | the DSMC fields |

`runInflow.py` prints the face count, mean |U|, T and the density range. Under
the standard dialect it also prints the **area-weighted mean number density** to
put in `constant/dsmcProperties` → `FreeStreamCoeffs`, because standard
`FreeStream` takes one number density per species: the plume's angular
*direction* survives (per-face U and T), its angular *density profile* does not.
See [Solver compatibility](../solver-compatibility.md).

**Judging the output.** Success here means: the three boundary files exist and
are non-empty (Allrun checks this), the face count matches the inflow patch,
and |U| and T are what the configured — frozen — model implies. With the
Pointwise mesh restored, `pytest tests/regression` confirms bit-exact agreement
with the pre-refactor code; that proves *unchanged*, never *right*.

Post-processing is **not wired up**: the shipped `system/sampleDict` was copied
from the archived 2d-wedge case and samples from r = 1.005 m, outside this
case's 0.5 m inflow surface (RP-17). The `sampling:` block in `case.yaml` is
validated but nothing reads it.

## Unresolved settings: HA-06 and HA-07

Two settings are inconsistent with what this case's own descendants later did.
Neither has been changed, because confirming the consequence needs a DSMC run.

**HA-06 — `coordinateSystem dsmcAxisymmetric` on a 3D mesh.**
`constant/dsmcProperties` sets

```
coordinateSystem   dsmcAxisymmetric;
axisymmetricProperties { maxRadialWeightingFactor 1000; revolutionAxis "x"; polarAxis "z"; }
```

while the mesh is fully three-dimensional, the inflow model is φ-dependent and
there is a `sym` patch. It was inherited from the archived `2d-wedge` case by
directory copy (commit `51a75fe`); `iss_solar_panels_2`, same lineage, switched
to `dsmcCartesian` (commit `115d4f5`). If radial weighting is applied on a 3D
mesh, output from this case and from `iss_solar_panels` is affected. Verify
before trusting historical output.

**HA-07 — `type Unspecified;` on the inflow patch.** The Pointwise export's
`constant/polyMesh/boundary` declares
`inflow { type Unspecified; nFaces 2044; startFace 192243; }`. `Unspecified` is
a Pointwise V18.5R2 export artefact, not a registered OpenFOAM patch type;
`iss_solar_panels_2` and all three `Fnum_*` cases carry `type patch;` on the
same mesh. Either the historical solver tolerated it, or those cases ran from a
locally corrected `boundary` that was never committed. A generated mesh emits
`type patch;` and avoids the question.

## Frozen defects that matter most

Every `LEGACY` value in `case.yaml` reproduces pre-refactor behaviour bit for
bit; the full list and confidence classification are in
[The plume source-flow model](../source-flow-model.md).

| ID | Frozen behaviour |
|---|---|
| SM-01 | `p0_pa` is 475 psi (3.275 MPa). The original set 34 500 Pa (5 psi) and never read it |
| SM-02 | density uses `legacy.rhoN_T0_K` = 300 K while velocity and temperature use `stagnation.T0_K` (identical here; a factor of 10 apart in the resolution-study cases) |
| SM-03 | the model runs on N₂ (γ = 1.4, M = 28.0134) while the solver simulates argon: v_ℓ = 788.16 m/s instead of 558.89 m/s |
| SM-04 | θ is measured from +z, not from the +x plume axis, and the angular dependence is a separable product of two angles in different planes |
| SM-05 | the normalisation constant normalises a different angular function from the one the density uses |
| AD-03 | `cylinder` and `plate` entries are written into every field file; this mesh has neither (a warning on load) |

## Files in the case

| Path | |
|---|---|
| `case.yaml` | every physical and mesh parameter — the only file to edit |
| `Allrun`, `Allclean`, `Allmesh` | run, reset, regenerate the mesh (destructive, opt-in) |
| `runInflow.py` | writes the `0/` boundary fields (`python runInflow.py [case_dir]`) |
| `constant/dsmcProperties.mnf`, `system/dsmcInitialiseDict.mnf` | the dictionaries in the MNF fork's schema, kept beside the standard ones |
| `system/blockMeshDict.orig` | a pre-refactor dictionary kept for provenance; `Allmesh` writes `system/blockMeshDict` |
| `mesh/*.pw` | binary Pointwise project files; not reproducible without a licensed interactive session |
| `processInflowData.py.legacy` | the pre-refactor script, kept for provenance — **do not run** |
| `getMeshStats.py`, `printInflow.py` | superseded by `plumetools`; left untouched |
| `centroids.png`, `norms.png`, `sym_theta.png` | historical inflow-surface diagnostics |

The pre-refactor `clean.sh` (the unguarded `rm -r 0* …` of finding RP-01) is in
history: `git show 7997ccc:cases/3d-inflow/clean.sh`. The tag
`pre-refactor-baseline` that older text names for commit `7997ccc` is not
present in this repository.

## Resources

No measured run time or memory for this case is recorded in the repository.
Configured values: `numberOfSubdomains` in `system/decomposeParDict`; the
generated mesh at the committed settings has 32 272 cells (measured 2026-08-02,
recorded in `case.yaml`). The comment in `case.yaml` warns that a
`refinement_level` above 4 needs a smaller `deltaT` in `controlDict`.

## Variants and safe customisation

| Change | Where | Effect |
|---|---|---|
| Finer inflow surface | `mesh.refinement_level` (and `deltaT`) | surface cell halves per level; re-mesh |
| O-grid instead of snappy | `mesh.type: block_mesh_ogrid`, `mesh.projection: searchable_sphere` | hexahedral mesh, worse quality floor ([Meshes](../mesh.md)) |
| Absorbing outer boundary | `mesh.outer_patch_type: patch` | only valid with the MNF fork (`output.dialect: mnf`) or a patch-selected inflow |
| Correct a frozen defect | a `legacy:` switch | changes results; breaks the regression golden by design |

## Figure placeholder

!!! figure "Figure placeholder — inflow number density on the hemisphere"
    **Subject:** `0/boundaryNumberDensity_Ar` coloured on the `inflow` patch, viewed along −x.
    **Source run:** `./Allrun --no-solve` on the restored Pointwise mesh (`git checkout f83233f^ -- cases/3d-inflow/constant/polyMesh`).
    **Plotting command:** none verified in this repository; render with ParaView from `open.foam`.
    **Destination:** `docs/assets/figures/3d-inflow-boundary-density.png`
    **Suggested caption:** "Legacy source-flow number density on the 3d-inflow hemisphere (frozen SM-01…SM-05 defects included)."

## Related

- [The plume source-flow model](../source-flow-model.md) — equations, frozen defects, open questions
- [Meshes](../mesh.md) — snappy vs O-grid, verification, which polyMesh is tracked
- [Solver compatibility](../solver-compatibility.md) — why the outer boundary reflects
- API: [`plumetools.config`](../reference/api/legacy/config.md), [`plumetools.sourceflow`](../reference/api/legacy/sourceflow.md), [`plumetools.inflow`](../reference/api/legacy/inflow.md), [`plumetools.foamio.fields`](../reference/api/foamio/fields.md)
- [Troubleshooting](../guides/troubleshooting.md#3d-inflow-and-the-legacy-cases)

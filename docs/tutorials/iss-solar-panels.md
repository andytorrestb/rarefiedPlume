# ISS solar panels: three legacy impingement cases

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Cases** | `cases/iss_solar_panels`, `cases/iss_solar_panels_2`, `cases/iss_solar_panels_3` |
    | **Purpose** | Notional plume impingement on a solar panel (commit `115d4f5`, 2023-02-28) |
    | **Model** | The legacy source-flow model, as a **per-case copy** of the pre-refactor script `processInflowData.py` — not `plumetools` |
    | **Evidence** | None recorded. No test covers these cases; their committed images are historical output, not validated results |
    | **Status** | Labelled *active — under refactor* in [`cases/ARCHIVE.md`](../../cases/ARCHIVE.md), but **not runnable with the supported toolchain** ([blockers](#blockers)) |
    | **Start in** | the case directory |

## Purpose and lineage

These three cases were copied from `3d-inflow` to study a plume impinging on a
panel. They share its hemispherical inflow surface (radius 0.5 m, 2044
triangular faces from the same Pointwise export) and its analytical source-flow
model, but they were never converted to `plumetools`: each carries its own copy
of `processInflowData.py` with the inputs typed into the script.

The model, its conventions (plume axis +x, θ measured from +z) and its frozen
defects are described in [The plume source-flow model](../source-flow-model.md);
the defects SM-01 (pressure ignored) and SM-02 (density at 300 K) apply here in
their original, unfrozen form.

## Variants

Values read from the tracked files; nothing here was re-run.

| | `iss_solar_panels` | `iss_solar_panels_2` | `iss_solar_panels_3` |
|---|---|---|---|
| `coordinateSystem` (`dsmcProperties`) | `dsmcAxisymmetric` | `dsmcCartesian` | `dsmcCartesian` |
| `nEquivalentParticles` | 4e14 | 1e16 | 1e16 |
| Patches (`constant/polyMesh/boundary`) | `inflow` (Unspecified, 2044) · `sym` (symmetry, 2444) · `vacuum` (patch, 7986) | `inflow` (patch, 2044) · `vacuum` (patch, 17378) | `inflow` (Unspecified, 2044) · `panel` (Unspecified, 1686) · `vacuum` (patch, 15692) |
| `To` in `processInflowData.py` | 300 K (hard-coded in two functions) | 700 K | 700 K |
| `Po` set in the script | 6894 Pa ("5 psi") | 16894 Pa | 16894 Pa |
| `Po` actually used for density | 475 psi, hard-coded (SM-01) | same | same |
| `controlDict` | `dsmcFoam+`, `deltaT` 2e-5 s, `endTime` 0.05 s, `writeInterval` 5e-3 s | same | same |
| `numberOfSubdomains` | 8 | 8 | 8 |

`iss_solar_panels_3` also contains a nested copy of `iss_solar_panels_2`
(`cases/iss_solar_panels_3/iss_solar_panels_2/`), including stale
`__pycache__` files — a directory-copy artefact, not a fourth variant.

Two unresolved questions inherited from `3d-inflow` apply here: whether
`dsmcAxisymmetric` on a 3D mesh (variant 1) applied radial weighting, and whether
the `type Unspecified;` patches (variants 1 and 3) were accepted by the solver
that produced the history ([HA-06 and HA-07](3d-inflow.md#unresolved-settings-ha-06-and-ha-07)).

## Blockers

These are the exact reasons no working command sequence can be given today.

1. **The solver is the MNF fork.** Every `system/controlDict` names
   `application dsmcFoam+;`, and the dictionaries use the fork's schema. The
   supported toolchain is standard OpenFOAM v2512 `dsmcFoam`; no source for the
   fork build that produced the history has been located
   ([Environment](../environment.md)).
2. **There is no run script.** The cases ship `clean.sh` and a gnuplot
   `monitor`, but no `Allrun`; the order in which the original runs called
   `checkMesh`, `topoSet`, `dsmcInitialise+`, `processInflowData.py` and the
   solver is not recorded in these directories. (`3d-inflow`'s `Allrun`
   documents why that order matters.)
3. **`processInflowData.py` has prerequisites it does not create.** It reads
   `log.checkMesh` (so `checkMesh` must have run, with its output redirected
   there), `constant/polyMesh/sets/inflow` (written by `topoSet`), and the mesh
   files, from the current directory. It takes no arguments.
4. **The inputs are in the script.** Changing `To` or `Po` means editing
   Python, and `Po` is ignored anyway (SM-01).

!!! danger "Do not run `clean.sh`"
    Each case's `clean.sh` is `rm -r 0* boundaries fieldMeasurements postProcessing processor*`
    with no `cd` into its own directory. Run from anywhere else, it deletes
    matching paths **there** — from the repository root, that includes any
    top-level path starting with `0`. The `0*` glob would also take a `0.org`.
    This is finding RP-01; `3d-inflow/Allclean` is the guarded replacement.

## What you can do safely

Read-only inspection, from the repository root:

```bash
git show HEAD:cases/iss_solar_panels_3/constant/polyMesh/boundary   # patches and face counts
diff cases/iss_solar_panels_2/processInflowData.py cases/iss_solar_panels_3/processInflowData.py
```

Converting a variant to `plumetools` would mean writing a `case.yaml` for it and
adopting `3d-inflow`'s `Allrun`/`Allclean`/`runInflow.py`. That is application
work, not documentation, and has not been done; see
[Known gaps](../development/known-gaps.md). The `panel` patch footprint of
variant 3 is used as a geometry fixture by the `nozzleplume` unit tests
(`tests/unit/test_nozzleplume_geometry.py`), which is the only place these cases
are exercised.

## Outputs

Committed artefacts are historical: `centroids.png`, `norms.png`,
`sym_theta.png` (inflow-surface diagnostics drawn by the legacy scripts) and a
`VTK/nonOrthoFaces` export in variants 2 and 3. No solver output is tracked, and
no record says which solver build or inputs produced the images. Treat them as
provenance, not results.

## Resources

Configured only: 8 subdomains, `deltaT` 2e-5 s to `endTime` 0.05 s (2500 steps).
No measured run time or memory is recorded.

## Figure placeholder

!!! figure "Figure placeholder — panel impingement pressure"
    **Subject:** surface pressure on the `panel` patch of `iss_solar_panels_3`.
    **Source run:** none possible with the supported toolchain ([blockers](#blockers)); requires a `plumetools` conversion or an MNF `dsmcFoam+` build.
    **Plotting command:** none verified.
    **Destination:** `docs/assets/figures/iss-solar-panels-3-panel-pressure.png`
    **Suggested caption:** "Surface pressure on the notional panel, iss_solar_panels_3 (state the solver build and inputs actually used)."

## Related

- [3d-inflow tutorial](3d-inflow.md) — the converted parent case
- [Solar-panel particle-resolution study](particle-resolution.md) — the `Fnum_*` sweep on variant 2's mesh
- [Solver compatibility](../solver-compatibility.md) — `dsmcFoam` vs `dsmcFoam+`
- [Case catalog](../reference/case-catalog.md)

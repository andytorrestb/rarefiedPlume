# cai2012-hybrid: a staged 3-D CFD → DSMC free plume

!!! abstract "At a glance"

    | | |
    |---|---|
    | **Purpose** | Cai & Wang's circular argon source (D = 0.2 m, S₀ = 2, T₀ = 300 K) computed as a continuum precursor, a fixed continuum-breakdown surface, and DSMC beyond it — one-way coupling, full 3-D Cartesian |
    | **Model** | `rhoCentralFoam` (with the custom `libhybridPlumeThermo.so`) → breakdown diagnostics → interface selection and assessment → `dsmcFoam` + `plumeFieldInflow` on the CFD mesh minus the retained continuum (`plumetools.hybridplume`) |
    | **Evidence** | Analytical pieces unit-tested (half-range moments, flux signs, collisionless reference to 1e-13). Smoke runs **executed** on a laptop, all intentionally under-resolved, every interface verdict *unsuitable*. Verification and production profiles **prepared, not run**. Nothing converged or validated |
    | **Runs on** | Python for `setup` and every analysis step; OpenFOAM v2512 with **both** custom libraries for the solver stages |
    | **Start in** | the repository root; cases live under `cases/cai2012-hybrid/runs/` (gitignored) |

The method, every decision, and the measured smoke numbers are in
[`plumetools.hybridplume`, explained](../hybridplume.md). Configuration
attribution: Cai & Wang, *J. Spacecraft and Rockets* 49(1), 2012 — the physical
source is reproduced in 3-D; Cai's axisymmetric numerics are not.

## Geometry and coordinates

The source disk is at x = 0 on the +x axis, inside a box of
`domain.x_max_over_D` × 2·`domain.half_width_over_D` × 2·`domain.half_width_over_D`
diameters. The loader **refuses** `mesh.axisymmetric`, `mesh.wedge`,
`dsmc.axisymmetric`, `dsmc.radial_weighting` and `domain.symmetry`: this family
is full 3-D, with no wedge, symmetry plane, empty patch or radial weight. The
source patch is a staircase of whole cells, and its area error is reported and
never rescaled (+1.86 % at D/10 in the coarse smoke profile, from `setup`).

## Prerequisites

| Step | Needs |
|---|---|
| `setup` and every `python -m plumetools.hybridplume` analysis step | Python with `plumetools` (numpy, scipy, PyYAML; matplotlib for `report`'s plots) |
| `cfd` stage | `blockMesh`, `topoSet`, `createPatch`, `checkMesh`, `rhoCentralFoam`, and **`libhybridPlumeThermo.so`** |
| `dsmc` stage | `subsetMesh`, `checkMesh`, `dsmcInitialise`, `dsmcFoam`, and **`libplumeDsmcBoundaryModels.so`** |
| parallel | `decomposePar`, `mpirun`, `reconstructPar`, `foamDictionary` |

Build both libraries **in the environment that will run the case** (WSL and a
cluster each need their own build):

```bash
applications/thermophysicalModels/Allwmake    # libhybridPlumeThermo.so
applications/dsmcBoundaryModels/Allwmake      # libplumeDsmcBoundaryModels.so
```

`./Allsetup` builds only the second. `libhybridPlumeThermo.so` instantiates one
extra `psiThermo` — `hePsiThermo<pureMixture<tabulatedTransport<…hConstThermo<perfectGas>…>>>`
— because v2512 compiles tabulated transport for `rhoThermo` but not for
`psiThermo`, which `rhoCentralFoam` uses; the μ(T) and κ(T) tables are written
from the same gas object as the DSMC dictionaries, so the continuum transport
follows the VHS power law. See [the library contracts](../reference/configuration.md#custom-openfoam-libraries).

`Allrun` puts `$FOAM_APPBIN` first on `PATH`, because a `decomposePar` from
another OpenFOAM build on the owner's laptop aborts on every case.

## Inputs and profiles

`cases/cai2012-hybrid/hybrid.yaml` (schema `hybridplume/1`) holds every input,
each labelled [PAPER], [REPO], [BIRD], [ASSUMPTION] or [NUMERICAL]. A profile
is a partial overlay deep-merged onto `base`; `setup` writes the merged result
to `<case>/hybrid.resolved.yaml`, which every later step reads.

| `base` section | What it sets |
|---|---|
| `source` | `diameter_m`, `speed_ratio`, `temperature_K`, `knudsen`, `mean_free_path_convention` (`vhs`) |
| `gas` | argon VHS: `mass_kg`, `diameter_m`, `t_ref_K`, `omega`, `prandtl`, `provenance` |
| `domain`, `mesh` | extents in D; `core_cell_over_D`, `max_outer_expansion` (profiles add `core_x_over_D`, `core_half_over_D`, `outer_cells`) |
| `cfd` | `rhoCentralFoam`; background pressure ratio (a numerical floor, 0 < r < 0.1) and temperature, background-sensitivity ratios, viscous Δt factor, Courant cap, end time in transits, writes per transit, reconstruction limiter |
| `breakdown` | thresholds on Kn_GLL and B_exp (0.05), sensitivity factors (0.5, 1, 2), lip exclusion, overlap ceiling |
| `interface` | `mode` (`automatic` \| `manual` \| `none`), suitability `criteria` (returning mass fraction ≤ 0.05, subsonic-normal area ≤ 0.25, no outer-boundary contact, ≥ 1 overlap layer), `allow_unsuitable`, the `manual` surface (`sphere` \| `cylinder` \| `box`) |
| `dsmc` | backend `dsmcFoam`, collisions, VHS, parcels per densest cell, Courant and Δt/τ caps, startup and sample transits, writes per transit, seed, Lagrangian moments, crossing window |
| `resources` | `max_ranks` (never "all cores"), `max_memory_GB`, bytes per parcel/cell, seconds per parcel-step and per CFD cell-step (cost-model constants) |

Committed profiles (`default_profile: smoke-kn0p01`):

| Profile | Kn | Mesh core | Interface | Purpose | Where |
|---|---|---|---|---|---|
| `smoke-kn0p01` | 0.01 | D/20 | automatic | hybrid machinery, under-resolved | laptop |
| `smoke-kn0p01-coarse` | 0.01 | D/10 | automatic | the same, cheaper | laptop |
| `smoke-kn0p01-coarse-full` | 0.01 | D/10 | `none` (full DSMC) | overlap reference on the same mesh | laptop |
| `smoke-kn1em3` | 0.001 | D/20 | automatic | extension machinery | laptop |
| `smoke-collisionless` | 100 | D/20 | `none`, collisions off | DSMC vs the independent collisionless reference | laptop |
| `verification-collisionless` | 100 | D/40 | `none`, collisions off | r/D = 10 reachable | HPC (64 ranks) |
| `verification-kn0p1` | 0.1 | D/40 | `none` | Cai Kn 0.1, full DSMC | HPC (96 ranks) |
| `verification-kn0p01` | 0.01 | D/100 | automatic | Cai Kn 0.01, hybrid | HPC (192 ranks) |
| `production-kn1em3` | 0.001 | D/100 | automatic | extension | HPC (192 ranks) |
| `production-kn1em4` | 0.0001 | D/133 | automatic | extension | HPC (192 ranks) |

Schema details: [configuration reference](../reference/configuration.md#hybridyaml).

## Cheap check: `setup` only (no OpenFOAM)

```bash
python -m plumetools.hybridplume setup cases/cai2012-hybrid/runs/kn0p01c \
    --config cases/cai2012-hybrid/hybrid.yaml --profile smoke-kn0p01-coarse --np 2
```

It resolves the profile, estimates resources (CFD and DSMC cells, parcels,
memory, Δt, steps, runtime, Δx/λ), and writes `hybrid.resolved.yaml`,
`results/setup.json` and the `cfd/` OpenFOAM case — in about ten seconds on a
laptop, touching nothing outside the case directory. `--background-ratio R`
writes a background-floor sensitivity copy `cfd_bg<R>/` instead of `cfd/`. Read
the estimates before running anything: a dense case is never called affordable
because its mesh fits in memory, and a smaller weight does not repair
Δx/λ > 1.

## The staged run

Each stage is one `Allrun` call; `Allrun` runs OpenFOAM tools and the Python
steps write inputs or read outputs.

```bash
C=cases/cai2012-hybrid/runs/kn0p01c
cases/cai2012-hybrid/Allrun $C setup --profile smoke-kn0p01-coarse -np 4
cases/cai2012-hybrid/Allrun $C cfd -np 4         # blockMesh, topoSet, createPatch, checkMesh, check-mesh, rhoCentralFoam
cases/cai2012-hybrid/Allrun $C post-cfd          # cfd-report, breakdown, interface -> results/interface.json
cases/cai2012-hybrid/Allrun $C dsmc -np 4        # dsmc-setup, subsetMesh, checkMesh, dsmc-fields, dsmcInitialise, dsmcFoam
cases/cai2012-hybrid/Allrun $C post-dsmc         # dsmc-post, reference, report -> results/summary.json, *.png
cases/cai2012-hybrid/Allrun $C crossings -np 4   # optional: exact exits over a per-step window
```

`all --profile NAME` runs setup → cfd → post-cfd → dsmc → post-dsmc. Other
`Allrun` options: `--config FILE`, `--subdir cfd|cfd_bg…` (for the `cfd` stage),
`--allow-unsuitable`.

**Profile choice.** For a full-DSMC profile (`interface.mode: none`) the `cfd`
stage builds the mesh but solves no precursor, and `post-cfd` only records
`full_dsmc`. To compare a hybrid with full DSMC on one mesh, run
`smoke-kn0p01-coarse` and `smoke-kn0p01-coarse-full` into two case directories,
then:

```bash
python -m plumetools.hybridplume overlap cases/cai2012-hybrid/runs/kn0p01c \
    --reference cases/cai2012-hybrid/runs/kn0p01c_full
```

**Interface suitability.** `post-cfd` always records a verdict in
`results/interface.json`; an *unsuitable* verdict is a result, not a failure.
`dsmc-setup` then **refuses** to build the DSMC region unless the interface was
accepted with `--allow-unsuitable`, in which case every report is labelled a
machinery smoke test. Every committed smoke run is unsuitable (returning
fraction 0.25–0.68, large subsonic-normal fractions) — see the explanation's
§10.

### Resume behaviour

Every solver stage resumes from the newest time complete on every rank
(`python -m plumetools.foamio.checkpoints plan`); a stage at `endTime` is not
re-run; rerunning the same command after an interruption continues it.
**Nothing is deleted**: to start over, remove the case directory. Two cautions
recorded by the owner: the `dsmcFoam` restart insertion burst is unexplained
(see the explanation's limitations), and an `Allrun` should not be edited while
a stage is running — bash reads scripts incrementally.

The `crossings` stage extends a **finished** DSMC run by
`dsmc.crossing_window_steps` steps written every step (the main `controlDict`
is kept as `system/controlDict.main`) and then counts every exit exactly.

### On a cluster

```bash
python -m plumetools.hybridplume hpc cases/cai2012-hybrid/runs/<case> --scheduler pbs --np 192
```

writes `<case>/hpc/run_<stage>.pbs` (and `results/hpc_scripts.json`), sized from
the resource estimates: walltime = 3 × estimate, at least 1 h, at most 120 h;
queue `longq`, `nodes=1:ppn=192` copied from the Kn extension's `runCase.sh`.
The scripts source `cases/cai2012-kn-extension/hpc_env.sh` — the owner's
cluster environment — unless `HYBRID_HPC_ENV` names another. `qsub` them in
stage order; a walltime kill is resumed by resubmitting the same script.
`--scheduler slurm` (the CLI default) writes scripts with `__SET_ME__`
partition and account placeholders that make `sbatch` fail until set: there
are no SLURM settings in this repository. `--scheduler none` writes plain
scripts. Compile both libraries on the cluster first.

## Outputs

All under `<case>/results/` (JSON written by the step named):

| File | Step | Contents |
|---|---|---|
| `setup.json` (`setup_<subdir>.json` for a sensitivity copy) | `setup` | profile, mesh, source discretisation and area error, resource estimates, ranks |
| `mesh_check_cfd.json` | `check-mesh` | built source patch against the predicted staircase (count, area) |
| `cfd_report_cfd.json` | `cfd-report` | precursor budgets: mass residual, exit-plane numerical inflow, background sensitivity |
| `breakdown.json` | `breakdown` | Kn_ρ, Kn_T, Kn_U, Kn_GLL, B_exp on the CFD cells; retained volume per threshold factor |
| `interface.json` | `interface` | the surface, its assessment against `interface.criteria`, the verdict (`full_dsmc` for mode none) |
| `dsmc_setup.json`, `dsmc_fields.json` | `dsmc-setup`, `dsmc-fields` | DSMC run plan; mapped inflow per face and its flux error |
| `dsmc_post.json` | `dsmc-post` | moments, measured fluxes, budgets, comparisons |
| `crossing_window.json` | `crossings` | exact exits over the per-step window |
| `reference.json` | `reference` | the independent collisionless reference at the report stations |
| `overlap.json` | `overlap` | hybrid against full DSMC, cell by cell |
| `hpc_scripts.json` | `hpc` | the written batch scripts and their sizing |
| `summary.json`, `*.png` | `report` | every step's result collected, and the figures |

Number density is n in m⁻³ and normalised by the source n_e where stated;
temperatures in the DSMC moments distinguish parallel and perpendicular
components (the collisionless smoke run shows T∥/T⊥ from 0.94 to 3.75 on the
axis), so a scalar temperature is not the whole state.

**Judging the output.** Check `interface.json` first: an unsuitable verdict
means the one-way hybrid is not a valid approximation there, whatever the
downstream fields look like. Then the budgets (mass residual, numerical inflow
at x = 0 — the background floor leaks 1–4 % of the source flow inward in the
smoke runs), the measured against expected gross injection, and the noise
floor quoted beside each comparison. A smoke result is mechanics, not physics.

## Resources

Estimates come from `setup` using the `resources` constants: 4 × 10⁻⁶ s per
parcel-step per rank by default, measured on the owner's laptop (WSL) and to be
re-measured on a cluster; profiles cap ranks (`max_ranks`) and memory
(`max_memory_GB`, 8 GB for laptop smoke profiles, up to 1200 GB for
production). The laptop runs were made on a 32 GB host with a 15 GB WSL VM,
never on all cores.

## Variants and safe customisation

| Want | Do |
|---|---|
| A different Kn | a new profile overlaying `source.knudsen` and a mesh fine enough for it; `setup` reports Δx/λ |
| A manual surface | `interface.mode: manual` and `interface.manual` (`sphere`/`cylinder`/`box`, size and centre in D) |
| Background-floor sensitivity | `setup … --background-ratio 1e-5`, then `Allrun $C cfd --subdir cfd_bg1e-05` |
| Fewer ranks | `-np N` on each stage; never more than the machine has |

## Figure placeholder

!!! figure "Figure placeholder — breakdown surface on the meridional plane"
    **Subject:** Kn_GLL and B_exp on the y = 0 plane of the precursor, with the selected interface outlined.
    **Source run:** `smoke-kn0p01` (D/20) after `cfd` and `post-cfd`.
    **Plotting command:** `python -m plumetools.hybridplume report <case>` draws the summary figures (verified in source; which panels it draws depends on the steps that ran).
    **Destination:** `docs/assets/figures/hybrid-breakdown-surface.png`
    **Suggested caption:** "Continuum-breakdown indicators and the fixed interface, Kn = 0.01 smoke (intentionally under-resolved; verdict unsuitable)."

## Related

- [`plumetools.hybridplume`, explained](../hybridplume.md) — backend decision, physics, breakdown, interface criteria, transfer, measurements, verification record, limitations
- [plumeFieldInflow](../plume-field-inflow.md) — the DSMC inflow the interface uses
- [Restart and recovery](../guides/restart-and-recovery.md) · [PBS jobs](../guides/pbs.md)
- API: [`hybridplume.workflow`](../reference/api/hybridplume/workflow.md), [`config`](../reference/api/hybridplume/config.md), [`physics`](../reference/api/hybridplume/physics.md), [`interface`](../reference/api/hybridplume/interface.md), [`halfrange`](../reference/api/hybridplume/halfrange.md)

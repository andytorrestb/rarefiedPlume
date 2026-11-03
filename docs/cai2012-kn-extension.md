# The Cai 2012 Knudsen extension

Cai, C., and Wang, L., *Numerical Validations for a Set of Collisionless Rocket
Plume Solutions*, Journal of Spacecraft and Rockets, Vol. 49, No. 1, 2012.

This document says what the extension is, what it reuses, what it adds, what
every compared quantity means, and what is not yet known. Running instructions
are in the [Kn-extension tutorial](tutorials/cai2012-kn-extension.md);
the parent study is documented in [`docs/cai2012-case.md`](cai2012-case.md).

> **Current setup (2026-10-26): seven-decade full suite prepared, not run.**
> `full_decades` generates Kn_e = 100, 10, 1, .1, .01, .001, .0001 into
> `Cases_full_suite/` with `manifest-full-suite.yaml`, preserving prior trees.
> The dictionary-only smoke passed 84 checks for seven cases and confirmed
> idempotence and scaling registration. No mesh, DSMC or scheduler job was run.
> Four cases are pre-run candidates; .01/.001/.0001 remain non-production.
> The wrong very-low-Kn direct-exit inflow surface is explicitly retained for
> future work. Current commands are in the study README. Historical validation
> records below remain labeled with their original date.

---

## 1. Scope

| Fixed | Value | Source |
|---|---|---|
| Nozzle exit diameter `D` | 0.2 m | [PAPER] |
| Exit speed ratio `S0 = U0/√(2RT0)` | 2 | [PAPER] |
| Exit temperature `T0` | 300 K | [ASSUMPTION], as in `cases/cai2012` |
| Exit condition | uniform drifting Maxwellian, `plumeFieldInflow` | [PAPER] / repository |
| Gas | argon, VHS `d = 4.17e-10 m`, `ω = 0.74`, `T_ref = 273 K` | repository model |
| Mean-free-path convention | VHS, `λ = (T0/T_ref)^(ω−½)/(√2 π d² n)` | `cases/cai2012` |
| Collisions | on, at every Kn (NTC, standard `dsmcFoam` v2512) | [PAPER] |
| Geometry | full 3-D Cartesian, `0 ≤ X/D ≤ 10`, `|Y|,|Z|/D ≤ 10` | `cases/cai2012` |

**Varied:** Kn_e in {100, 10, 1, .1, .01, .001, .0001}, on the
nozzle diameter. There is one particle multiplier, np1x. The configured exit-cell
parcel target is 2000 for Kn>=1 and 40 for Kn<=.1; it is verified per case,
not treated as a physical invariant. Density follows Kn -> lambda=Kn*D -> n0
using the existing VHS convention. No number density is entered by hand.

**Policy:** `meta.knudsen_suite: full_decades` explicitly permits all seven in
generation, manifests, case-level guards and direct-exit post-processing.
Unmarked historical studies, manifests and cases retain the four-Kn policy and
paper-Kn exclusions. The template's physical settings remain inherited.
External-study results are not imported. Every active case is registered in
scaling through the study's manifest/cases_dir, with missing data reported.

**Source limitation:** at Kn<=.001 the direct-exit inflow surface is a known
limitation deferred by the user. Metadata, the manifest, KNOWN-LIMITATIONS,
HPC logs and scaling suitability retain the warning. This setup does not replace
it with a downstream source, or weaken mesh/time-step requirements.

**Not in scope:** axisymmetric or planar geometry, impingement plates,
source-surface sensitivity.

---

## 2. What was reused, and what was added

The audit of `plumetools/cai2012`, `cases/cai2012` and `cases/cai2012-health`
found almost everything already present. The project's standing rule is that
the paper-specific packages and their case families stay untouched as the
record of what was run, so **nothing in them was edited**: the extension
imports them, and its new code is a separate package, `plumetools/cai2012ext/`.

| Need | Existing capability | Use here |
|---|---|---|
| Study schema, Kn × particle levels | `plumetools.cai2012.study` (`load_study`, `apply_case_parameters`, `manifest_entry`) | **reused unchanged**; `study.yaml` is in that schema |
| Kn → n0, U0, exit state | `plumetools.cai2012.gas`, `inflow` | reused |
| Mesh plan, budget coarsening | `plumetools.cai2012.mesh.plan` | reused; an unbudgeted plan is also computed to report the criterion's cell count |
| Weight and time step | `inflow.derive_run_settings` (`weight = n0 V_exit / ppc / multiplier`) | reused, at multiplier 1 |
| Pre-run resolution checks | `plumetools.nozzleplume.checks` (names, definitions) and `nozzleplume.config.ChecksConfig` (default thresholds) | the same checks, evaluated on the cai2012 plan and run settings (§5) |
| Post-run resolution diagnostics | `plumetools.rarefaction.run.diagnose`; `nozzleplume.diagnostics.plot_centreline`; R-4D `viz.yaml`'s diagnostic entries | reused unchanged, fed from a cai2012 `case.yaml` (§7a) |
| Cai startup (10⁴ t0 at Kn = 0.01) | `dsmc.transient_basis: cai`, `resolution.reference_knudsen` | reused, switched on |
| Dictionaries, meshing, running | `cases/cai2012/baseCase` scripts | copied **byte-identical** (test-pinned) |
| Checks | `plumetools.cai2012.checks` | reused, and read as production criteria |
| Tree validation | `plumetools.cai2012.verify.verify_tree` | reused, plus independent checks |
| Centreline n, U, T and n contours | `plumetools.cai2012.post` via `postProcess.py` | reused unchanged |
| Off-axis analytical n, **u (3-D)**, T | `plumetools.cai2012.analytical.moments` | reused — **no new analytical formula was needed** |
| Occupancy audit, per-frame centreline error | `plumetools.cai2012.audit`, health `Allaudit.py` | reused byte-identical |
| Cost model, completed-run test, wall-clock record | `plumetools.cai2012.health`, health `AllrunCases` | reused, per case, by `hpc.py` from each case's `runCase.sh` |
| Averaged-field imagery, pinned ranges, sampling window | `plumetools.viz`, health `viz.yaml` | reused; window and n/n0 now from the manifest |

| Added (`plumetools/cai2012ext/`) | Why it was missing |
|---|---|
| `policy.py` | fixed matrix, exclusions, inheritance check, the manifest's Kn gate |
| `provenance.py` | independent Kn recomputation; ρ0, p0, startup in three units, code version |
| `feasibility.py` | production verdicts from the nozzle-plume checks, criterion cell count, cost/disk |
| `hpc.py` | what a case's PBS job does around `./Allrun`: the excluded-Kn guard, the `NON-PRODUCTION` label, the finished test, the wall clock recorded under a file lock |
| `diagnostics.py` | the adapter that feeds `plumetools.rarefaction` from a cai2012 `case.yaml` and adds the exit-plane check |
| `generation.py` | the parent generator deletes and re-copies un-run cases; this one never touches a matching case |
| `verification.py` | recomputes what the extension promises (§10) |
| `comparison.py` | angular profiles, velocity projections, pressure, masked metrics |
| `convergence.py` | far-field stationarity — the parent has only the per-frame centreline error |
| `figures.py`, `summary.py`, `vizspec.py` | fixed-scale figures, study table, per-Kn viz specs |

---

## 3. The seven active configurations

Derived from the current `manifest-full-suite.yaml`; generated inputs are under
`Cases_full_suite/`. These are configuration checks, not simulation results.

| Case | Kn_e | n_e [m^-3] | Exit-cell parcel target | Pre-run status |
|---|---:|---:|---:|---|
| `Kn100_np1x` | 100 | 6.62008e+16 | 2000 | production-candidate |
| `Kn10_np1x` | 10 | 6.62008e+17 | 2000 | production-candidate |
| `Kn1_np1x` | 1 | 6.62008e+18 | 2000 | production-candidate |
| `Kn0p1_np1x` | 0.1 | 6.62008e+19 | 40 | production-candidate |
| `Kn0p01_np1x` | 0.01 | 6.62008e+20 | 40 | unsuitable-for-production |
| `Kn0p001_np1x` | 0.001 | 6.62008e+21 | 40 | unsuitable-for-production |
| `Kn0p0001_np1x` | 0.0001 | 6.62008e+22 | 40 | unsuitable-for-production |

`nEquivalentParticles = n_e V_exit / configured_target`, derived separately
for each case. The reference startup and averaging windows remain shared.

---
## 4. Settings inherited from `cases/cai2012-health`

The template is `cases/cai2012-health/baseCase/case.yaml` with exactly five keys
changed. `policy.inheritance_report` compares the two templates key by key at
every generation; 62 keys are inherited, and a difference not declared below
stops generation.

| Setting | Value | Source | Status |
|---|---|---|---|
| `resolution.target_particles_per_cell` | **40** | health `study.yaml`, weights `ppc040` | override of the health template's 20 — its generator writes the row value per case |
| `dsmc.sampling_domain_transits` | **4.5** | health `study.yaml`, samplings `s4p5` | override of the template's 1.5, same reason |
| `dsmc.transient_basis` | **`cai`** | `cases/cai2012` template default | **deliberate override** of health's `transits` (2 transits), which that family documents as far-field-limiting |
| `exit.knudsen` | per case (template: 10) | — | health pins 100; the full profile now includes it |
| `dsmc.n_subdomains` | **192** | — | override of the health template's 12 (its development machine's cores): one HPC node of 192 cores, each case its own PBS job; `./Allrun` reads it from `system/decomposeParDict`, `simple` `(4 6 8)` |
| `dsmc.write_interval_s` | 3.536e-4 s | health template | inherited; rounded to whole steps: 258 (3.5420e-4 s) at Kn ≥ 1, 496 (3.5374e-4 s) at Kn = 1e-3, 483 (3.5365e-4 s) at Kn = 1e-4 |
| `dsmc.courant_target` | 0.2 | health template | inherited; `Δt` recomputed per mesh |
| `mesh.*` (`target_cell_over_mfp` 1, `max_core_cell_over_D` 0.05, core 2 D × 1.5 D, 40/34 outer cells, expansion cap 40, `max_cells` 4e6) | | health template | inherited; the mesh is re-planned per Kn |
| `checks.*` | | health template | inherited; read as production criteria |
| `post.*` (centreline tube 0.10 D, 200 bins, `xz` plane) | | health template | inherited |
| `run_time_multiplier`, `output_frequency_multiplier` | 1, 1 | — | pinned: the parent study's 3 and 4 would multiply the window and cadence |

The health family's **generated** `ppc040/s4p5/case.yaml` is not available
locally (its `Cases/` tree is not in this checkout), so the comparison is
against its template plus its `study.yaml` rows, which is what generation
produces. The health family's `Kn = 100` mesh, `Δt` and weight are **not**
copied: all three are recomputed for each Kn.

---

## 5. Resolution and feasibility, before meshing

| Kn_e | Planned cells | cell/lambda_e | dt/tau_e | Pre-run status |
|---:|---:|---:|---:|---|
| 100 | 4,000,000 | 0.001 | 5.4758e-05 | production-candidate |
| 10 | 4,000,000 | 0.01 | 0.00054758 | production-candidate |
| 1 | 4,000,000 | 0.1 | 0.0054758 | production-candidate |
| 0.1 | 1,310,720 | 0.5 | 0.027379 | production-candidate |
| 0.01 | 3,726,000 | 2.6941 | 0.14602 | unsuitable-for-production |
| 0.001 | 3,918,213 | 26.17 | 1.4223 | unsuitable-for-production |
| 0.0001 | 3,726,000 | 269.51 | 14.602 | unsuitable-for-production |

`τ0 = λ0/c̄` is the exit mean collision time, the densest state in the domain.
Courant is a molecule at `U0 + 3σ` against the smallest cell.

`feasibility.assess` reads each case against the **nozzle-plume family's
pre-run checks** (`plumetools.nozzleplume.checks`) — the same names and
definitions, so an extension case and an R-4D case are judged by one set of
resolution rules — plus two of this family's:

| check | threshold | source |
|---|---|---|
| `dsmc_checks` (no hard failure: Courant > 1, nozzle faces and area) | — | `plumetools.cai2012.checks` |
| `mesh_budget` (mesh not coarsened to fit) | `mesh.max_cells` (4e6) | template |
| `cell_over_mean_free_path` | ≤ `checks.max_cell_over_mfp` (1) | template; Cai's `dx/λ0 = 1` |
| `time_step_collision_time` (`Δt/τ0`) | ≤ 1 | nozzle-plume `ChecksConfig.max_time_step_over_collision_time` default |
| `time_step_courant` | ≤ `checks.warn_courant` (0.5) | template |
| `transient_length` | averaging starts after one flight of the jet across the domain (no plate) | nozzle-plume `RunSettings.settling_time_s` definition |
| `internal_energy_exchange` | no internal modes frozen by the collision model | nozzle-plume check; argon, 0 dof |
| `particles_per_cell` (exit) | ≥ `checks.min_particles_per_cell` (5) | template |
| `parcel_count` | ≤ 3e7 | nozzle-plume `ChecksConfig.warn_parcels` default |

The two thresholds the cai2012 `checks` section has no key for are imported
from `plumetools.nozzleplume.config`, not retyped.

**Kn = 1 and 10** pass all of them (candidates). **Kn = 1e-3 and 1e-4** fail
`cell_over_mean_free_path`, `mesh_budget` and `time_step_collision_time`. The
physical inputs of those cases are still exact — no Kn was raised, the nozzle
was not enlarged, the gas was not changed. The mesh was coarsened by the
existing budget rule and the verdict says so; **more parcels would change none
of the numerical criteria** (asserted in
`tests/unit/test_cai2012ext_policy.py`). Startup passes `transient_length` 17.7
times over in all seven.

---

## 6. Startup, averaging and writes

**Startup** is Cai's: 10⁴ reference collision times, with the reference the
**Kn = 0.01 exit state** for every case —
`t_ref = λ_ref / c̄ = 0.002 m / 398.9 m/s = 5.0143e-6 s`, so sampling starts at
`5.01435e-2 s`, identical in all seven active cases. `validate_cases.py` recomputes
this from `λ_ref = 0.01 D` and `c̄ = √(8kT0/πm)` without the gas module.

| | Kn = 10 | Kn = 1 | Kn = 1e-3 | Kn = 1e-4 |
|---|---|---|---|---|
| reference collision times | 10⁴ | 10⁴ | 10⁴ | 10⁴ |
| the case's own exit collision times | 10 | 10² | 10⁵ | 10⁶ |
| domain transits | 17.7 | 17.7 | 17.7 | 17.7 |

The reference is held deliberately — it is what the paper specifies and what
makes the seven startups the same instant. In a supersonic plume into vacuum
the time to steady state is set by transits (the slow tail of the exit
Maxwellian filling the far field, `docs/cai2012-case.md` §5.1), and 17.7
transits is nine times the health family's 2. That it is only 10 local
collision times at Kn = 10 is irrelevant to filling there, because the flow is
nearly collisionless.

**Averaging** is separate: `dsmc.sampling_domain_transits = 4.5` after the
startup, i.e. `1.273e-2 s`. The end time is snapped up to a whole write
interval, so the achieved window is 4.53–4.56 transits (recorded per case). No
study-wide multiplier touches either.

**Writes** keep the health cadence, 1/8 transit. Through a 17.7-transit startup
that is **178 writes per case, 141 of them transient** — no `*Mean` field, read
by nothing. `controlDict` has one `writeInterval`, so the cadence cannot differ
between startup and sampling without changing `plumetools.cai2012.dictionaries`,
which is frozen. The consequence is the disk estimate in the README (~212 GB
per rarefied case, ~2.7 TB for all four): nothing is pruned, so every write's
fields and parcel cloud are kept under `processor*/` and again reconstructed.
A two-cadence `controlDict` would be the real fix (§11).

Old results are not touched: nothing in this family reads or rewrites
`cases/cai2012` or `cases/cai2012-health` outputs. If results from those
families are ever supplied for comparison, their sampling histories are
incompatible with this one — 2 transits of startup and 1.5 (or 0.5/4.5) of
averaging, against 17.7 and 4.5 here — and must not be pooled.

---

## 7. What is compared, and how

Everything is computed from the **time-averaged** fields (`rhoNMean`,
`rhoMMean`, `momentumMean`, `linearKEMean`, `dsmcRhoNMean`), and every
analytical value comes from `plumetools.cai2012.analytical` — the
free-molecular disk integral already in the repository, re-derived there from
first principles. **The reference throughout is Cai's collisionless analytical
solution.** No external-study DSMC fields are imported; Cai's reported
errors at Kn = 100 / 0.1 / 0.01 are not used as tolerances, and the study table
has no column for them.

**Coordinates.** `x` along the axis from the exit centre; `r = √(y²+z²)`
cylindrical; `R = √(x²+r²)` spherical; `θ = atan2(r, x)` from the axis.

**Normalisation.** `n/n0`; velocities as `u√β0` (local speed ratios);
`T/T0`; `p/p0`.

**Matched estimators.** A DSMC value is an average over cells; the analytical
solution is evaluated at the **same cell centres** and averaged with the **same
weights** (volume for `n`, `p`, `F`; molecule count for `u`, `T`), and errors
are taken against that. The exact curve on the line, sphere or plane is drawn
alongside. This is `post.centerline`'s method extended to every output.

| # | Output | Definition | Sampling region | Cai |
|---|---|---|---|---|
| 1 | Centreline `n/n0`, `U_x√β0`, `T/T0` | `post.centerline`, unchanged; plus `p/p0` | cells with `r ≤ 0.10 D`, 200 bins in `0 ≤ X/D ≤ 10` | Eqs. 18–21, Figs. 19–21 (as cited in `docs/cai2012-case.md`; closed forms re-derived there) |
| 2 | `n/n0` contours | `post.density_plane` / `analytical_plane`, levels 0.1, 0.01, 0.001 | the two cell layers straddling `y = 0`, `|Z/D| ≤ 5` | Eq. 5, Figs. 6–8 (as above) |
| 3 | Angular `n/n0` at r/D = 10 | shell average per 2.5° polar bin | cells whose radial extent reaches `9.5 ≤ R/D ≤ 10.5`, all azimuths, each shared among the bins its angular footprint overlaps | **not verified** (§8) |
| 4 | Angular mass flux at r/D = 10 | `F(θ) = n u_R R² / (Γ0 π R0²)` per steradian; `Γ0` the exit one-way Maxwellian flux | as 3 | **not verified** |
| 5 | Off-axis `u_x√β0`, `u_r√β0` | `u_r = (y u_y + z u_z)/r`, cylindrical, positive away from the axis | centre-plane cells, as 2 | **not verified** |
| 6 | Off-axis `T/T0`, `p/p0` | `3RT = ⟨v²⟩ − |u|²` (local translational); `p = nkT`, `p/p0 = (n/n0)(T/T0)` | as 5 | **not verified** |
| 7 | Per-case metrics, study table | §8; `results/comparison.yaml`, `results/study-table.csv` | all of the above | — |
| 8 | Knudsen comparison | the active Kn cases on one figure (`knudsen-comparison.png`): centreline `n`, `U`, `T`, angular `n`, `F`, and the centreline `Δmax/λ` from §7a beside them. Normalised, the collisionless curve is the same for every Kn, so the spread is the collisional departure | as 1, 3, 4 | — |

**The 3-D velocity in a meridional plot.** Each cell's Cartesian velocity is
projected: `u_x` (axial), `u_r` (cylindrical radial), `u_φ = (y u_z − z u_y)/r`
(swirl — zero for this flow, reported as a symmetry diagnostic), and the
spherical `u_R = cos θ u_x + sin θ u_r`, `u_θ = −sin θ u_x + cos θ u_r`. `u_r`
and `u_R` are **different quantities**: the plane shows `u_r`, the sphere uses
`u_R` because it carries the mass through it. On the axis the projections are
0/0; the axisymmetric limit `u_r = u_φ = 0` is used (verified analytically:
`u_r` is even in the lateral coordinate and vanishes on the axis).

**Pressure** is the scalar kinetic pressure, a third of the trace of the
pressure tensor, from the **local** temperature. Never `n k T0`. In a
collisionless expansion the tensor is anisotropic; this is its mean, and the
same definition is applied to the analytical solution.

**Mass flux and the area factor.** Integrating `F` over the hemisphere needs
`dΩ = 2π sin θ dθ`. Collisionless, every molecule leaves the disk forward and
crosses the r/D = 10 sphere once, so `∫ F 2π sin θ dθ = 1` exactly. The
analytical moments reproduce this to 8e-6 at 0.25° bins (1.0008 at the default
2.5°, the midpoint-rule error) — an independent check of `Γ0`, the `R²` and the
area factor (`tests/unit/test_cai2012ext_comparison.py`). The DSMC integral is
reported beside it; collisions that scatter molecules back through `x = 0`
lower it, and the meshed nozzle area (+0.59% at Kn ≥ 1) raises it.

**Near the nozzle lip** (`x = 0, r = R0`) the analytical density is
discontinuous and a cell-centre value is not a cell average. Points within
0.1 D of the lip are excluded from the plane metrics and counted
(`plane_excluded_nozzle_lip`). In the exit plane itself (`x ≤ 0`) the disk
integral is singular; every cell centre has `x > 0`, and analytical grid
points there are NaN rather than an error.

**Coverage of the sphere.** At r/D = 10 the graded mesh's outer cells are
0.4–0.5 D — as wide as a 2.5° bin there — so each cell is shared among the bins
its footprint `V^(1/3)/R` overlaps, and the matched analytical value goes
through the same sharing. On every production mesh all 36 bins are covered
(≥ 228 cells per bin at Kn ≥ 1, ≥ 783 on the dense meshes). Empty bins, if any,
are counted, and the collisionless integral is also reported over the covered
bins only.

### 7a. Resolution diagnostics, after the run

`Alldiagnose.py` (step 2 of `./AllpostCases`) runs the **nozzle-plume family's
rarefaction and mesh-resolution diagnostics** —
[`plumetools.rarefaction`](rarefaction-diagnostics.md), the pipeline R-4D's
`./Allpost` runs — on each case. `plumetools/cai2012ext/diagnostics.py` is the
adapter: it takes the gas, the mean-free-path convention and `D` from the
cai2012 `case.yaml`, the thresholds from the case's `diagnostics.yaml` (the
nozzle-plume `diagnostics:` section; the cai2012 schema rejects new sections),
and calls the same `diagnose`. Every field, statistic and caveat is as
documented there:

* volume fields in the diagnosed time — `mfp` (VHS at the local T), `KnD`,
  `KnCellMax`/`KnCellV` and `cellOverMfpMax`/`cellOverMfpV` (resolution),
  `KnGradN`, `KnGradT`, `cellSize*`, `mfpResolved` (Δ_max ≤ λ/3 as a map),
  `parcelsPerStep`, `sampledParcels`, `diagValid*`, `diagStatus`;
* `diagnostics_summary.yaml` — validity, Δ/λ statistics, the worst cells, and
  lower bounds on a mesh that met λ/3. The mesh is graded Cartesian, not an
  octree, so no refinement levels are listed and the estimate is a lower bound
  on a lower bound (a graded mesh refines whole lattice planes);
* the centreline profile on **this family's** tube and stations
  (`post.centerline_*`), so a station is the same cells as in
  `centerline.csv`, drawn with the nozzle-plume function;
* images of every field from `viz-diagnostics.yaml` — R-4D's entries, with
  `KnCellMax` and `mfpResolved` on pinned ranges so the active Kn cases compare;
* one block of this family's, `exit_plane`: the measured `KnD` and Δ_max/λ in
  the exit cells beside the target Kn and the pre-run cell/λ0. It is the check
  that what generation planned is what the run produced. On a real
  `blockMesh` mesh with synthetic exit-state fields (§10) it measured
  `KnD = 10.12` against Kn = 10 and Δ/λ = 0.01236 against a planned 0.0125
  (the exit cells' `T` sits a little below `T0` because only the outgoing half
  of the Maxwellian is present).

**Validity at the baseline weight.** The diagnostics trust a cell only with
more than 20 parcels per step (`min_sampled_particles`, the nozzle-plume
default, `per_step` basis). With 40 parcels in the exit cell that holds near
the exit, where the density is above ~half of `n0`, and in coarse outer cells
whose volume makes up for the density; elsewhere the values are written and
drawn, flagged `diagValid = 0` and left out of every statistic. On the
synthetic check it was 1.2% of the cells and 10.8% of the volume. That is the
diagnostics reporting the statistical reach of one weight honestly, not a
defect; `docs/rarefaction-diagnostics.md` explains what lowering the threshold
or counting `accumulated` samples would cost.

---

## 8. Error metrics

`comparison.masked_metrics` reports, per quantity: valid and excluded counts;
max / mean / RMS **relative** error `|m − e|/|e|` where it is well
conditioned; max / mean / RMS **absolute** error `|m − e|/scale` everywhere it
is defined.

| quantity | relative where | absolute scale | excluded |
|---|---|---|---|
| `n/n0` (centreline, plane) | analytical `n/n0 ≥ 1e-3` (Cai's lowest contour) | 1 | below the floor; lip zone. An **empty cell inside the plume is scored** (100% deficit) |
| `n/n0`, `F` (r/D = 10) | analytical value ≥ 1e-3 × the sphere's peak | 1; peak `F` | below the floor; bins without cells. (At r/D = 10 the whole plume is below `n/n0 = 1e-3` past ~40°, so the plane's absolute floor would leave almost nothing.) |
| `u_x√β0`, `T/T0`, `p/p0` | as density, and DSMC molecules present | `S0`; 1; 1 | also cells with no sampled molecules — counted, not scored as zero |
| `u_r√β0`, `u_R√β0` | **never** — passes through zero on the axis | `S0` | as above |

Relative errors are fractions in the YAML and percents in the table.

**What these are not.** Not tolerances, not Cai's numbers, and at Kn ≤ 1 not
expected to be small: the collisionless solution is *supposed* to depart from a
collisional flow. The comparison measures that departure; at Kn = 10 it is
the near-collisionless check.

**Verification of the Cai figure/equation mapping.** The repository's
documentation cites Cai Eq. 5 and Figs. 6–8 for the density field and Eqs.
18–21 and Figs. 19–21 for the centreline, as pointers to where re-derived
results appear. **The paper itself is not in the repository**, and no
accessible primary copy was available here, so the correspondence of the
angular, off-axis velocity, temperature and pressure outputs to specific Cai
figures or equations is **not verified** and none is claimed. This blocks one
thing only — labelling those figures with Cai's numbers. It does not block
computing them: every formula used is a moment of the same free-molecular
distribution `analytical.moments` already integrates, and the definitions
(projections, `p = nkT`, the flux normalisation) are standard and stated
above. No digitised data are used anywhere.

---

## 9. Convergence: what counts as evidence

A flat total parcel count is **not** evidence of far-field convergence — it is
dominated by the near field and plateaus while the far centreline is still
filling (measured in `cases/cai2012`: flat at 1.32 transits, far centreline
still 5.6% low after 2). So:

* **`Allaudit.py`** (health, unchanged): measured parcels per cell per region,
  and the centreline error against Cai at every sampled frame — the running
  mean settling.
* **`Allconverge.py`** (new): the centreline tube in X/D bands (0–1, 2–3, 4–6,
  7–10) and the r/D = 9–10 shell in polar sectors (0–15°, 15–45°, 45–90°). Per
  region, (a) the running mean `rhoNMean` must have moved less than 2% over the
  second half of the window, and (b) the instantaneous `rhoN` must show no
  linear trend beyond max(2%, 2σ) across the window. A region is `stable`
  only if the trend test could resolve a 4% change; otherwise `inconclusive` —
  never a pass. The domain total is printed as context and excluded from the
  verdict.

Verdicts: `stable-window-demonstrated`, `not-demonstrated`, `inconclusive`,
`insufficient-frames`. The 2% tolerance is an analysis parameter
(`--tolerance`), chosen as the size of the far-field deficit the short
startup produced in the parent study.

---

## 10. Historical checks (2026-10-21, legacy four-case study)

The table below records the earlier four-case configuration and is not the
validation record for the new seven-case profile. Current preparation produced
70 dictionaries for seven cases, passed 84 independent verification checks, and
confirmed runtime-policy acceptance, idempotence and scaling registration;
see `results/full-suite-preparation/smoke.json`. No HPC execution is claimed.

| Check | How | Outcome |
|---|---|---|
| Exactly 4 configurations, the four Kn at the baseline | generator + `validate_cases.py` | pass |
| A second particle level is refused | policy test | pass |
| No excluded Kn generated, listed or run | policy tests; a hand-edited manifest with Kn = 0.01 → `AllmeshCases`/`AllpostCases` refuse; a case at Kn = 0.01 → `hpc guard` (`runCase.sh`, `runMesh.sh`) exits 2 | pass |
| D, S0, T0, gas identical in all 4 | manifest `held_fixed`; `validate_cases.py` | pass |
| `n0` reproduces Kn (independent VHS formula) | all 4 | < 1e-12 relative |
| Weight = `n0 V_exit / 40` | recomputed; `verify_tree` on `dsmcProperties` | pass |
| Startup = 10⁴ Kn = 0.01 collision times, same in all 4 | recomputed from `λ_ref = 0.01 D` | pass |
| Averaging 4.5 transits after startup | recomputed | 4.53–4.56 (snapped to a write) |
| Inheritance from the health template | key-by-key | 62 inherited, 5 declared overrides, 0 undeclared |
| Pre-run checks are the nozzle-plume checks, with its thresholds | unit test against `plumetools.nozzleplume.config.ChecksConfig` | pass |
| Idempotent generation | re-run: 4 × "keep (identical)"; a template change shows as `stale` until `--refresh` | pass |
| Conflicts reported, results never overwritten | unit tests: input edit → refused; with results → refused under `--regenerate`; only `--overwrite` replaces | pass |
| Dictionaries generate for all 4 | `./AllmeshCases --dict-only` (no OpenFOAM) | pass |
| `validate_cases.py` (inputs + dictionaries) | | 44 checks, 0 failures |
| Analytical limits | `u_r` even and zero on axis; `u_φ = 0`; banded quadrature = direct call to 1e-10; hemisphere flux = 1 to 5e-5 | pass |
| Estimators | a field equal to the analytical solution scores zero error through every estimator; a 10% deficit scores 10% | pass |
| Resolution diagnostics on a cai2012 case | unit tests: a uniform exit state gives `mfp = Kn D` and `KnD = Kn` exactly; the exit-plane check, the thresholds from `diagnostics.yaml`, this family's centreline tube, 40 vs 10 parcels against the validity rule, excluded Kn refused | pass |
| The whole per-case chain on a real mesh | **non-scientific** smoke case: `Kn10_np1x` coarsened to 58,080 cells, meshed by OpenFOAM's `blockMesh`/`topoSet`/`createPatch` (`Allmesh`, nozzle faces as predicted), synthetic analytical fields at OpenFOAM's own cell centres; then `Allpost`, `Alldiagnose.py`, `Allcompare.py`, `Allaudit.py`, `Allconverge.py` | all run; 18 diagnostic fields; particle-weight identity holds; exit `KnD` 10.12 vs Kn 10 |
| Diagnostic imagery | the smoke case through `viz-diagnostics.yaml` with ParaView 6.2 (offscreen) | 18 images |
| Unit tests | `tests/unit` (WSL, the repository `.venv`) | 1,812 passed, 6 skipped; 169 of them the extension's (`test_cai2012ext_*.py`) |

**Not checked here, because it needs the solver:** `dsmcFoam` accepting the
dictionaries at these densities (the format is the parent family's, which has
run), meshing the production-size meshes (the smoke case used the same
`Allmesh` on a coarsened copy), and the per-Kn flow-field imagery over a
sampled series (the specs are validated with the repository's own spec
loader). **No result in this document is a scientific validation;** that
requires completed HPC runs and the convergence evidence of §9.

---

## 11. Open problems

1. **Kn = 1e-4 and 1e-3 cannot be resolved in 3-D Cartesian at any budget this
   project has.** `dx = λ0` needs 1.8e13 and 1.9e10 cells. Even `dx = 10 λ0`
   at Kn = 1e-3 is ~2e7 cells in the core alone, and `Δt ≤ τ0` would then also
   need halving. Options: refinement confined to the
   near-exit region where `λ` is smallest (the plume rarefies within a few
   diameters); accepting under-resolution and labelling the run as such —
   what is done here: both are run, each with a `NON-PRODUCTION` marker
   written by its `runCase.sh`; or an axisymmetric solver, which is out of
   scope. A fourth, removing the dense core and injecting at a downstream
   sphere, is the interface track: [`cai2012-kn-interface.md`](cai2012-kn-interface.md).
2. **Write cadence through the startup** (§6). A controlDict with a coarse
   interval until `average_start_s` and the health cadence after it would cut
   the candidates' disk by ~4×; that is a change to the frozen
   `plumetools.cai2012.dictionaries`, so it is reported rather than made.
3. **Cost model** has no collision work. The first completed candidate will
   replace it with a measurement.
4. **Cai figure numbers** for outputs 3–6 (§8).
5. **One parcel weight.** With no sweep, statistical adequacy rests on each
   run's own evidence — the occupancy audit, `diagValid`, the stationarity
   verdict — and at 40 parcels per exit cell the diagnostics' own validity rule
   (more than 20 per step) holds near the exit only (§7a). Whether that is
   enough is a judgement on the first completed run, not one the configuration
   can make.

### Observed in the frozen packages (reported, not fixed)

* `cases/cai2012/generate_cases.py` with `--knudsen` / `--particles` writes a
  manifest containing **only the selected cases**, so a filtered regeneration
  silently drops the rest from every runner's schedule.
* `plumetools.cai2012.study.clone_cases` deletes and re-copies any existing
  case without a top-level time directory — including a meshed one, or one whose
  results are only under `processor*/`. This is why the extension has its own
  idempotent generation.
* The `mapfile -t CASES < <(python ...) || die` pattern in the parent runner
  scripts cannot see a failure inside the process substitution.
* `docs/cai2012-case.md` §10 item 6 says off-axis analytical velocity and
  temperature are not implemented; `analytical.moments` implements them.

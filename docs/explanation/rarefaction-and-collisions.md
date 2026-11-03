# Rarefaction, collisions and the criteria built on them

The families define a Knudsen number, a mean free path and their resolution
criteria in slightly different ways, for reasons each records. This page puts
the definitions side by side with their code locations, so a number from one
family is not read with another's convention. It summarises; each family's own
page has the reasoning.

## Mean free path

All families start from the hard-sphere expression at the VHS reference
diameter d and differ in the temperature correction:

$$
\lambda_{\mathrm{HS}} = \frac{1}{\sqrt2\,\pi d^2 n},
\qquad
\lambda_{\mathrm{VHS}} = \lambda_{\mathrm{HS}}\left(\frac{T}{T_{\mathrm{ref}}}\right)^{\omega - 1/2}.
$$

| Where | Conventions | Default | Notes |
|---|---|---|---|
| `cai2012.gas.mean_free_path` (G1) | `vhs`, `hard_sphere` | `vhs` | `gas.mean_free_path_convention`; the factor is 1.0229 for the committed argon at 300 K. Cai does not state GRASP's convention: an assumption |
| `nozzleplume.gas.mean_free_path` | `vhs`, `hard_sphere` | `vhs` | r-4d uses ω = ½, where the two coincide |
| `rarefaction.core.mean_free_path` | `vhs`, `hard_sphere` | `vhs` | element by element; NaN (not ∞) for n ≤ 0, because there it means *no data* |
| `hybridplume.physics.Gas.mean_free_path` | `vhs`, `hard_sphere`, `viscosity_simple` | `vhs` | `viscosity_simple`: λ = (μ/p)√(πRT/2), from the gas's VHS viscosity law |
| `cai2012ext.interface.mean_free_path` | as passed | — | the closure state on candidate spheres |
| `markelov1999.resolution.mean_free_path` | hard sphere only | — | deliberately uncorrected: it feeds a warning threshold, and T varies by orders of magnitude across the domain (its docstring writes the omitted factor inverted; see [known gaps](../development/known-gaps.md)) |

The mean collision time used by the time-step criteria is τ = λ / c̄ with
c̄ = √(8RT/π) (`hybridplume.physics.Gas.collision_time`; the Cai and nozzle
families' `collision_time` functions follow the same form).

## Knudsen numbers

| Quantity | Definition | Where |
|---|---|---|
| exit Kn (Cai) | λ₀ / D, `exit.characteristic_length: diameter` (Cai's statement); `radius` halves λ₀ and doubles n₀ | `cai2012.gas.mean_free_path_from_kn`, `number_density_from_kn` |
| exit Kn (nozzle plume) | λ_exit / D, derived from the given exit state, never an input unless `exit.knudsen` is chosen | `nozzleplume.exitstate.ExitState.knudsen` |
| reference state of Cai's criteria | Δx/λ₀ and Δt/t₀ are referred to the **Kn = 0.01** exit properties for every case (`resolution.reference_knudsen`) | `cai2012.inflow` |
| local `KnD` | λ / D per cell | `rarefaction.core` |
| local `KnCellMax`, `KnCellV` | λ / Δ_max, λ / V^(1/3) — mesh resolution | `rarefaction.core` |
| local `KnGradN`, `KnGradT` | λ \|∇n\|/n, λ \|∇T\|/T — rarefaction on the flow's own scales | `rarefaction.core` |
| hybrid `Kn_rho`, `Kn_T`, `Kn_U`, `Kn_GLL` | λ\|∇ρ\|/ρ, λ\|∇T\|/T, λ\|∇\|U\|\|/max(\|U\|, c); Kn_GLL is their maximum (Boyd, Chen & Candler 1995) | `hybridplume.breakdown` |

The cell and gradient numbers measure different things: the first measures the
mesh, the second the flow. None of them, and no threshold on one, says a DSMC
solution is valid — the time step and the parcels per cell matter as much
([Rarefaction diagnostics](../rarefaction-diagnostics.md)).

## Resolution and time-step criteria

| Criterion | Value | Where it is enforced |
|---|---|---|
| cell / λ₀ ≤ 1 (Cai's Δx/λ₀ = 1) | `checks.max_cell_over_mfp` 1.0 — a warning | Cai and nozzle-plume checks |
| λ / Δ_max ≥ 3 for a resolved cell | `diagnostics.cells_per_mean_free_path` 3.0 | `rarefaction` (`mfpResolved`) |
| Courant on U₀ + 3σ through the smallest cell | `dsmc.courant_target` 0.2 sets Δt; `checks.max_courant` 1.0 is a **hard error** | Cai, nozzle-plume derivation and checks |
| Δt / τ_exit ≤ 1 | `checks.max_time_step_over_collision_time` 1.0 | nozzle plume; `cai2012ext.feasibility` |
| ≥ 5 parcels per cell (floor), 20 targeted in one region | `checks.min_particles_per_cell`, `resolution.target_particles_per_cell` | Cai, Markelov |
| > 20 sampled parcels for a valid diagnostic cell | `diagnostics.min_sampled_particles` | `rarefaction` |

Standard `dsmcFoam` has **one** particle weight for the whole domain
(`DSMCCloud::nParticle_`), so an occupancy target is met in exactly one region
and every other region follows from the density and cell-volume ratio. The
families report the occupancy that implies everywhere, and the audits measure
it after the run (`cai2012.audit`, `markelov1999.resolution.audit`).

## Collisionless and collisional

- **Cai's analytical solution is collisionless** (free molecular): the
  distribution at a point is the exit Maxwellian restricted to the directions
  whose backward ray meets the disk (`cai2012.analytical`).
- **DSMC runs with collisions on at every Kn**, including Kn = 100: switching
  them off would make agreement true by definition rather than demonstrating
  that DSMC approaches the collisionless limit (`dsmc.collisions_enabled`).
  The hybrid's `smoke-collisionless` and `verification-collisionless` profiles
  do switch them off, on purpose, to test the DSMC machinery against an
  independent collisionless reference (`hybridplume.collisionless`).
- **The virtual-collision model** adds reduced, mutually thinned
  representative pair events to the collisionless field, with OpenFOAM's VHS
  kernel; it is a correction model under assessment, not a DSMC substitute.
- **The scaling maps** deform the collisionless field's coordinates only; they
  conserve nothing and carry no collision physics.

## Breakdown and interface criteria

Two families decide where continuum or source models stop being trusted.

| Track | Parameter | Definition | Threshold |
|---|---|---|---|
| Kn-extension interface | Bird's B, `general` | (2λU/(r c̄)) · M²/(M²−1) — the exact −(λ/(ρ c̄)) U·∇ρ on a radial streamtube; NaN at M ≤ 1 | `bird_threshold` 0.05 (validation set 0.02, 0.03, 0.05) |
| | Bird's B, `radial_approx` | Stewart's Eq. (3), without the acceleration factor | |
| Hybrid | `B_exp` | −τ_c U·∇ln ρ, signed (positive in an expansion), with a nominal τ_c = λ/c̄ | `breakdown.thresholds.B_exp` 0.05 |
| | `Kn_GLL` | as above | `breakdown.thresholds.Kn_GLL` 0.05 |

The hybrid then judges the chosen surface on what crosses it: returning mass
fraction ≤ 0.05, subsonic-normal area fraction ≤ 0.25, no contact with the
outer boundary, at least one overlap layer (`interface.criteria`). A one-way
Maxwellian inlet cannot carry viscous stress or heat flux, and absorbs the
returning flux; where that flux is large the hybrid is unsuitable and is
labelled so.

## Conservation checks

| Check | What is conserved | Where |
|---|---|---|
| inflow flux verification | mass, momentum, energy through the meshed inflow, two independent routes | `markelov1999.flux` |
| interface flux balance | the closure's ∫ j dΩ against the exit mass flow; meshed against ideal sphere flux (`meshed_flux_tolerance` 2 %) | `cai2012ext.interface` |
| plate loads | force = jet momentum flux + diffuse recoil; heat = jet energy flux − re-emitted energy, for a collisionless jet onto a plate that catches all of it | `tests/openfoam/test_nozzleplume_openfoam.py` |
| hybrid budgets | control-volume mass, momentum and energy of the precursor, the DSMC region and the interface | `hybridplume.budgets` |
| injected flow | parcels injected per step against the Maxwellian flux prediction | Cai and nozzle-plume OpenFOAM tests; hybrid `dsmc-post` |

The staircase nozzle's area error is a direct multiplier on every injected and
continuum flux: it is measured, reported and checked against a limit, never
rescaled.

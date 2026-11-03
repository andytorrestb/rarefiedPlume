# The Cai 2012 extension's interface track

Removing the dense core from the DSMC domain at `Kn_e = 1e-3, 1e-4`, and
injecting a continuum continuation of Cai's exit on a sphere downstream.

Running instructions are in the tutorials: [direct-exit cases](tutorials/cai2012-kn-extension.md)
and [this interface track](tutorials/cai2012-interfaces.md). The extension itself is
[`docs/cai2012-kn-extension.md`](cai2012-kn-extension.md). This document is the
method and the findings of the interface track.

> **Status (2026-10-22): implemented and smoke-tested; nothing validated.**
>
> * The overlap validation **cannot run as declared**. At `Kn_e = 0.01` Bird's
>   B on the plume axis is already 0.062 at the equivalent anchor sphere and
>   grows outward, so B = 0.02, 0.03 and 0.05 admit no sphere at any core
>   angle (§5). All three validation cases are generated and recorded as
>   `NO_ADMISSIBLE_INTERFACE`.
> * The two production cases are blocked on two decisions this document does
>   not make for you: the core angle and the source-validity radius (§4).
> * Markelov's angular law is a sonic-orifice law. Imposed on Cai's
>   supersonic exit it puts 0.36–0.45 times Cai's collisionless density on
>   the axis (§6). This is measured and reported, not tuned.
> * One tiny development case ran end to end in OpenFOAM v2512 on 1 and 2
>   ranks (§11). That checks the mechanics only.

---

## 1. What this increment is

| | |
|---|---|
| Physical exit | Cai's: `D_e = 0.2 m`, `S0 = 2`, `T_e = 300 K`, argon VHS, `Kn_e = λ_e/D_e` with the repository's VHS convention. None of these change. |
| Injection | A sphere centred on the exit, `R_interface` chosen by Bird's B; patch `nozzle` (same name, so `plumeFieldInflow` and `dsmcProperties` are unchanged), role: **continuum interface**. |
| State injected | Markelov's angular law with a finite-temperature radial-streamtube closure (§3), per face, in patch order. |
| Where | Only `Kn_e = 1e-3, 1e-4` (production), `Kn_e = 0.01` (validation, against the existing full-domain benchmark), and development smoke cases. `Kn_e = 10, 1` and Cai's three cases stay full-domain. |
| What is compared | Every case still against **Cai's analytical solution launched from the physical exit** (`Allcompare.py`). Validation cases also against the **full-domain DSMC run** (`Alloverlap.py`). |

Not in this increment: an angle-dependent surface, CFD coupling, AMR, plates,
forces.

---

## 2. Audit (Stage 1)

**Branch.** `cai2012-kn-extension-pbs` at `de071b6`, the inspection baseline:
no commits since. Uncommitted before this work: `cases/r-4d/baseCaseHPC/case.yaml`
(and its two logs), the build products of `plumeFieldInflow`, and the prompt
file. None of them was touched.

**Solver.** Standard OpenFOAM v2512 `dsmcFoam` with
`libplumeDsmcBoundaryModels.so` (`plumeFieldInflow`), confirmed in the
generated `controlDict`/`dsmcProperties` and by running it (§11). Not MNF
`dsmcFoam+`.

**The Kn = 0.01 benchmark** (`cases/cai2012/Cases/Kn0p01_np1x`) is **not in
this checkout**: `cases/cai2012` holds its template and study only. It may be
on the cluster. Nothing here reads it except `Alloverlap.py`, read-only.

| Need | Status | Where |
|---|---|---|
| Markelov angular law, `θ_L`, `I`, `A_M`, `V_inf` | **reused, frozen** | `plumetools.markelov1999.sourceflow` |
| Half-range drifting-Maxwellian face fluxes | **reused** as the cross-check (tests) | `plumetools.markelov1999.flux` |
| Kn → n_e, exit state, VHS mean free path | **reused** | `plumetools.cai2012.{inflow,gas}` |
| Graded box lattice, blockMeshDict | **reused** (re-planned around the sphere) | `plumetools.cai2012.mesh` |
| Sphere carving, snappyHexMeshDict | **reused** | `plumetools.foamio.primitives` |
| Run settings (Δt, weight, startup, writes) | **reused**, fed the interface state | `plumetools.cai2012.inflow.derive_run_settings` |
| Field writers, measurement fields | **reused** | `plumetools.cai2012.foamfields`, `markelov1999.foamfields` |
| Staging, idempotent generation, conflicts | **reused** | `plumetools.cai2012ext.generation` |
| Cai analytical solution | **reused** | `plumetools.cai2012.analytical` |
| Per-face inflow | **no C++ change** | `plumeFieldInflow` already reads per-face n, U, T |
| Equivalent-source mapping, closure, B, sphere selection | **missing → added** | `cai2012ext/interface.py` |
| Truncated mesh composition, patch verification | **missing → added** | `cai2012ext/interface_mesh.py` |
| One derivation path, assessment, identity, inflow | **missing → added** | `cai2012ext/interface_case.py` |
| Interface study, tree, manifest | **missing → added** | `cai2012ext/interface_study.py`, `interfaces.yaml` |
| Overlap metric with uncertainty | **missing → added** | `cai2012ext/overlap.py` |

`plumetools.sourceflow` (the legacy model) is not used. `cases/2d-wedge/calcCaseParams.py`
was not used either; its spherical formulas are superseded by §3.

---

## 3. The model

### 3.1 What is published, and what is assumed

`[PAPER]` Lumpkin, Stewart & Markelov, AIAA 99-3455, p. 3: `θ_L`, `f(θ)`, `I`,
`A_M`, `V_inf`, `ρ_M`. Stewart & Lumpkin (2012), Eqs. (1)–(4): Bird's B, the
spherical area-Mach relations, the simplified radial B and its hard-sphere
inversion. Cai & Wang (2012): the exit.

`[ASSUMPTION]` — this project's coupling, not a published result:

| | |
|---|---|
| CP-1 | Source at the exit centre, `x_s = x_e`, on the exit axis. Stored separately from the exit origin. |
| CP-2 | Conserved targets: the exit's **gross** outward molecular mass and energy fluxes. |
| CP-3 | Mass per steradian has Markelov's shape at every radius, normalised by the full plume's `I`. |
| CP-4 | Independent isentropic radial streamtubes: `ρ U r² = j(θ)`, `c_p T + U²/2 = h_target`. No transverse pressure gradient. |
| CP-5 | Anchor: `U(r_a) = U_e`, `ρ(r_a, 0) = ρ_e`. An equivalent sphere, not a nozzle exit. |

### 3.2 The equations

Exit (`a = √(2RT)`, `s = u_n/a`, `H = erfc(−s)/2`, `E = e^{−s²}/(2√π)`):

```
J1 = u_n H + a E           J3 = (u_n³ + 1.5 u_n a²) H + a (u_n² + a²) E
mdot = A_e m n_e J1        Edot = A_e (m n_e/2) [J3 + 2 R T_e J1]       h = Edot/mdot
```

For `s < 0` both moments are evaluated with `erfcx`, so there is no
cancellation down to exact underflow (`s ≈ −27`).

Source and closure (`γ = 5/3`, `θ_L = π/2`, `k = (γ+0.41)/(γ−1) = 3.115`):

```
j(θ) = mdot f(θ) / (2π I)          T_t = h/c_p      T_a = T_t − U_e²/(2 c_p)
M_a = U_e/√(γ R T_a)               r_a = √(mdot / (2π I ρ_e U_e))
F(M) = (1/M)[(2/(γ+1))(1+(γ−1)M²/2)]^((γ+1)/(2(γ−1)))     F(M(r)) = F(M_a)(r/r_a)²
T = T_t/(1+(γ−1)M²/2)   U = M√(γRT)   ρ = j(θ)/(U r²)   U_vec = U r̂
```

The Mach number is solved on the supersonic branch only (vectorised bisection
in `ln M`).

Bird (signed, positive in an expansion; `λ` is the case's own VHS convention
at the local T):

```
B_general       = (2 λ U / (r c̄)) M²/(M² − 1)      default; exact on a radial streamtube
B_radial_approx = 2 λ U / (r c̄)                    Stewart Eq. (3)
```

### 3.3 The numbers (all three Kn share every shape quantity)

| | Kn = 0.01 | 1e-3 | 1e-4 |
|---|---|---|---|
| `n_e` [m⁻³] | 6.620e20 | 6.620e21 | 6.620e22 |
| `mdot` gross [kg/s] | 9.750e-4 | 9.750e-3 | 9.750e-2 |
| gross / net mass | 1.000245 | ← | ← |
| `h` [J/kg] | 4.0599e5 (net 4.0607e5) | ← | ← |
| `T_t`, `T_a` [K] | 779.85, 299.85 | ← | ← |
| `M_a`, `r_a/D`, `r*/D` | 2.1915, 0.7173, 0.5444 | ← | ← |
| `V_inf` [m/s] | 901.10 | ← | ← |
| `p_t r_o²` [Pa m²] | 0.4306 | 4.306 | 43.06 |
| B on the axis at `r_a`, general / radial | **0.0624** / 0.0494 | 0.00624 / 0.00494 | 0.000624 / 0.000494 |

B is exactly proportional to `Kn_e`: everything else in the mapping is
independent of the density.

### 3.4 The 1999 terminal model, compared

Mapped through the same `p_t r_o²`, Markelov's printed density **is** the
closure's `U → V_inf` limit: `ρ_M = j/(V_inf r²)`, checked to 1e-12 against
the frozen `markelov1999.sourceflow.mass_density`. The closure is denser by
`V_inf/U(r)`: 1.139 at 1 D, 1.046 at 2 D, 1.013 at 5 D, 1.005 at 10 D. The
terminal model's own local pair `(V_inf, T_t)` holds `c_p T_t + V_inf²/2 = 2h`:
twice the available energy, which is why it is compared and never injected.
Every `interface-report.yaml` carries this comparison
(`terminal_model_comparison`).

---

## 4. Decisions this document does not make

**Core angle `θ_c`.** No `Kn_e = 0.01` DSMC profile is available locally, so no
angle is proposed for a scientific run: `core_half_angle_deg` is `null` in
every production and validation entry, and generation reports
`CONFIGURATION_INCOMPLETE`. The smoke case uses 30°, marked provisional.

What the choice costs, before any DSMC result (`generate_interfaces.py
--report`, `decision_table`, source-validity radius at its floor `r_a`):

| θ_c | injected mass | Kn 1e-3, R/D at B = 0.02 / 0.03 / 0.05 | Kn 1e-4, R/D at B = 0.02 / 0.03 / 0.05 |
|---|---|---|---|
| 10° | 6.1 % | 1.43 / 1.88 / 2.70 | 7.37 / 9.5* / 9.5* |
| 20° | 22.6 % | 1.30 / 1.70 / 2.43 | 6.62 / 8.91 / 9.5* |
| 30° | 44.7 % | 1.10 / 1.43 / 2.03 | 5.49 / 7.39 / 9.5* |
| 45° | 76.0 % | 0.75 / 0.96 / 1.32 | 3.47 / 4.66 / 6.77 |
| 60° | 94.2 % | none / none / none | 1.63 / 2.15 / 3.09 |

\* limited by the 10 D box (and Cai's r/D = 10 sphere), not by B. Enlarging
the domain and the diagnostic radius is an explicit `geometry.*` /
`diagnostic_radius_over_D` change, recorded as such.

The tail beyond `θ_c` is **not injected** and **not redistributed** into the
core; the injected fraction is `∫₀^θc f sin θ dθ / I` with the full plume's `I`.

**Source-validity radius.** No paper certifies a minimum radius for this
mapping. `source_validity_min_radius_over_D` is `null` for every scientific
case; the smoke case uses 1 D, marked provisional.

**The thresholds.** 0.02, 0.03 and 0.05 are this study's sensitivity choices.
§5 shows why they cannot be used at `Kn_e = 0.01`.

---

## 5. The Kn = 0.01 overlap validation is infeasible as declared

The selector searches outward from `r_min = max(r_a, validity radius)`. On the
axis — which every core contains — `B_general` is 0.0624 at `r_a` and only
grows outward (0.110 at 1 D, 0.312 at 2 D, 2.89 at 10 D). Hence:

* B = 0.02, 0.03, 0.05: **`NO_ADMISSIBLE_INTERFACE` for every core angle and
  every validity radius.** The exit state of Cai's own `Kn = 0.01` case is
  already past these thresholds by Bird's criterion.
* Under Stewart's simplified form the axis value at `r_a` is 0.0494, just
  under 0.05 — but any core of more than a few degrees exceeds 0.05 at its
  edge, so the conclusion stands under both definitions.

Nothing was changed to make a sphere fit. The three validation cases are
generated, refuse to mesh, and say why in `results/interface-report.yaml` and
`interface-manifest.yaml` (`validation.B_at_smallest_admissible_radius`). The
production cases therefore cannot acquire `coupling_validation: passed`, and
`runCase.sh` refuses them unless `run_mode: exploratory`.

Ways forward, each a study decision for you, none taken here:

1. Validate at thresholds the Kn = 0.01 plume admits (the smallest is ≈ 0.062
   on the axis, ≈ 0.098 at a 30° core edge, at `r_a`) — a different, looser
   continuum criterion, and to be labelled as such.
2. Validate at a lower Kn, which needs a resolved full-domain benchmark there
   — at `Kn_e = 1e-3` the direct exit needs ~2e10 cells.
3. Accept exploratory production runs, labelled NON-PRODUCTION, and treat the
   coupling as unvalidated.

---

## 6. The angular law does not fit this exit

Markelov's `f(θ)` is the law of a **sonic** orifice: `θ_L = 90°` at γ = 5/3.
Cai's exit is at `M_a = 2.19`, whose Prandtl–Meyer turning limit is 64.3°.
Imposing the sonic law (CP-3) spreads the mass much more than Cai's
collisionless solution does (`angular_shape_comparison`, at r/D = 10):

| inside θ | 5° | 10° | 15° | 20° | 30° | 45° | 60° |
|---|---|---|---|---|---|---|---|
| Markelov | 1.6 % | 6.1 % | 13.3 % | 22.6 % | 44.7 % | 76.0 % | 94.2 % |
| Cai collisionless | 4.1 % | 15.2 % | 30.8 % | 47.8 % | 75.9 % | 95.1 % | 99.4 % |

On the axis the closure's density is **0.45×** Cai's collisionless value at
2 D, **0.38×** at 5 D and **0.36×** at 10 D. Cai reports his `Kn = 0.01` DSMC
within ~12 % of the collisionless centreline, so an overlap test of this
closure at the 5 % target would be expected to fail by a wide margin. The
mismatch is stated in every assessment (`angular_law_axis_density`, not
blocking: the collisionless solution is not the truth at low Kn — the overlap
study is) and `f` is not tuned. The smoke run's centreline sits 55 % below
the collisionless value at 2–3 D, consistent with this.

The mapping itself is conservative: on the ideal sphere the core injects its
share of the exit's gross mass and energy to 7e-6 (the quadrature of the
reported sphere), and gross and net differ by 0.02 % at the exit but only
~1e-15 at the interface, where `s ≈ 5`.

---

## 7. Mesh, inflow and numerical settings

**Mesh.** `plumetools.cai2012.mesh.plan` with the uniform core re-sized to
enclose the sphere plus `min(0.25 D, ½ domain_clearance)`, the core cell set to
`target_cell_over_mfp × λ(R, 0)` (capped as the parent caps it), and the outer
cell counts limited so no outer cell is smaller than the core cell. The budget
counts **retained** cells (lattice centres outside the sphere). Then
snappyHexMesh at level 0 removes `r < R` and snaps to the analytic sphere. No
topoSet/createPatch, no wall, no symmetry plane; `x = 0` outside the sphere
stays the open `upstreamVacuum`.

**Patch verification** (`interface_case verify-mesh`, from the patch's own
geometry): vertices within 2 % of R, least-squares centre at the source, area
within `checks.max_nozzle_area_error` of 2πR², every face oriented with its
OpenFOAM area vector toward the source, enough faces in the core, one region
(checkMesh). OpenFOAM's right-hand-rule orientation is computed explicitly:
`plumetools.geometry.normal` uses the opposite, legacy orientation (§12).

**Inflow** (`runInflow.py`, per face in patch order): `n` in the core, exactly
0 beyond it; `U(r) r̂` and positive `T(r)` on every face (`plumeFieldInflow`
aborts on a zero `boundaryT`); fields read back and compared in patch order;
gross and net mass and energy integrated over the actual faces against the
core targets. A scientific case whose meshed gross flux misses by more than
`meshed_flux_tolerance` (2 %) is refused — a failed coupling, not something to
rescale.

**Δt and weight.** The declared reference cell is the core cell on the +x axis
half a cell beyond R: the densest retained cell. `nEquivalentParticles =
n_ref V_ref / 40`. `Δt = min(Courant 0.2 × smallest cell / (|U| + 3√(RT))
at the interface, the reference cell's collision time)`. The startup is still
Cai's 10⁴ reference collision times and the averaging window 4.5 transits of
`x_max/U_e`, unchanged from the direct cases.

**Identity.** `Allmesh` writes `interface-identity.yaml` (hashes of `case.yaml`
and the `interface:` section, R, core cell, divisions, weight). `runInflow.py`
refuses a mesh that is not this configuration's; the PBS guard refuses a case
with solver output whose identity is missing or different. Old direct-exit
parcels can never be resumed on an interface mesh.

What the production cases would be at an illustrative 30° core and 1 D
validity radius (not a proposal):

| | Kn 1e-3 | Kn 1e-4 |
|---|---|---|
| R/D, limited by | 2.033, B | 9.5, domain |
| retained cells (budget 4e6) | 3.92e6, coarsened | 3.90e6, coarsened |
| core cell / λ(R,0) | 5.0 | 8.8 |
| cells for cell = λ | 1.5e8 | 2.7e9 |
| parcels, steps | 2.4e7, 55 892 | 5.2e6, 19 402 |

Both stay **unsuitable for production** inside the inherited 4e6-cell budget.
The interface cuts the criterion's cell count by ~130× and ~7000× against
the direct exit (1.9e10, 1.8e13) — not enough on its own.

**Memory.** Bytes per cell and per parcel are not calibrated for `dsmcFoam`
at these densities, so the 1.5 TB ceiling (`resources.memory_ceiling_bytes`,
in bytes) is reported as `UNCALIBRATED` until you set `bytes_per_cell` /
`bytes_per_parcel` from a measured run. The repository's PBS scripts request
no memory, so no site syntax for it is assumed (`memory_ceiling_pbs: null`).
Walltime stays the existing queue defaults.

---

## 8. The overlap comparison

`plumetools.cai2012ext.overlap`, run by `Alloverlap.py` in a validation case.

* **Reference:** the full-domain DSMC run at the same Kn, read-only. Not
  Cai's analytical solution.
* **Region Ω, declared before any result** (`overlap.region`): exit-centred
  sector `4 ≤ R/D ≤ 9.5`, `θ ≤ 15°`, bins of 0.5 D × 2.5°, all azimuths. Refused
  if it starts within `interface_buffer_over_D` (1 D) of any interface, leaves
  a domain, or exceeds the injected core angle.
* **Estimator** (mesh-independent): cells assigned to bins by centre; `n`
  volume-weighted, `|U|`, `u_R`, `u_θ`, `T` molecule-weighted per bin.
  Empty bins in either run are excluded and counted — never scored as zero.
* **Metric:** `E_q = √(Σ w (q_if − q_full)² / Σ w q_full²)`, `w` the exact
  sector volume, for `n`, `|U|`, `T`; targets 0.05, 0.05, 0.10 (project
  targets, configurable). `u_R` and `u_θ` are reported beside them.
* **Uncertainty:** block means recovered exactly from fieldAverage's running
  means, `M_ab = ((t_b−t_s)M_b − (t_a−t_s)M_a)/(t_b−t_a)`, six blocks; per-bin
  standard errors; a linearised, noise-corrected interval for `E`. **Pass** only
  if the upper bound is below the target; **fail** only if the lower bound is
  above it; otherwise — or with fewer than three blocks, or a noise floor of
  half the target — **inconclusive**.
* Only the block-boundary frames are read, one at a time, reduced to Ω's
  cells.

Output: `results/overlap.yaml` (verdicts, intervals, excluded and low-signal
counts, the region and definitions) and `results/overlap_bins.csv`.

---

## 9. Gates and labels

| status | meaning | runCase.sh |
|---|---|---|
| `interface-infeasible` | no admissible sphere, or a decision missing | guard refuses (exit 2) |
| `unsuitable-for-production` | a production case misses a numerical criterion | label refuses (exit 3) unless `run_mode: exploratory` |
| `pending-coupling-validation` | numerically sound, no passed overlap at this threshold | label refuses unless exploratory |
| `production-candidate` | everything met and `validation_evidence` passed for the same coupling (below) | runs, unmarked |
| `interface-validation` | an overlap run | runs if numerically sound (or exploratory), marked NON-PRODUCTION |
| `development-smoke` | mechanics only | runs, marked NON-PRODUCTION |

`numerical.status` and `coupling_validation.status` are reported separately.
The direct-exit cases' gates are unchanged.

**What evidence counts** (`plumetools/cai2012ext/interface_evidence.py`, since
2026-10-27). `Alloverlap.py` records in `overlap.yaml` the coupling it tested,
and a production case may rest on a `pass` only if every one of these is
recorded and equal to its own: the model, `purpose: interface-validation`, the
Bird definition and threshold, the core half-angle, the tail treatment, the
angular exponent offset, and the source problem (D, S0, T0, the gas and its
collision model). The Knudsen number, R, the mesh, the time step and the weight
are recorded but never compared: validating at Kn = 0.01 and carrying the
verdict to 1e-3 and 1e-4 is the design. A pass for another coupling reports
`incompatible` (with the items), and evidence that does not say which coupling
it tested -- anything written before 2026-10-27 -- reports
`insufficient-provenance`; re-running `Alloverlap.py` on the validation case
rewrites it from the same fields. Before 2026-10-27 a `pass` with a matching
threshold was enough. The core angle is the one conservative choice here: the
overlap region is confined to the injected core and the excluded tail shapes
what reaches it, so a pass at one angle is not taken to validate another.

---

## 10. Commands

The laptop checks, the development smoke run and the cluster sequence moved to
the [interface-track tutorial](tutorials/cai2012-interfaces.md), which is now the
one place they are maintained. §11 below records what was actually run on
2026-10-22.

---

## 11. Checks performed (2026-10-22, laptop)

| kind | what | outcome |
|---|---|---|
| source inspection | §2 | as tabled |
| Python, physics (`test_cai2012ext_interface.py`, 66) | θ_L = π/2; ∫ j dΩ = mdot; core fraction = 1 − cos^(k+1) θ_c; anchor = exit state; ρUr² = j and h conserved at 20 (r, θ); area-Mach and isentropic density; supersonic branch; U → V_inf; terminal density = j/(V_inf r²) through the frozen module; Stewart Eq. (2) and Eq. (4) inversion; finite-difference d ln ρ/dr = B_general; ratio M²/(M²−1); J1, J3 vs quadrature for s ∈ [−8, 12], limits, stability to s = −1e4; face fluxes = `markelov1999.flux`; Kn scaling; selector: root, edge breakdown, no root at Kn 0.01, insufficient domain, domain-limited, invalid, incomplete; truncation without renormalisation; config validation | pass |
| Python, integration (`test_cai2012ext_interface_case.py`, 26) | generation beside `Cases/`, idempotent, conflicts/stale; track policy; retained-cell count vs brute force; snappy dictionary; infeasible case writes no mesh; per-face fields in patch order on a synthetic hemisphere; wrong winding/radius/identity refused; guard, label, evidence | pass (1 skipped on native Windows, runs under WSL) |
| Python, overlap (`test_cai2012ext_overlap.py`, 8) | config, region, exact sector volumes, metric, interval verdicts, empty-bin exclusion, exact block recovery | pass |
| Python, existing extension tests | generation, policy, diagnostics, comparison, convergence, hpc | see the handoff |
| solver smoke, 1 rank | `smoke_Kn0p001_B0p05`: 36,946 cells, one region, area −0.10 %, R snapped to 1e-9; meshed gross flux −0.96 % of the core target; 213 steps in 26 s; injection 4.27e-3 kg/s from the first 4 steps vs 4.31e-3 predicted; 280 parcels/step vs 279; steady 18,802 parcels vs 1.90e4 estimated | ran |
| solver smoke, 2 ranks | same, fresh: decomposed inflow = reconstructed exactly (1513 faces); 18,792 parcels; post chain (Allpost, Alldiagnose, Allcompare, Allconverge, collect) ran; overlap mechanics on serial vs 2-rank: E_n 0.019 against a noise floor of 0.016 | ran |
| **not done** | any production-size mesh; any overlap validation; any scientific run; HPC; imagery for interface cases (`runViz.sh` refuses) | — |

---

## 12. Limitations and observations

* **Unmet:** overlap validation (§5); scientific core angle and validity
  radius (§4); the angular law (§6); mesh budget for both production cases
  (§7); memory calibration.
* The closure ignores transverse pressure gradients and finite-exit structure
  (CP-4); the sphere's face centroids lie up to 0.15 % inside R on the smoke
  mesh; the per-face core cut is by face centroid.
* The bin temperature uses cell values, so it carries no intra-bin velocity
  spread.
* `Allaudit.py` (the direct-exit occupancy regions) is not run for interface
  cases; Alldiagnose's exit-plane check is skipped there.
* **In shared code, reported, not changed:** `plumetools.geometry.normal`
  returns `(p0−p2)×(p0−p1)`, the opposite of OpenFOAM's face area vector; the
  interface code computes OpenFOAM's orientation itself.
  `markelov1999.sourceflow.normalization_integral` is accurate to ~1e-11 (scipy
  `quad` defaults), not 1e-14 as its docstring says.
* **Local environment:** the WSL user bin's `decomposePar` links hyStrath's
  `libdecompose+.so` and fails on every case, OpenFOAM's cavity tutorial
  included; `$FOAM_APPBIN` first on `PATH` fixes it. The cluster environment
  was not checked for the same shadowing.

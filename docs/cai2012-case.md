# The Cai & Wang 2012 case family

Cai, C., and Wang, L., *Numerical Validations for a Set of Collisionless Rocket
Plume Solutions*, Journal of Spacecraft and Rockets, Vol. 49, No. 1, 2012.

This document separates **what the paper states** from **what this
implementation had to decide**, and gives the reasoning for each decision. The
running instructions are in
[`cases/cai2012/README.md`](../cases/cai2012/README.md).

> **This is not a reproduction of Cai's DSMC results, and must not be described
> as one until the numerical comparison supports it.** Cai used an axisymmetric
> DSMC code (GRASP) on a uniform grid. This is three-dimensional Cartesian
> OpenFOAM `dsmcFoam` v2512, at a different resolution, with a different
> collision-partner selection implementation and different molecular constants.
> What is being validated is a **trend** — see §8.

---

## 1. Values stated by the paper

| Quantity | Paper | In this implementation |
|---|---|---|
| Nozzle diameter `D` | 0.2 m | `nozzle.diameter_m` |
| Nozzle radius `R0` | — | **derived** as `D/2` = 0.1 m; there is no radius key |
| Exit speed ratio `S0` | 2.0 | `exit.speed_ratio` |
| Knudsen numbers | 100, 0.1, 0.01 | `study.yaml`, one case each |
| Characteristic length for `Kn` | the nozzle **diameter** | `exit.characteristic_length: diameter` |
| Gas | argon | `gas.species_name: Ar` |
| Collision model | Variable Hard Sphere | `dsmc.binary_collision_model: VariableHardSphere` |
| Collision selection | No-Time-Counter | see §3.3 — nothing to configure |
| Exit condition | uniform drifting Maxwellian | `0/boundaryU`, `0/boundaryT`, `0/boundaryNumberDensity_Ar` |
| Sampling domain | `0 ≤ X/D ≤ 10` | `geometry.x_max_over_D`, `post.centerline_x_over_D_max` |
| Contour levels | `n/n0` = 0.1, 0.01, 0.001 | `post.contour_levels` |
| Resolution | `Δx = Δy = λ0` | `mesh.target_cell_over_mfp: 1.0` — see §4 |
| Time step | `Δt = t0` | see §5 |
| Transient | at least `10⁴ t0` before sampling | `dsmc.transient_collision_times: 10000` |
| Reference state for `λ0`, `t0` | the **Kn = 0.01** exit properties | `resolution.reference_knudsen: 0.01` |
| Reported max centreline density error | 0.07% / 3.63% / 12.03% | reference values only — §8 |

---

## 2. The derivation chain

Nothing physical is stored as a number that could have been derived. Four inputs
produce the entire case:

```
[PAPER]      Kn, D          →  λ0 = Kn·D                        [DERIVED]
[ASSUMPTION] T0             →  n0 from λ0                       [DERIVED]
[PAPER]      S0             →  U0 = S0·√(2RT0)                  [DERIVED]
                            →  particle weight from n0 and the cell volume
                            →  Δt from the cell size and U0 + 3σ
```

`cases/cai2012/manifest.yaml` records every one of these per case, so a results
directory can always be traced back to the Knudsen number that produced it. A
case named `Kn0p1` **cannot** run another case's density, because no density is
stored anywhere; `study.yaml` is forbidden from overriding `exit.knudsen`.

### 2.1 `Kn → n0`

For a VHS gas, carrying the speed-dependent cross-section through the
equilibrium collision rate and dividing the mean thermal speed by it gives

```
λ_VHS = (T/T_ref)^(ω − 1/2) / (√2 · π · d_ref² · n)                        (G1)
```

The `Γ(5/2 − ω)` factors that appear in the collision rate and in the mean
relative speed cancel exactly, which has two consequences worth stating:

* at `T = T_ref` the VHS mean free path **is** the hard-sphere one,
  `λ = 1/(√2 π d² n)` — the same expression
  `plumetools.markelov1999.resolution` already uses;
* away from `T_ref` the two differ by `(T/T_ref)^(ω−1/2)`, which for argon at
  300 K against 273 K is **1.0229**. Small, but not round-off.

Inverting (G1) gives the density each case needs:

```
n0 = (T0/T_ref)^(ω − 1/2) / (√2 · π · d_ref² · λ0)                         (G2)
```

| Case | `λ0 = Kn·D` | `n0` |
|---|---|---|
| Kn = 100 | 20 m | 6.620 × 10¹⁶ m⁻³ |
| Kn = 0.1 | 0.02 m | 6.620 × 10¹⁹ m⁻³ |
| Kn = 0.01 | 0.002 m | 6.620 × 10²⁰ m⁻³ |

**[ASSUMPTION]** Cai does not state which mean-free-path convention GRASP used.
This package defaults to the full VHS form (G1); `gas.mean_free_path_convention:
hard_sphere` drops the temperature factor. The choice moves every density by
2.3%, uniformly across all three cases, so it cannot change the *trend* the
study exists to test.

### 2.2 `S0 → U0`

Cai's definition exactly, with `R = k_B/m` for argon (208.24 J/(kg K)):

```
U0 = S0 · √(2 R T0) = 2 × 353.48 = 706.95 m/s
```

No γ and no gas-dynamic relation is involved. `β0 = 1/(2RT0)` is Cai's
normalising constant, so `U√β0` is simply the local speed ratio.

---

## 3. Assumptions and implementation choices

### 3.1 Exit temperature — 300 K **[ASSUMPTION]**

Cai does not appear to print an absolute `T0`. The quantities he validates are
normalised — `n/n0`, `U₁√β0`, `T/T0` — and depend on `S0` and `Kn`, not on `T0`.
300 K is recorded in `case.yaml` as an explicit assumption rather than presented
as a paper value.

Changing it rescales `U0` and `n0` together and leaves every validation curve
unchanged, which is the test of whether an assumption is load-bearing: this one
is not.

### 3.2 Argon VHS constants **[ASSUMPTION / existing repository model]**

```
mass     = 6.63e-26 kg
diameter = 4.17e-10 m   at   T_ref = 273 K
omega    = 0.74
internal degrees of freedom = 0
```

These are the values `cases/3d-inflow` already runs argon with, kept because the
case specification asks to start from the repository's model rather than
introduce new molecular parameters. One difference from the literature is worth
recording: **Bird's tabulated argon has ω = 0.81**, not 0.74, at the same
4.17 × 10⁻¹⁰ m / 273 K reference. 0.74 is the OpenFOAM tutorial's N₂ exponent.

ω enters the study in exactly one place, the `(T0/T_ref)^(ω−1/2)` factor in (G2).
Between 0.74 and 0.81 that factor moves from 1.0229 to 1.0269 — 0.4% on every
density, identically across the three cases. It is recorded here rather than
silently changed because the repository's argon is the stated starting point.

`diameter_m` and `t_ref_K` are a **pair**: the cross-section is `π d_ref²` *at*
`T_ref` and scales from there, so 4.17 × 10⁻¹⁰ m at 273 K and the same number at
300 K are different gases. The configuration keeps both.

### 3.3 No-Time-Counter — nothing to configure

Cai states NTC. There is no setting for it: NTC is the only collision-selection
algorithm `DSMCCloud::collisions()` implements in standard `dsmcFoam`, so the
requirement is met by construction. This is recorded rather than left as an
unexplained absence in the generated dictionary.

### 3.4 Collisions stay on at Kn = 100

`NoBinaryCollision` exists and would make the Kn = 100 case reproduce the
collisionless solution by fiat. That proves nothing. The validation claim is that
DSMC *approaches* the collisionless solution as the density falls, and that
requires the collision machinery to be running and finding almost nothing to do —
which is exactly what `log.dsmcFoam` reports (`No collisions` on most steps at
Kn = 100).

`dsmc.collisions_enabled: false` exists as a diagnostic, and the generated
`dsmcProperties` says **COLLISIONLESS** in capitals when it is used.

### 3.5 Domain extents **[ASSUMPTION]**

`0 ≤ X/D ≤ 10`, `|Y|/D ≤ 10`, `|Z|/D ≤ 10`. The `X` range is Cai's sampling
range. The lateral extent is an assumption: the plume must reach a boundary that
deletes molecules rather than one that reflects them, and 10 D puts that boundary
where `n/n0 < 10⁻³` everywhere along it — outside the lowest contour Cai draws.

Every extent is expressed in nozzle diameters, so `D` is the only length in the
configuration that sets a scale.

### 3.6 The boundaries

| Patch | Geometric type | Behaviour |
|---|---|---|
| `nozzle` | `patch` | injects; the disk `r ≤ R0` at `x = 0` |
| `upstreamVacuum` | `patch` | the rest of `x = 0`; **deletes** |
| `vacuum` | `patch` | the other five faces; **deletes** |

The types are load-bearing, not labels. `particle::hitBoundaryFace` finds no
handler for a plain `patch` and sets `keepParticle = false` — the molecule is
deleted, which is what "expands into vacuum" means. A `wall` would run
`DSMCParcel::hitWallPatch` and reflect.

`upstreamVacuum` is deliberately **not** a symmetry plane. The configuration is
symmetric about the axis, not about `x = 0`; a molecule that scatters back
through the exit plane has left the domain, and calling that boundary "symmetry"
would reflect it back in.

This layout is only possible because the case runs **`plumeFieldInflow`**, which
takes an explicit patch list. Stock `FreeStream` injects on every `patch`-type
boundary with no selection list, so all six boundaries would become inlets; the
only way to stop that is to declare them walls, which replaces the vacuum with a
reflecting box. See [`docs/solver-compatibility.md`](solver-compatibility.md).

---

## 4. The mesh — the main departure from Cai

### 4.1 What a literal translation would cost

Cai's grid is uniform at `Δx = λ0(Kn = 0.01) = 2 mm`. Translated literally into
the 3-D Cartesian box above:

```
1000 × 2000 × 2000 = 4.0 × 10⁹ cells
```

and at least as many parcels — several hundred gigabytes of parcel storage alone,
before the mesh is counted. `plumetools.cai2012.mesh.uniform_cost` computes this
for whatever configuration is loaded, and **`./Allmesh` prints it before writing
anything**. The impracticality is a reported result, not something discovered
halfway through a mesh.

Cai can afford his grid because his solver is **axisymmetric**: the mesh is a 2-D
`(x, r)` sheet. That is the whole difference.

### 4.2 What is built instead

A tensor-product multi-block box: a uniform fine **core** wrapped around the
nozzle and the near plume, and one geometrically expanding segment on each
outward side — 2 × 3 × 3 = 18 blocks.

```
  z
  ^   +--------+---------------------------+
  |   |        |                           |   expanding
  |   +--------+---------------------------+
  |   | core   |  expanding in x           |
--+---#========+===========================+--> x
  |   | (fine, |                           |
  |   +--------+---------------------------+
  |   |  uniform)                          |   expanding
      +--------+---------------------------+
     x=0    core_x                       x_max
```

Two properties are enforced rather than hoped for:

* **the exit disk lies entirely inside the uniform core**, so every injecting
  face is the same size (`mesh.core_half_over_D > 0.5` is a configuration
  error);
* **the first cell of each expanding segment matches the core cell**, because
  the expansion ratio is *solved for* rather than typed in. A step in cell size
  at the core boundary would appear in the sampled density as a feature that
  looks physical. `mesh.outer_expansion` caps the total ratio; if the cap binds,
  the mesh reports the jump instead of hiding it.

### 4.3 Resolution actually achieved

The core cell is `min(target_cell_over_mfp · λ0, max_core_cell_over_D · D)`,
then coarsened if `mesh.max_cells` requires it — which is reported every time.

The geometric cap is what makes Kn = 100 meshable at all: there `λ0` is 20 m, a
hundred nozzle diameters, and the collision criterion alone would allow a cell
larger than the nozzle. `max_core_cell_over_D: 0.05` puts 20 cells across the
exit diameter.

| Case | `λ0` | core cell | `Δx/λ0` | cells | meets Cai's `Δx/λ0 = 1`? |
|---|---|---|---|---|---|
| Kn = 100 | 20 m | 10 mm | 0.0005 | 1 310 720 | yes, trivially |
| Kn = 0.1 | 20 mm | 10 mm | 0.5 | 1 310 720 | yes |
| Kn = 0.01 | 2 mm | 5.39 mm | **2.69** | 3 726 000 | **no** |

**Kn = 0.01 is the one case where Cai's resolution criterion is not met.** A
2 mm core in this domain is roughly 3 × 10⁷ cells, past `mesh.max_cells`. The
achieved ratio is in `manifest.yaml`, the DSMC checks warn about it on every
mesh and every run, and the effect is a biased collision rate in the densest
region — the region where that case's departure from the collisionless solution
comes from. **Its density error should therefore be read as a lower bound on
what a properly resolved run would give, not as a measurement.**

Raising `mesh.max_cells` and re-running is a configuration change, not a code
change.

### 4.4 The nozzle patch is a staircase

A circle has no exact representation on a Cartesian grid. The `nozzle` patch is
the set of `x = 0` faces whose centres satisfy `y² + z² ≤ R0²`, selected by
`topoSet`'s `cylinderToFace` and promoted by `createPatch`.

`plumetools.cai2012.mesh.predicted_nozzle_faces` applies the same rule in Python
before the mesh exists; `./Allmesh` compares the prediction with the built mesh
and fails if they disagree, because a mesh that is not the mesh that was checked
is not worth solving on.

| Case | faces | area error vs `πR0²` |
|---|---|---|
| Kn = 100, Kn = 0.1 | 316 | +0.59% |
| Kn = 0.01 | 1092 | −0.24% |

The injected molecule flow is proportional to the patch area, so this error is a
direct multiplier on the whole solution. `checks.max_nozzle_area_error` (5%) is a
hard failure, and the OpenFOAM tests assert that the error *falls* under mesh
refinement rather than merely being small.

---

## 5. The time step — Cai's `Δt/t0 = 1` is not reconstructible

Cai states `Δt/t0 = 1`, with `t0` referred to the Kn = 0.01 exit properties, but
does not define `t0`.

Taking the most common reading — the time to travel one mean free path at the
mean thermal speed, which is also `1/ν` for the collision rate (G1) was derived
from:

```
t0 = λ_ref / c̄ = 0.002 / 398.9 = 5.014 × 10⁻⁶ s        [ASSUMPTION]
```

At Cai's own cell size (`Δx = λ_ref = 2 mm`) a molecule at `U0 + 3σ ≈ 1457 m/s`
covers **3.6 cells** in one `t0`. A DSMC step that lets molecules skip cells
without being offered a collision partner in them is not a valid step, and this
implementation fails a case on it (`checks.max_courant`).

So the time step is **derived from an explicit Courant target** instead:

```
Δt = courant_target · Δx_min / (U0 + 3σ)          courant_target = 0.2
```

and Cai's value, and the Courant number it would give on the mesh in hand, are
computed and printed for every case. That is the honest position: the criterion
is represented, not claimed to be reproduced.

| Case | `Δt` used | Courant | Cai's `Δt = t0` would give |
|---|---|---|---|
| Kn = 100, Kn = 0.1 | 1.373 × 10⁻⁶ s | 0.20 | 0.73 |
| Kn = 0.01 | 7.322 × 10⁻⁷ s | 0.20 | 1.37 (fails) |

Three standard deviations rather than the mean thermal speed, because the mean
would let the fast tail cross several cells per step unnoticed — and it is the
tail that breaks the collision sampling.

### 5.1 The transient

`dsmc.transient_basis` selects what sets the start of sampling:

* **`cai`** (the `baseCase` default) — `10⁴ t0` = 0.0501 s, Cai's own statement.
  That is 18 transits of the domain at `U0`, about 36 000 timesteps.
* **`transits`** — a configurable number of domain transits.

`cases/cai2012/study.yaml` overrides the three cases to `transits` with 2
transits of transient and 1.5 of sampling. The flow is a supersonic beam into a
vacuum with no recirculation and nothing slow to equilibrate, so it is steady
quickly: **measured** on the Kn = 100 case, the parcel count reaches its
predicted steady value of 1.155 × 10⁶ at 1.32 transits and is flat thereafter.
1.5 transits of sampling is roughly 3000 timesteps of averaging over 1.2 × 10⁶
parcels.

The override is in `study.yaml` where it can be seen rather than in the template,
both durations are printed whichever is used, and `manifest.yaml` records the
resulting `average_start_s`.

#### Continuing a run

`system/controlDict` carries `startFrom latestTime`, so re-running `./Allrun`
resumes the case from its newest time directory rather than restarting it, and
`fieldAverage` continues accumulating from the state stored in that directory's
`uniform/functionObjects/`. `./Allrun` skips the two steps that would otherwise
destroy what it is resuming from — `dsmcInitialise` and `decomposePar -force` —
and `./Allmesh` refuses to re-mesh a case that has results. Details in
[`cases/cai2012/README.md`](../cases/cai2012/README.md).

One OpenFOAM detail forced a change here. `writeControl runTime` computes its
write index as `(value − startTime)/writeInterval`, and on a resume `startTime`
is the **resume point** — so the schedule begins again there, and a resumed leg
shorter than one write interval produces no time directory at all. Measured:
resuming `Kn100` at `t = 8.487 × 10⁻³` s with a `1.980 × 10⁻³` s interval ran the
remaining 1030 steps to `endTime` and wrote nothing.

The generated `controlDict` therefore uses `writeControl timeStep`. The step
index is **global** — each time directory stores it in `uniform/time` and it
continues across restarts — and the end time is snapped to a whole number of
write intervals, so a write lands exactly on `endTime` from any resume point.

#### What the short transient costs, measured

"Steady" in the parcel count is not the same as "steady everywhere". The parcel
count is dominated by the near field, and it does plateau at 1.32 transits. But
the **far** centreline is fed by the *slow tail* of the exit distribution, and
those molecules are still in flight.

On the far centreline the visible cone is narrow — 2.95° at `X/D = 10` — so the
axial velocity distribution there is essentially the exit Maxwellian,
`U0 = 707 m/s` with `σ = 250 m/s`, and a steady-state density weights it as
`f(v_x)` (the exit *flux* goes as `v_x f`, and slower molecules linger
proportionally longer). A molecule reaching `x` by the time sampling starts at
`t_start` needs `v_x > x/t_start`:

| `X/D` | `x` | `v_min = x/t_start` | in σ | missing | measured deficit |
|---|---|---|---|---|---|
| 3 | 0.6 m | 106 m/s | −2.4 σ | ~0.8% | ~1% |
| 7 | 1.4 m | 247 m/s | −1.8 σ | ~3% | 3.2% |
| 10 | 2.0 m | 354 m/s | −1.4 σ | ~8% | 5–6% |

The measured Kn = 100 centreline (`Cases/Kn100/results/centerline.csv`) agrees
with the analytical solution to **~1% out to X/D ≈ 3** and then falls
progressively below it, reaching −5.6% at `X/D = 9.7`. That is not a collisional
effect and not a mesh effect: it is the slow tail of the exit distribution not
having arrived yet, and it is exactly what Cai's `10⁴ t0` transient exists to
avoid. (The figures in the table above shrank by a third when the sampling
window was extended from 1.0 to the configured 1.5 domain transits by simply
re-running `./Allrun` — which is the resume workflow doing its job.)

**So the far-field half of the Kn = 100 comparison is transient-limited, not
converged.** Setting `dsmc.transient_basis: cai` in `study.yaml` fixes it at
about three times the run time. This is recorded rather than trimmed out of the
comparison range.

---

## 6. The analytical comparison solution

`plumetools/cai2012/analytical.py` implements Cai's collisionless solution for a
circular exit. **The closed forms there are re-derived from the underlying
free-molecular integral, not transcribed from the paper.** Equation numbers are
cited as pointers to where the same results appear in Cai (Eq. 5 for the density
field, Eqs. 18–21 for the centreline), not as a claim that the algebra is
character-for-character his. Nothing is digitised from a plot.

With `t0 = cos θmax = x/√(x² + R0²)`, `E0 = e^(−S0²(1−t0²))`,
`A = 1 + erf S0` and `A0 = 1 + erf(S0 t0)`:

```
n/n0                    = ½ [A − t0·E0·A0]

(n/n0)·U1√β0            = e^(−S0²)(1−t0²)/(2√π) + (S0/2)[A − t0³·E0·A0]

(n/n0)·β0⟨v²⟩           = S0·e^(−S0²)(1−t0²)/(2√π)
                          + 2(3/8 + S0²/4)·A
                          − 2(3/8·t0 + S0²/4·t0³)·E0·A0

T/T0                    = (2/3)[β0⟨v²⟩ − (U1√β0)²]
```

The derivation is written out in the module docstring. Every limit that can be
established independently is asserted in the tests:

| Limit | Value | Note |
|---|---|---|
| `x → 0`, any `S0` | `n/n0 → ½(1 + erf S0)` | **not exactly 1** |
| `x → 0`, `S0 = 2` | `n/n0 = 0.99766` | 0.23% below 1 |
| `x → 0`, `S0 = 0` | `n/n0 = ½` | only the outgoing half is present |
| `x → 0`, `S0 = 0` | `U√β0 = 1/√π` | `⟨|v_x|⟩ = √(2kT/πm)` |
| `x → 0`, `S0 = 0` | `T/T0 = 1 − 2/(3π) = 0.78779` | half-Maxwellian |
| `x → ∞` | `n/n0 → 0` like `x⁻²` | point source |
| `x → ∞`, `S0 = 2` | `U√β0 → 2.44` | **rises** downstream |

The exit-plane density deserves a note, because it is the one place where the
obvious expectation is wrong. `n/n0` does **not** approach 1 as `x → 0`: only
the forward-moving half of the exit distribution is present in the plume, so the
limit is `½(1 + erf S0)`. At Cai's `S0 = 2` that is 0.99766 — indistinguishable
from 1 on a plot, and exactly 1/2 with no drift.

The rising centreline velocity is the signature of a collisionless plume: the
visible cone narrows downstream and the slow molecules leave it first, and there
are no collisions available to bring the mean back to the bulk speed.

The full field (Cai Eq. 5) is evaluated by Gauss–Legendre quadrature over the
exit disk. It is the same solution as the closed form, and the test that they
agree on the axis to 1 × 10⁻¹⁰ is what establishes that.

---

## 7. What is measured, and how

Everything comes from the **time-averaged** fields `fieldAverage` writes:
`rhoNMean`, `rhoMMean`, `momentumMean`, `linearKEMean`. A single-timestep `rhoN`
is one step's worth of parcels — `DSMCCloud::resetFields` zeroes it every step —
and comparing that with an analytical curve measures shot noise.

`dsmcFoam` writes moments, not derived quantities, so:

```
u      = momentumMean / rhoMMean
⟨v²⟩   = 2 · linearKEMean / rhoMMean
3 R T  = ⟨v²⟩ − |u|²
```

which is the same definition the analytical solution uses, so the two
temperatures are the same quantity rather than two things with the same name.

**Centreline.** Cells within `post.centerline_radius_over_D` (0.10 D) of the axis
are binned in `x`. The density is a volume-weighted mean — total molecules over
total volume — and the velocity and temperature are weighted by molecule count.
On a graded mesh an unweighted mean lets one huge outer cell outvote a hundred
core cells, and a cell with no molecules carries no information about the
velocity there.

A DSMC "centreline" is unavoidably an average over a **tube**, and the plume is
peaked on the axis, so that average sits below the on-axis analytical value. The
bias is not small — 6.8% in the density at `X/D = 1` for a 0.25 D tube, 1.1% at
0.10 D — and it is the same size as the collisional departure this study exists
to measure. So the analytical solution is put through **the same estimator**:
`plumetools.cai2012.analytical.moments` is evaluated at the same cell centres and
averaged with the same weights. That is what `metrics.yaml` scores against; the
on-axis curves, which are what Cai plots, are carried alongside as
`*_analytical_axis` and are what the figures show.

**Density plane.** The two cell layers straddling the centre plane, at the cell
centres, with no resampling. A graded Cartesian mesh has no natural grid to
interpolate onto, and a resampled field would hide where the mesh actually is —
which, on a contour spanning three decades, is the thing most worth being able to
see. `matplotlib.tricontour` draws it directly.

Cell centres and volumes come from OpenFOAM's own function objects
(`postProcess -func writeCellCentres` / `writeCellVolumes`), which `./Allpost`
runs.

---

## 8. How to read the numbers

Cai reports these maximum centreline analytical-vs-DSMC density errors:

| Kn | Cai |
|---|---|
| 100 | 0.07% |
| 0.1 | 3.63% |
| 0.01 | 12.03% |

**These are reference values, not pass/fail tolerances.** Reproducing them
exactly is not the claim and would not be meaningful: different solver, different
mesh topology, different resolution, different molecular constants, and a
different definition of "maximum" (over which stations, at what sampling noise).

What is being validated is the **trend**:

```
Kn = 100    agreement should be extremely strong -- this is the collisionless limit
Kn = 0.1    a visible collisional departure
Kn = 0.01   the largest departure from the collisionless model
```

Two things will systematically bias this implementation's numbers relative to
Cai's, both upward:

* the **Kn = 0.01 mesh is 2.7 mean free paths per cell** (§4.3), which biases the
  collision rate low in the densest region;
* the **staircase nozzle** carries a sub-percent area error (§4.4), which is a
  direct multiplier on the injected flow and therefore on `n/n0` near the exit.

DSMC results are random variables. Cai puts the worst-case centreline scatter at
roughly 1%. `results/metrics.yaml` reports maximum, mean and RMS relative error
per case, and the maximum is the noisiest of the three — the mean and RMS are the
ones to compare across cases.

---

## 8a. What was actually measured

Two of the three cases have been run to completion at the settings in
`study.yaml` (2 domain transits of transient, 1.5 of sampling, 12 ranks).
`Kn = 0.01` is generated and meshed but **has not been run** — 3.7 × 10⁶ cells
and ~7.5 × 10⁶ parcels over 13 523 steps is several hours on this machine.

Centreline relative error against the tube-matched analytical solution:

| Case | density max | density mean | density RMS | velocity max | temperature max | Cai's reported max density |
|---|---|---|---|---|---|---|
| Kn = 100 | 5.59% | **1.07%** | 1.46% | 1.50% | 4.16% | 0.07% |
| Kn = 0.1 | 17.62% | **6.81%** | 8.72% | 3.49% | 54.90% | 3.63% |
| Kn = 0.01 | — not run — | | | | | 12.03% |

**The trend is the expected one**: the departure from the collisionless solution
grows by a factor of 6.4 in the mean density error from Kn = 100 to Kn = 0.1, and
by a factor of 22 in the temperature. The magnitudes are larger than Cai's, and
§8 says why they would be; the Kn = 100 figure is additionally inflated by the
short transient (§5.1) — over `X/D ≤ 3`, where the transient has no effect, the
Kn = 100 density agrees with the analytical solution to about **1%**.

The Kn = 0.1 departure has the right **sign** in all three quantities, which is a
stronger statement than the magnitudes:

* **density** is *above* the collisionless solution near the exit (+4% at
  `X/D = 0.7`) and *below* it downstream (−18% at `X/D = 9.7`) — collisions
  collimate the near plume and then deplete the far centreline;
* **velocity** is 2–3% *higher* downstream — collisions convert thermal energy
  into directed motion, which free-molecular flow cannot do;
* **temperature** is up to 55% *lower* — a collisional expansion keeps cooling,
  where the collisionless solution freezes at `T/T0 → 0.28`.

The `Kn = 100` density contour reproduces all three of Cai's levels (0.1, 0.01,
0.001) to within the sampling noise. The `Kn = 0.1` contour visibly pulls in from
the collisionless one at the 0.001 level, closing near `X/D ≈ 8` where the
analytical contour reaches ≈ 10 — the same depletion the centreline shows.

Reproduce with:

```bash
cd cases/cai2012
./generate_cases.py && ./AllmeshCases && ./AllrunCases && ./AllpostCases
```

---

## 9. Differences from Cai's setup, collected

| | Cai 2012 | This implementation |
|---|---|---|
| Solver | GRASP, axisymmetric DSMC | OpenFOAM v2512 `dsmcFoam`, 3-D |
| Geometry | axisymmetric `(x, r)` | full 3-D Cartesian, no symmetry assumed |
| Mesh | uniform, `Δx = λ0` | graded multi-block, uniform core + expansion |
| `Δx/λ0` | 1 | 0.0005 / 0.5 / 2.69 |
| Nozzle exit | exact circle | staircase, area within 0.6% |
| `Δt` | `t0` | Courant-limited; `t0` reported alongside |
| Transient | `10⁴ t0` | 2 domain transits (`study.yaml` override) |
| Collision selection | NTC | NTC — the only one `dsmcFoam` implements |
| VHS argon | not stated | `d = 4.17e-10 m`, `ω = 0.74`, `T_ref = 273 K` |
| `T0` | not stated | 300 K, assumed |
| Particle weighting | not stated | one global weight, 20 parcels in the exit cell |

---

## 10. Known limitations

1. **Kn = 0.01 is under-resolved** against Cai's own criterion (§4.3). This is
   the single largest caveat on the study and it affects the case the whole
   comparison is most interested in.
2. **The transient is shorter than Cai's** (§5.1) — 2 domain transits against
   his 18-equivalent — and it **measurably limits the far field**. Beyond
   `X/D ≈ 5` the sampled centreline density falls progressively below the
   analytical solution, reaching −5.6% at `X/D = 9.7`, because the slow tail of
   the exit distribution has not arrived. The arithmetic is in §5.1 and matches
   the measurement. `dsmc.transient_basis: cai` fixes it at ~3x the run time.
3. **`T0` is assumed** (§3.1). Harmless for the normalised comparison, but it
   means no dimensional quantity here can be compared with the paper.
4. **ω = 0.74 is the repository's argon, not Bird's** (§3.2). A 0.4% uniform
   shift in every density.
5. **The mean-free-path convention is assumed** (§2.1). A 2.3% uniform shift.
6. **Off-axis analytical velocity and temperature are not implemented.** Only the
   density field is available off the axis (Cai Eq. 5), which is what the contour
   comparison needs; the velocity and temperature comparisons are centreline
   only, as §9 of the case specification asks.
7. **Nothing here has been compared with an independent DSMC code.** The
   comparison is against the analytical collisionless solution, which is exact
   for Kn → ∞ and is *supposed* to disagree at Kn = 0.01.

---

## 11. Related documents

* [`cases/cai2012/README.md`](../cases/cai2012/README.md) — how to run it
* [`docs/plume-field-inflow.md`](plume-field-inflow.md) — the custom inflow model
* [`docs/solver-compatibility.md`](solver-compatibility.md) — why stock
  `FreeStream` cannot express this case
* [`docs/markelov1999-case.md`](markelov1999-case.md) — the other validation
  family, and the design this one follows

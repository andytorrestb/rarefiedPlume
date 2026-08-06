# The AIAA 99-3455 case family

Lumpkin, F. E., Stewart, B. D., and Markelov, G. N., *Study of 3D Rarefied Flow
on a Flat Plate in the Wake of a Cylinder*, AIAA 99-3455, 1999.

This document separates **what the paper prints** from **what this
implementation had to decide**, and gives the reasoning for each decision. The
running instructions are in
[`cases/markelov1999/README.md`](../cases/markelov1999/README.md).

---

## 1. Values printed by the paper

| Quantity | Paper | SI |
|---|---|---|
| Orifice diameter | 0.8255 mm | 8.255 × 10⁻⁴ m (radius 4.1275 × 10⁻⁴ m) |
| Cylinder diameter | 6 in | 0.1524 m (radius 0.0762 m) |
| Cylinder length | 18 in | 0.4572 m |
| Cylinder centre, x | 11.75 in | 0.29845 m |
| Plate width, y | 6 in | 0.1524 m |
| Plate height, z | 15 in | 0.3810 m |
| Cylinder base → plate face | 6 in | 0.1524 m |
| Reservoir pressures | 5, 25, 100, 475 psi | 34 473.79 / 172 368.93 / 689 475.73 / 3 275 009.71 Pa |
| Gas | nitrogen | γ = 1.4, M = 28.0134 g/mol |
| Collision model | VHS + Larsen-Borgnakke | `LarsenBorgnakkeVariableHardSphere` |
| Rotational collision number | Z_R = 5 | `relaxationCollisionNumber 5` |

Derived, never written as a literal:

```
cylinder upstream  x = 0.29845 − 0.0762            = 0.22225 m
cylinder downstream x = 0.29845 + 0.0762           = 0.37465 m
plate upstream face x = 0.37465 + 0.1524           = 0.52705 m
plate downstream    x = 0.52705 + thickness        = 0.53975 m
```

`./Allmesh` prints this arithmetic, and
`plumetools.markelov1999.verify` measures every one of them back out of the
generated mesh rather than reading them from the config.

---

## 2. The source-flow model

Model identifier: **`markelov1999_axisymmetric`**, in
`plumetools/markelov1999/sourceflow.py`.

With `θ` the angle from the `+x` plume centreline and `r` the distance from the
orifice centre:

```
P1   θ_L  = π/2 · (√((γ+1)/(γ−1)) − 1)
P2   V    = √(2γkT₀ / ((γ−1)m))
P3   f(θ) = cos(πθ / (2θ_L)) ^ ((γ + 0.41)/(γ − 1))
P4   A    = 0.5·√((γ−1)/(γ+1)) / ∫₀^θ_L sin θ · f(θ) dθ
P5   ρ(r,θ) = (2Ap₀/V²) · (2/(γ+1))^(1/(γ−1)) · (r_e/r)² · f(θ)
P6   U(r,θ) = V · r̂
```

At γ = 1.4, for nitrogen at T₀ = 300 K:

| | |
|---|---|
| θ_L | 2.276853 rad = 130.454° |
| V | 789.4849 m/s |
| angular exponent | 4.525 |
| ∫ sin θ f(θ) dθ | 0.357656358 |
| A | 0.570727014 |
| (2/(γ+1))^(1/(γ−1)) | 0.633938145 |

`ρ` is a mass density. The conversion to number density happens in exactly one
place, `n = ρ/m`, with the same `m` that P2 uses.

`θ_L > π/2` exactly when `γ < 5/3`, so for nitrogen the cone reaches 130° while a
hemispherical inflow surface spans at most 90°. The clipping in P3 therefore
never fires on this geometry — which is why it is **checked** on the meshed
patch rather than assumed.

### Differences found between the paper and the legacy code

`plumetools/sourceflow.py` is frozen bit-for-bit to reproduce the pre-refactor
results, defects included. It is not the same model. Six differences, none of
them a rounding detail:

| | Legacy | Paper-faithful |
|---|---|---|
| **SF-1** angular dependence | a *product* `f(θ)·f(φ)` of two angles in different planes (`angular_separable`) | one angle, as printed |
| **SF-2** angle convention | `θ = arccos(z/r)`, from `+z`, patched inside the angular function by an `abs(θ − π/2)` shift | measured directly from the `+x` plume axis |
| **SF-3** orifice radius | `throat_radius_m: 0.0041275` — **ten times** the paper value, so density is **100×** too large (P5 goes as `r_e²`) | 4.1275 × 10⁻⁴ m |
| **SF-4** reservoir pressure | `calculateRhoN` hard-coded 475 psi and ignored each case's own setting (SM-01); the archived sweep therefore varied nothing (AR-02) | a required argument with no default; exactly linear in `p₀` |
| **SF-5** normalisation γ | a separate `integrand_gamma` defaulting to 1.4 whatever γ was passed, *and* normalising the θ-only function while the density applied the separable one — so `A` did not normalise what it scaled (SM-05) | one γ throughout, normalising exactly the `f` that P5 applies |
| **SF-6** beyond θ_L | `f_θ` takes `abs()` of the cosine so the profile folds and rises again; `f_φ` does not, and returns **NaN** (SM-06) | exactly zero, and finite everywhere |

A seventh, not in that module but of the same kind: the reference case evaluates
an N₂ model while its solver simulates argon (SM-03). This family is nitrogen
throughout, and `plumetools.config` cross-checks `gas.species_name` against the
generated `dsmcProperties`.

The legacy functions are **not** removed. The regression suite freezes them
deliberately.

### The `0.41`

`(γ + 0.41)/(γ − 1)` is the exponent printed in the paper. The `0.41` is an
empirical constant for which the paper gives no derivation. It is exposed as
`angular.exponent_offset` so a future correction is a configuration change.

---

## 3. Assumptions and implementation choices

Everything in this section is a decision this implementation made. Each is tagged
`[ASSUMPTION]` in `baseCase/case.yaml`.

### 3.1 Inflow hemisphere radius — 0.1524 m (6 in)

The paper's physical source is the 0.8255 mm orifice. The hemisphere is a
**computational-domain** choice: the surface on which the analytical source-flow
solution is handed to the particle solver. DSMC cannot start at the orifice —
the density there is far too high for any affordable weighting, and the flow is
still collisional.

0.1524 m is the inflow arc radius the archived `wake-cylinder` cases intended
(recorded in [`cases/ARCHIVE.md`](../cases/ARCHIVE.md), finding AR-01). The
*value* is reused; that lineage's implementation is not — it divided quadrilateral
vertex sums by a literal `3.0`, putting every centroid at 4/3 of its true radius.

Geometric consistency is checked, not assumed:

```
clearance = 0.22225 − 0.1524 = 0.06985 m  (2.75 in)
```

`MarkelovGeometry.validate` refuses a configuration where this is ≤ 0, and
`plumetools.markelov1999.checks` reports the distance on every run.

**The source-flow equations use the orifice radius, never this one.** Confusing
them scales density by `(R/r_e)² ≈ 1.4 × 10⁵`, so `plumetools.config` requires
`stagnation.throat_radius_m < geometry.sphere_radius_m`.

### 3.2 Plate thickness — 0.0127 m (0.5 in)

The paper does not clearly establish a thickness. 0.5 in is thin against the 6 in
width, and thick enough to mesh as a solid body with two distinct faces rather
than as a zero-thickness baffle.

It is configurable, and it moves only the downstream face: the gap is measured to
the *upstream* face, so changing the thickness does not change the paper
dimension.

### 3.3 Domain extents

| | Extent | Why |
|---|---|---|
| x | [0, 0.9] m | The plate's downstream face is at 0.53975 m, leaving 0.36 m — 2.4 cylinder diameters — for the plate wake before particles are deleted. |
| y | [0, 0.3] m | Half domain; `y = 0` is the symmetry plane. 0.3 m is 3.9 cylinder radii, well outside the plume cone at the bodies. |
| z | [−0.35, 0.35] m | The cylinder spans ±0.2286 and the plate ±0.1905, leaving 0.12 m past the cylinder end caps. |

**Only `y = 0` is a symmetry plane.** `z = 0` is not: the cylinder end caps and
the plate edges are inside the domain, and halving `z` would delete the
three-dimensional wake this case exists to study. `plumetools.markelov1999.mesh`
refuses a domain whose `z` bounds are not strictly outside both bodies.

The `x = 0` plane outside the inflow cavity is `upstreamVacuum` — an **open
patch**, not a symmetry plane. Naming it `symmetry` would reflect back every
particle that scattered upstream.

### 3.4 Mesh resolution

Background 0.025 m → 36 × 12 × 28 = 12 096 cells. Coarse deliberately: the far
field is near-vacuum, and its cells only have to *find* the surfaces.

| surface | level | cell | resolves |
|---|---|---|---|
| inflow | 2 | 6.25 mm | 24 cells across the radius |
| cylinder | 2 | 6.25 mm | ~77 cells around the circumference |
| plate | 2 | 6.25 mm | 2 cells across the 12.7 mm thickness |

Measured against OpenFOAM v2512: **34 346 cells**, `checkMesh` **Mesh OK**, max
non-orthogonality 37.6 / mean 10.2, max skewness 0.69, max aspect ratio 3.42.

The plate at level 3 (3.125 mm) also meshes correctly — 56 274 cells, `Mesh OK`,
and `verify` measures the plate's dimensions exactly either way. Level 2 is
shipped because the plate sits where the plume is ~6× less dense than at the
cylinder, and halving the cell divides the occupancy there by eight: 3.6
particles/cell at level 2 against 0.4 at level 3. Neither meets the target;
0.4 is noise.

Only the three surfaces are refined. **No wake, shock or gap refinement regions**
— `snappyHexMeshDict` carries an empty `refinementRegions` block as the place to
add them.

### 3.5 Molecular properties

The paper prints the gas, the collision-model family and Z_R, but not the VHS
coefficients modern OpenFOAM needs. These come from OpenFOAM's own N₂ model, and
are recorded rather than claimed as a reproduction:

| | Value | Source |
|---|---|---|
| `diameter` | 4.17 × 10⁻¹⁰ m | `tutorials/discreteMethods/dsmcFoam/wedge15Ma5` |
| `omega` | 0.74 | as above |
| `internalDegreesOfFreedom` | 2 | as above — rotation only, no vibration |
| `Tref` | 273 K | as above |
| `mass` | 4.651735 × 10⁻²⁶ kg | **derived** as `M/(1000·N_A)` |

The mass is derived rather than taken from the tutorial's rounded 46.5 × 10⁻²⁷
kg (0.04% away) so the analytical model and the solver use one `m`. Two masses
that nearly agree are harder to debug than one.

`gas.gamma` and `internal_degrees_of_freedom` are two statements of the same
physics, so `plumetools.config` refuses a configuration where they disagree —
finding SM-03 in a new costume.

### 3.6 Stagnation temperature — 300 K

Not recorded here as a paper value. Used by P2 and written as the inflow
temperature.

### 3.7 Time step — 2 × 10⁻⁷ s

At `V + 3σ_thermal` ≈ 1685 m/s a particle crosses **0.054** of the smallest cell
per step. `plumetools.markelov1999.checks` makes crossing more than one cell a
**hard error**: particles that skip a cell are never offered a collision partner
in it, which invalidates the collision sampling.

`end_time_s` = 4 × 10⁻³ s is roughly 3.5 domain transits at 789 m/s, with
`fieldAverage` starting at 2 × 10⁻³ s. **These are assumptions, not a converged
result.** Check the particle count in `log.dsmcFoam` before reading a pressure.

### 3.8 Particle weighting

Derived per case so `resolution.sizing_region` reaches
`target_particles_per_cell`, and written back into each generated `case.yaml` so
it is auditable rather than implicit. Because density is linear in `p₀`, the
weight is too — and the total particle count is therefore the same for all four
pressures, ~3.4 × 10⁶.

**Standard `dsmcFoam` supports no variable weighting.** `DSMCCloud::nParticle_`
is a single scalar. There is no radial or adaptive scheme to configure, so the
occupancy shortfall away from the sizing region is a limitation to report, not a
parameter to tune.

### 3.9 Pressure-averaging windows

±15° in azimuth about the windward (180°) and leeward (0°) generators, ±0.0381 m
(1.5 in) about mid-span. A window rather than a face because a one-face reading
depends on where snappyHexMesh happened to put that face, and would change with
the refinement level rather than with the physics.

---

## 4. Inflow flux verification

The solver samples a **half-range drifting Maxwellian**, not `ρUA`. The number
flux it accumulates is Bird eqn 4.22:

```
F = n · c_mp / (2√π) · [ exp(−s²) + √π · s · (1 + erf s) ],   s = (U·n̂)/c_mp
```

The drift term `n·U·n̂` is the `s → ∞` limit. Verifying against `ρUA` would not be
verifying the boundary condition that is applied, and at `s = 0` it reports no
inflow at all while a real surface effuses `n·c_mp/(2√π)`.

Every moment is computed three ways:

1. closed-form half-range moments;
2. deterministic Gauss-Legendre quadrature, using no closed form;
3. Monte Carlo over the drifting Maxwellian.

Measured on the shipped 5 psi case:

```
quantity           analytical   from 0/ fields   rel. error
number [1/s]     1.285506e+20     1.285506e+20     8.6e-12
mass [kg/s]      5.979834e-06     5.979834e-06     8.6e-12
momentum x [N]   3.876458e-03     3.876458e-03     8.1e-12
momentum y [N]   2.102529e-03     2.102529e-03     1.2e-11
momentum z [N]   8.622568e-08     8.622568e-08     2.0e-10
energy [W]       3.725918e+00     3.725918e+00     8.4e-12

independent quadrature agrees to 6.0e-15 (no closed form used)
drift-only n·U·A gives 99.9% of the true half-range flux
```

The field round-trip error is the `%.10g` formatting of the written values, not a
model discrepancy. The drift-only figure is 99.9% here because this plume is
highly directed (`s ≈ 1.87`) — small, but it is a property of this case rather
than a reason to skip the thermal terms.

---

## 5. DSMC quality checks

Most **warn**: a cell coarser than the mean free path degrades the collision
statistics, but that is a judgement about accuracy, and silently redesigning
someone's mesh is worse than telling them.

Five conditions are **errors**:

* non-positive or non-finite density or temperature;
* a NaN in the generated inflow fields;
* the inflow surface intersecting the cylinder;
* a time step letting a particle cross more than one of the smallest cells;
* a required patch missing from the mesh.

Every result — pass, warn or fail — is recorded in `case-summary.json`.

---

## 6. Known limitations

* **No validation.** No comparison with DAC or experiment. `AllpostCases` says so
  in its own output.
* **The 20 particles/cell target is met in one region only.** See §3.8.
* **The cylinder cell is marginal against the mean free path at 475 psi.**
  Measured: the ratio is **1.36** there — a 6.25 mm cell against a 4.6 mm mean
  free path — so `checks` warns. It falls to 0.42 in the wake, 0.24 at the plate,
  and below 0.02 everywhere at 5 psi, so only the highest-pressure case is
  affected. One more refinement level at the cylinder fixes the ratio and
  multiplies the particle count by eight.
* **No wake or gap refinement.** Deliberately, per the scope; the configuration
  path is in place.
* **The 12 in gap cases are not implemented.** The geometry derives the plate
  position from `bodies.gap_m`, so adding them is a second `gap_dir` and no code
  change.
* **Steady state is assumed, not demonstrated.** `end_time_s` is a guess.
* **The estimator ignores the bodies.** It evaluates the undisturbed source-flow
  density, so the compressed stagnation region is denser than it reports — which
  makes the estimate conservative, and is another reason it is not a measurement.

---

## 7. Related documents

* [`cases/markelov1999/README.md`](../cases/markelov1999/README.md) — running it
* [`docs/plume-field-inflow.md`](plume-field-inflow.md) — the custom boundary model
* [`docs/solver-compatibility.md`](solver-compatibility.md) — what stock `dsmcFoam` cannot express
* [`docs/source-flow-model.md`](source-flow-model.md) — the **legacy** model, and why it is frozen
* [`cases/ARCHIVE.md`](../cases/ARCHIVE.md) — the frozen `wake-cylinder` lineage

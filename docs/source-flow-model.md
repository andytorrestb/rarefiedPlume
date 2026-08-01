# The plume source-flow model

Reconstructed from `cases/3d-inflow/processInflowData.py.legacy` (555 lines),
cross-checked against all eight copies in the 3D lineage and the full commit
history.

## Provenance — read this first

**The source of these equations could not be established.**

An exhaustive search of every comment, commit message, variable name and document
in this repository — grepping for `simons`, `boynton`, `roberts`, `south`,
`reference`, `doi`, `paper`, `et al`, `cite` — turned up exactly **one**
literature reference, in `cases/2d-wedge/mesh/generateBlockMeshDict.m`:

> This program creates a blockMeshDict file to create an axisymmetric wedge mesh
> for reproducing the isentropic expansion V&V case in Stuart and Lumpkin (2012).
>
> Jonathan S. Pitt, Aegis Aerospace, Inc., for LENS
> Created: 21 Sep 2021 · Updated: 10 Nov 2021

Plot legends elsewhere read "DAC (Stewart, Lumpkin)", so "Stuart" is a typo for
**Stewart**. That reference covers the **archived** 1d/2d-wedge
isentropic-expansion validation case and the digitised `n_nstar_radius.dat` /
`T_Tstar_radius.dat` data.

**It does not cover the 3D source-flow model below.** The limiting angle, the
`(γ + 0.41)/(γ − 1)` angular exponent, and the normalisation integral appear in
no citation anywhere. No source has been invented for them. The equations here
were read out of the code.

## Coordinate system and geometry

Measured from `cases/3d-inflow/constant/polyMesh`:

| Quantity | Measured |
|---|---|
| Cells | 99 240 tetrahedra (204 717 faces, 192 243 internal; exactly 4 faces/cell) |
| Points | 19 315 |
| Inflow patch vertex \|r\| | **exactly 0.5000000000 m** |
| Inflow face centroid \|r\| | 0.4989640 – 0.4997035 m |
| Inflow face vertex counts | all 3 |
| Inflow extent | x ∈ [−0.001, 0.4999], y, z ∈ [−0.5, 0.5] |

The inflow patch is a **hemisphere of radius exactly 0.5 m centred at the origin,
opening toward +x**. Face centroids fall 0.06–0.21% inside it purely from
faceting — the centroid of a triangle inscribed in a sphere lies inside the
sphere. That is why the hard-coded `r = 0.5` is the *correct* nominal radius and
why using each face's own `|c|` would be worse: it would inject faceting spread
into `f₃ = (r_e/r)²` as if it were physics.

**Angles as implemented** (`findSymmetryTheta`, `findSymmetryPhi`):

    θ = arccos(z / 0.5)     ∈ [0, π]      polar angle from **+z**
    φ = arctan2(y, x)       ∈ (−π, π]     azimuth from +x in the x–y plane

**The plume axis is +x**, inferred from the one-sided hemisphere, from
`revolutionAxis "x"` in `dsmcProperties`, and from the first sampling line. It is
stated nowhere in the original code. θ is therefore *not* measured from the plume
axis — see SM-04.

## Equations

| # | Equation | Implementation |
|---|---|---|
| **E1** | θ_ℓ = (π/2)·(√((γ+1)/(γ−1)) − 1) | `limiting_angle` |
| **E2** | v_ℓ = √(2γ k_B T₀ / ((γ−1) m)) | `limiting_velocity` |
| **E3** | f_A(θ) = cos^n(πθ / 2θ_ℓ), n = (γ+0.41)/(γ−1) | `angular_theta_only` |
| **E4** | f(θ,φ) = \|cos(π\|θ−π/2\| / 2θ_ℓ)\|^n · cos^n(πφ / 2θ_ℓ) | `angular_separable` |
| **E5** | A = ½√((γ−1)/(γ+1)) / ∫₀^θℓ sin x · f_A(1.4, x) dx | `normalization_coefficient` |
| **E6** | n = (2A p₀ / v_ℓ²) · (2/(γ+1))^(1/(γ−1)) · (r_e/r)² · f(θ,φ) · (1000 N_A / M_w) | `number_density` |
| **E7** | U = (v_ℓ cosφ sinθ, v_ℓ sinφ sinθ, v_ℓ cosθ) | `velocity` |
| **E8** | T = T₀, written as the vector `(T 0 0)` | — |

**E6 units.** `2Ap₀/v_ℓ²` has units Pa/(m² s⁻²) = kg m⁻³, a mass density; the
middle factors are dimensionless; `1000 N_A / M_w` converts kg m⁻³ to
molecules m⁻³ for M_w in g/mol. **Dimensionally consistent.**

**E7 note.** |U| = v_ℓ on every face, and because θ and φ derive from the same
centroid, U is radially outward: **U = v_ℓ r̂**. This part of the model is
self-consistent by construction and is the one piece unaffected by SM-04.

**E8 note.** `boundaryT` is a `volVectorField` with dimensions `[0 0 0 1 0 0 0]` —
dsmcFoam+ takes per-component translational temperature, and the model writes the
y and z components as zero.

## Parameters as frozen

| Symbol | `case.yaml` key | Value | Units | Status |
|---|---|---|---|---|
| γ | `gas.gamma` | 1.4 | – | ⚠ N₂; solver simulates Ar (γ = 5/3) |
| T₀ | `stagnation.T0_K` | 300 | K | used by velocity and temperature |
| T₀ (density) | `legacy.rhoN_T0_K` | 300 | K | ⚠ separate knob — SM-02 |
| M_w | `gas.molar_mass_g_per_mol` | 28.0134 | g/mol | ⚠ N₂; argon is 39.948 |
| m | — | 4.6672e−26 | kg | derived, `M_w × AMU_TO_KG` |
| p₀ | `stagnation.p0_pa` | 3 275 009.71275 | Pa | ⚠ 475 psi — SM-01 |
| r_e | `stagnation.throat_radius_m` | 0.0041275 | m | ⚠ "exit" or throat? — SM-10 |
| r | `geometry.sphere_radius_m` | 0.5 | m | ✅ exact for this mesh |
| θ_ℓ | — | 2.27686 | rad | 130.45° |
| v_ℓ | — | **788.164** | m/s | \|U\| on all 2044 faces |
| n exponent | `angular.exponent_offset` | 0.41 → n = 4.525 | – | ⚠ uncited — SM-07 |

Constants are the **legacy** values, kept to their original precision:
k_B = 1.3806e−23, N_A = 6.023e23, amu = 1.66605e−27. Substituting CODATA values
would change every number the model produces. See `plumetools/constants.py`.

## Confidence classification

### 1. Directly established by the code

- The spherical convention: θ = `acos(z/r)` from +z, φ = `arctan2(y, x)` from +x.
- E1–E8 exactly as written, including every hard-coded literal.
- |U| = v_ℓ uniformly; direction exactly radial.
- The face normal is computed, plotted, threaded into `face_data[1]`, bound to
  `norm` inside `calculateU`, and **never read**.
- The centroid divisor is the literal `3.0`, regardless of vertex count.

### 2. Inferred from naming, geometry, or configuration

- **Plume axis = +x.** Never stated.
- **`r_e` is a nozzle exit or throat radius**, 4.1275 mm — an 8.255 mm diameter
  halved. Since `f₂ = (2/(γ+1))^(1/(γ−1))` is the isentropic *sonic* density
  ratio ρ\*/ρ₀, it reads most consistently as the **throat** radius. The variable
  name says "exit". Unresolved.
- The `|θ − π/2|` shift compensates for θ being measured from +z: it equals
  `|arcsin(z/r)|`, the elevation above the x–y plane.
- The `½√((γ−1)/(γ+1))` numerator in E5 is a flux normalisation intended to
  recover the sonic-throat mass flux.

### 3. Cannot be verified from this repository

- **The provenance of the model.** See the top of this document.
- Whether the intended gas was N₂ or Ar. The model says N₂; the solver, the
  output filename, `calcCaseParams.py` (`R = 208 # ??`, argon's gas constant) and
  the archived validation plots (γ = 1.67) all say Ar.
- Whether the separable θ·φ product in E4 was intended, or a step toward a single
  off-axis angle that was never finished.
- Whether continuum-breakdown or transitional limits were ever checked for the 3D
  plume. Nothing equivalent to the archived `calcCaseParams.py` Bird-parameter
  calculation exists here.
- Whether `dsmcAxisymmetric` (in `3d-inflow` and `iss_solar_panels`) was
  intentional on a 3D, φ-dependent mesh — HA-06.

### 4. Known defects, frozen deliberately

| ID | Defect |
|---|---|
| **SM-01** | `calculateRhoN` hard-coded `Po = 475*6894.75729`. Every case set `physical_props['Po']` — 34 500, 6 894, 16 894, 744 634 Pa — and it was discarded. n ∝ p₀, so no case's stated pressure describes its run. |
| **SM-02** | From `iss_solar_panels_2` on, velocity and temperature used the configured T₀ (700 K, then 3000 K) while density kept 300 K. v_ℓ ∝ √T₀ and n ∝ 1/v_ℓ², so mass and energy flux disagree by a factor of 10 in the resolution study. |
| **SM-03** | The model runs on N₂ while `dsmcProperties` declares `typeIdList (Ar)`, mass 6.63e−26 and `rotationalDegreesOfFreedom 0`. 788.16 m/s against 558.89 m/s. |
| **SM-04** | θ is measured from +z, not the +x plume axis. E4 is a separable product of two angles **in different planes**, not a function of the single off-axis angle `arccos(x/r)`. Velocity is unaffected; the density angular profile is a different function from a single-angle source flow. |
| **SM-05** | The γ inside E5's integrand is hard-coded to 1.4, so the argument passed to `calculateNormCoeff` was ignored. Worse, A normalises E3 (θ-only) while the density applies E4 (θ·φ). A does not normalise the function it scales. |
| **SM-06** | No clipping at θ_ℓ. `f_θ` takes `abs()` of the cosine, so the profile folds and rises again past the limit. `f_φ` does not, so for \|φ\| > θ_ℓ the base is negative and the fractional exponent yields **NaN** — 0 of 2044 faces here, but 99 of 99 on the archived wake-cylinder meshes. |
| **SM-09** | `r = 0.5` is hard-coded in two unrelated places. Correct for this mesh, but it binds the code to one geometry, and applying a fixed radius to per-face centroids tilts U off radial by 3.1e−4. |
| **SM-10** | `r_e` — exit radius or throat radius? See above. |
| **GP-01** | The centroid divisor is `3.0` regardless of vertex count. **Fixed** in `plumetools.geometry.centroid`; bit-identical on triangles, 4/3 different on quads. |
| **AD-03** | `cylinder` and `plate` boundary entries are written into every field file; this mesh has neither. |

## Open questions

These need a human, or an environment this review did not have.

1. **What is the source of the source-flow model?** Specifically the
   `(γ + 0.41)/(γ − 1)` exponent — where does `0.41` come from?
2. **Is `r_e` the nozzle exit radius or the throat radius?** `f₂` implies the
   throat; the name says exit. The absolute density scales as `(r_e/r*)²`.
3. **Should the model use N₂ or argon?** The solver, output filename and
   archived validation all say argon.
4. **Was the separable θ·φ form intended**, or should the angular dependence be a
   function of the single off-axis angle from the plume axis?
5. **HA-06** — is `coordinateSystem dsmcAxisymmetric` valid on a 3D mesh, and did
   it affect the results already produced with it? Needs dsmcFoam+.
6. **HA-07** — does OpenFOAM v1706 load `type Unspecified;`? If not, the
   committed `3d-inflow` cannot be run as-is. Needs OpenFOAM.

## Model evolution

Reconstructed from the commit history. **Newer is not automatically more
correct**: two late changes were regressions.

| Quantity | Before | After | Commit | Verdict |
|---|---|---|---|---|
| θ_ℓ | `sqrt((g+1)/(g-1) - 1)` | `sqrt((g+1)/(g-1)) - 1` | `22a5fdb` | ✅ Fix. The old form gave 3.512 rad at γ=1.4 instead of 2.277 — every result before 2022-06-08 used a different limiting angle. |
| v_ℓ | `k = 1` | `k = 1.3806e-23` | `22a5fdb` | ✅ Fix. Boltzmann's constant was literally 1. |
| A | `sqrt((g-1) + (g+1))` | `sqrt((g-1) / (g+1))` | `22a5fdb` | ✅ Fix. `+` → `/`. |
| f(θ) | `cos(...)` (exponent 1) | `cos(...)**((g+0.41)/(g-1))` | `22a5fdb` | ⚠ Model change, unexplained. Introduced the uncited `0.41`. |
| φ | `atan(abs(y/x))` | `arctan2(y, x)` | `22a5fdb`+ | ⚠ Mixed. Correct quadrant, but extends φ to (−π, π] — **this created the NaN path**. |
| θ | `atan(abs(y/z))` | `acos(z/r)`, r = 0.5 | `ebc4588` | ⚠ Introduced the hard-coded radius; still measured from +z. |
| U | scalar `743` | full spherical → Cartesian | `4f1ac66` | ✅ Major improvement — a real 3D radial field. (`743` was the archived 2d-wedge inflow speed.) |
| T | `800` | `300` | `4f1ac66` | ⚠ Unexplained. `802.2867` was the 2d-wedge value. |
| Normals | `cross(A,B)` | `cross(B,A)` | `b22f573` | ➖ No effect — the value is never used. |
| Coord system | `dsmcAxisymmetric` | `dsmcCartesian` | `115d4f5` | ✅ Likely a fix — and evidence `3d-inflow` is misconfigured. |
| T₀ wiring | all literal | config in U and T only | `115d4f5` | ❌ **Partial fix = new bug.** Created SM-02. |
| `angularDependence` | θ·φ product | `return f_theta` | `7f6ec54` | ❌ Regression in the archived wake-cylinder cases. |
| — | — | added `input()` per face | `7f6ec54` | ❌ Regression. Debug code reached `main`. |

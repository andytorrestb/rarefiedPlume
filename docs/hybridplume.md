# `plumetools.hybridplume` — steady one-way CFD → DSMC free plume, full 3-D

A statistically steady, single-species argon plume from a circular exit into
vacuum, computed with a continuum precursor (OpenCFD `rhoCentralFoam`), a fixed
continuum-breakdown surface, and DSMC (standard `dsmcFoam` + this repository's
`plumeFieldInflow`) beyond it. Configuration attribution: Cai & Wang, *J.
Spacecraft and Rockets* 49(1), 2012, DOI 10.2514/1.A32046 — its physical source
is reproduced in 3-D; its axisymmetric numerics are **not**.

Everything is full three-dimensional Cartesian: no wedge, no symmetry plane, no
empty patch, no radial weighting. The configuration loader refuses those keys.

Case family: [`cases/cai2012-hybrid/`](../cases/cai2012-hybrid/README.md).

---

## 1. Scope of this increment

In: the 3-D precursor; breakdown diagnostics; automatic/manual/none interface
selection with an explicit unsuitable-interface verdict; the CFD → DSMC transfer
with measured flux accounting; DSMC moments, exact particle crossings, budgets;
an independent collisionless reference; resource estimates; resumable laptop and
HPC execution.

Out (by request): two-way coupling, moving interfaces, internal nozzle flow,
surfaces, chemistry, internal energy modes, turbulence models, axisymmetric
shortcuts.

Starting points recorded 2026-10-28:

| repository | ref | SHA |
|---|---|---|
| rarefiedPlume `main` (GitHub) | — | `a5eb48b` |
| rarefiedPlume local `refactor/plumetools-maintainability` (not on GitHub) | base of this work, by request | `8087bdb` |
| hyStrath `port/openfoam-v2512` | inspected, **not modified, not used** | `34008d0` |
| hyStrath `master` (checked out; uncommitted local edits) | untouched | `984e300` |

Work branch: `feature/hybridplume-3d` (local).

## 2. Backend decision

Standard `dsmcFoam` (v2512) with `plumeFieldInflow` is sufficient, so it is the
only backend. Read in the v2512 source:

* `VariableHardSphere::sigmaTcR` is Bird's VHS cross-section; NTC selection is
  the only one `DSMCCloud::collisions()` implements; `none` disables collisions
  for the collisionless verification (the type name in v2512 is `none` —
  `cai2012/dictionaries.py` still writes `NoBinaryCollision`, reported earlier);
* parcels are 3-D with 3-D velocities; one global weight;
* `plumeFieldInflow` injects across a named patch list with per-face `n`, `U`,
  `T`; any `patch`-type boundary deletes parcels (`DSMCParcel::hitPatch`
  returns false), which is exactly vacuum and interface absorption.

hyStrath `dsmcFoam+` was not needed: full 3-D removes the axis-weighting motive,
and nothing required (weighting, sampling, collisions) is missing from the
standard solver. Its inflow schema (vector `boundaryT` of Ttra/Trot/Tvib) is not
interchangeable with the scalar one used here.

## 3. The physical state (`physics.py`, the single source)

| quantity | value | label |
|---|---|---|
| `D` | 0.2 m | paper |
| `S0 = Ue / sqrt(2 R Te)` | 2 | paper |
| gas | argon | paper |
| `Kne = lambda_e / D` | 100, 0.1, 0.01 (paper); 1e-3, 1e-4 (extensions) | paper / extension |
| `Te` | 300 K | assumption (Cai validates normalised fields) |
| `m` | 6.63e-26 kg | repository |
| `d_ref`, `T_ref` | 4.17e-10 m, 273 K | repository = Bird (1994) Table A1 |
| `omega` | **0.81** | Bird Table A1 — the Cai family uses 0.74, the OpenFOAM N2 tutorial value |
| `gamma`, `Pr` | 5/3, 2/3 | monatomic |

`ae = sqrt(2RTe)` = 353.48 m/s is the thermal speed scale, not the sound speed;
`Ue` = 706.95 m/s; `Me = S0 sqrt(2/gamma)` = 2.191.

`ne = (Te/Tref)^(omega-1/2) / (sqrt(2) pi d^2 Kne D)` (the `vhs` diagnostic
convention, Bird eq. 4.76 with eq. 4.62's viscosity), then `rho_e = m ne`,
`pe = ne k Te`:

| Kne | ne [1/m³] | pe [Pa] | Re_D | tau_c [s] |
|---|---|---|---|---|
| 100 | 6.664e16 | 2.760e-4 | 0.027 | 5.01e-2 |
| 0.1 | 6.664e19 | 0.2760 | 27.4 | 5.01e-5 |
| 0.01 | 6.664e20 | 2.760 | 273.6 | 5.01e-6 |
| 1e-3 | 6.664e21 | 27.60 | 2736 | 5.01e-7 |
| 1e-4 | 6.664e22 | 276.0 | 27359 | 5.01e-8 |

**Transport consistent with the DSMC collision model**: `mu_ref =
15 sqrt(m k Tref/pi) / (2 (5-2w)(7-2w) d^2)` (Bird eq. 4.62, the first
Chapman–Enskog approximation for the VHS cross-section above) = 2.1154e-5 Pa s
— Bird's table gives 2.117e-5 (checked in the tests); `mu = mu_ref
(T/Tref)^omega`; `kappa = mu cp/Pr`. OpenFOAM v2512 cannot represent a power
law in `psiThermo` (Sutherland only), so
`applications/thermophysicalModels/hybridPlumeThermo` instantiates OpenFOAM's
own `tabulatedTransport` for `hePsiThermo` (≈20 lines, no new physics). The
table is 240 log-spaced points, 0.5–5000 K; its interpolation error is measured
at ≤ 2e-5 and recorded in every `setup.json`.

**Mean-free-path conventions** differ by constant factors and are reported side
by side: `hard_sphere/vhs` = 0.971 at 300 K; `viscosity_simple/vhs` =
`15 pi/(2(5-2w)(7-2w))` = 1.296. A threshold quoted under one convention must be
converted before it is compared with another.

## 4. Continuum precursor (`cfd.py`)

Laminar compressible NSF, Kurganov–Tadmor central-upwind fluxes. Source:
`fixedValue` p, T, U (Mach 2.19, every characteristic enters). `exitPlane`
(`x = 0`, `r > Re`) and `farField`: `waveTransmissive` p towards the background,
`inletOutlet` T and U — open outflow, consistent with the DSMC deleting parcels
on the same patches.

### Robustness: what was measured (2026-10-28, 4096-cell probes, 2 transits)

| floor `p_bg/pe` | `T_bg` | reconstruction | result |
|---|---|---|---|
| 1e-4 | 300 K | vanLeer, no viscous cap | T < 0 on step 1 (`e` solve) |
| 1e-4 | 300 K | Sutherland instead of table | same — not the tabulated library |
| 1e-2 | 300 K | vanLeer | T < 0 after ~1 transit |
| 1e-3 | 30 K | vanLeer / upwind / Minmod | fail in the startup shock |
| 1e-3 | 300 K | vanLeer | T < 0 at 0.4 transit |
| **1e-3** | **300 K** | **Minmod + viscous cap** | **2 transits, clean** |
| **1e-4** | **300 K** | **Minmod + viscous cap** | **2 transits, clean** |

Root cause of the first failure: `rhoCentralFoam` advances `div(tau.U)`
explicitly; its stability limit `dt < dx^2 rho/mu` is set by the background gas
(~1e-8 s at 1e-4, 300 K, 1 cm cells) and `maxCo` knows nothing about it.
`maxDeltaT` is therefore capped at `cfd.viscous_dt_factor dx^2 rho_bg/mu(T_bg)`
(0.25 by default; reported in `setup.json`). The cost is real: on D/20 cells at
Kne 0.01 the step is 4.8e-8 s, 33x the convective step.

### The background is a floor, and its effect is measured

The default is `p_bg/pe = 1e-3`. The CFD report shows it directly: through
`exitPlane` the solver's numerical mass flux `phi` is **inward** —
1.4% of the source mass flow at Kne 1e-3 and 4.4% at Kne 0.01 (D/10) — the floor
gas drawn in where the lip expansion falls below it. The reconstructed boundary
flux misses this (it uses the inletOutlet value), which is why both are written.
The breakdown step compares the retained core against any `cfd_bg*` copy run at
another floor (`setup --background-ratio 1e-4`, then `Allrun ... cfd --subdir
cfd_bg1e-04`).

## 5. Breakdown diagnostics (`breakdown.py`)

On the unaveraged 3-D CFD cells, weighted least-squares gradients over face
neighbours (OpenFOAM's `leastSquares`; exact for linear fields; one-sided at
boundaries):

* `Kn_rho = lambda|grad rho|/rho`, `Kn_T = lambda|grad T|/T`,
  `Kn_U = lambda|grad|U||/max(|U|, c)`, `Kn_GLL = max(...)`;
* `B_exp = -(lambda/c_bar) U.grad(rho)/rho`, signed (positive in expansion),
  a nominal-collision-time rate parameter — not interchangeable with other
  "Bird parameter" definitions.

Invalid cells (non-positive or non-finite state, unsolvable stencil) are NaN and
never continuum. The retained region is the 3-D face-connected component(s) of
continuum cells touching the source patch; disconnected pockets go to DSMC and
are counted. Cells within `lip_exclusion_over_D` of the lip circle are reported
separately: gradients there are mesh-dependent by construction. Thresholds
(`Kn_GLL`, `B_exp` = 0.05) are sensitivity settings; the report evaluates 0.5x,
1x and 2x.

## 6. Interface selection and suitability (`interface.py`)

The surface is a set of real internal mesh faces between retained and
non-retained cells, oriented CFD → DSMC. Per face: state, `Un`, `Mn = Un/c`,
indicators on both sides, NSF flux, Maxwellian forward and reverse half-range
fluxes. The region's whole boundary is checked to close (sum of outward area
vectors / boundary area: 5e-15 measured). Verdicts: `suitable`, `unsuitable`
(with reasons), `no_continuum_region`, `full_dsmc`. Criteria (all reported with
their values): returning mass fraction, subsonic-normal-Mach area fraction,
outer-boundary contact, DSMC-side overlap layers within 2x the thresholds,
invalid states; for a manual surface, retained non-continuum volume.

If the retained region does not cover the whole source disk, the uncovered
source faces stay in the DSMC mesh and inject the source state directly — the
source is never truncated.

A DSMC region on an `unsuitable` surface is built only with
`--allow-unsuitable`, and every report from it says `smoke_only: true`.

## 7. Transfer and flux accounting (`mapping.py`, `halfrange.py`)

Canonical interface fields (`interface_fields.npz`): `numberDensity`,
`massDensity`, `velocity`, `translationalTemperature`, `faceArea`,
`faceNormal`, `faceCentre`, `valid`, `cfdFace`, `cfdCell`. Face states are
linear (OpenFOAM weights) between the two CFD cells; `n = rho/m` is exact and
`p = n k T` is derived.

Two v2512 facts shaped the DSMC mesh step:

* `subsetMesh` creates a missing exposed patch as type **`empty`** — useless for
  3-D DSMC. `dsmc-setup` first appends a zero-face `interface` patch of type
  `patch` to `polyMesh/boundary`;
* `subsetMesh` writes **no** `faceMap`/`cellMap`. Faces and cells are matched
  geometrically (KD-tree) and verified: coincident centres, DSMC area vector =
  `-n A`, one-to-one. Unmatched or invalid faces are refused.

Four fluxes are kept distinct: the CFD net flux (NSF, with viscous and
conductive parts reported separately); the gross incoming half-range flux the
DSMC is asked to inject; the returning flux; and the actual DSMC net flux. The
half-range module is verified against quadrature (1e-14, including `s = -30`)
and forward − reverse reproduces the inviscid fluxes exactly; a Maxwellian
inlet carries **no** viscous stress or heat flux, and the report quantifies
that difference rather than claiming NSF flux matching.

## 8. DSMC measurements (`post.py`, `particles.py`)

* `fieldAverage` cell moments: `n`, `U = momentum/rhoM`, `T` from
  `linearKE` (never averaged itself);
* parcel-snapshot tensor moments: `C = P2/N - UU`, trace and directional
  temperatures, `T_parallel`/`T_perpendicular`; block standard errors;
* measured gross injection from the log's `Particles inserted` (exact count);
* **exact crossings**: `DSMCCloud::evolve` runs inflow → move → collide, so a
  parcel present at the end of step k and gone at k+1 left along
  `x_k + v_k t`; `Allrun ... crossings` writes the cloud every step for a window
  and every exit is ray-traced to its boundary (interface, source, exitPlane,
  farField) with mass `W m`, momentum `W m v_k`, energy `W m |v_k|^2/2`.

**Restart transient (measured, open issue).** Every dsmcFoam restart with
`plumeFieldInflow` inserts a burst for ~10 steps (smoke-kn1em3: 549, 983, 1075,
856, 735 … settling to ~590/step; identical inflow fields before and after;
with or without per-step writes). Its cause was not found. The crossing
analysis discards `dsmc.crossing_discard_steps` (12) post-restart steps and
reports the series. A production resubmission adds ~1300 parcels once per
restart — small against thousands of sampling steps, but it should be traced.

### Startup and stationarity

Cai's startup statement (at least 10^4 `t0`) uses the collision time of the
**Kne = 0.01** case (`t0 = lambda/c_bar` = 5.01e-6 s here) for every case; it is
not each case's own collision time, and at Kne = 100 it would be meaningless as
a collision count. This family therefore prescribes startup and sampling in
domain transits (`dsmc.startup_transits`, `dsmc.sample_transits`) and
**measures** stationarity: `dsmc-post` reports the relative change of the
molecule count over the last quarter of the run, and `cfd-report` the change of
domain mass and energy over the last transit. Cai's 10^4 `t0` (0.050 s) is
reported for comparison by the resource step, not imposed.

## 9. Resources (`resources.py`)

Before launching: cells, parcels (weight from `parcels_per_cell_target` in the
densest DSMC cell), memory, `dt` = min(Courant on `|U| + 3 sigma`, a fraction
of the smallest nominal collision time), steps, runtime, `dx/lambda`, step over
collision time, occupancy, densest cell. Measured on the laptop (WSL, `/mnt/c`):
4.0e-6 s per parcel-step per rank at 39k parcels (smoke-kn1em3, 4 ranks),
1.0–1.4e-5 at 20–30k parcels on 3 ranks — the default is 4e-6 and should be
re-measured on the cluster. Rank imbalance is measured from processor clouds
(1.02–1.39 in the smoke runs).

A dense case is never called affordable because its mesh fits in memory; a
smaller weight reduces noise but does not repair `dx/lambda > 1`.

## 10. Verification: hierarchy and what was executed

| level | status |
|---|---|
| analytical: half-range moments, NSF flux signs, collisionless reference vs Cai closed form | **executed** (unit tests; 1e-13 agreement) |
| collisionless full-3-D DSMC vs reference (`smoke-collisionless`) | **executed**, see `results/` of that case |
| Cai Kne 0.1, 0.01 full DSMC at verification resolution | prepared (`verification-*` profiles, HPC) |
| hybrid vs full DSMC on one mesh (`smoke-kn0p01-coarse[-full]`) | **executed**, smoke resolution, interface unsuitable |
| Kne 1e-3, 1e-4 hybrid extensions | smoke executed at 1e-3 (unsuitable, domain too short); production prepared |

Executed smoke results (2026-10-28, all **intentionally under-resolved**):

| quantity | smoke-kn1em3 (D/20) | smoke-kn0p01-coarse (D/10) |
|---|---|---|
| CFD mass residual / source flow | 1.8e-6 | −4.6e-5 |
| exitPlane numerical inflow / source | 1.4% | 4.4% |
| retained cells / covered source area | 38 308 / 100% | 554 / 75% |
| interface verdict | unsuitable: return 0.68, subsonic Mn 36%, far-field contact | unsuitable: return 0.25, subsonic Mn 87% (total Mach 2.8) |
| measured / expected gross injection − 1 | −3.7e-5 (±8.4e-4) | +1.2e-6 (±2.7e-3) |
| exact interface return (40 steps, 14 578 parcels) | 4.130e-3 kg/s vs CFD-Maxwellian 4.598e-3 (−10%) | — |
| injected (log, post-discard) vs expected | 6.704e-3 vs 6.717e-3 kg/s | — |
| DSMC net through interface vs CFD continuum net | 2.57e-3 vs 2.12e-3 kg/s (+21%) | — |
| snapshot return estimator vs exact | 3.93e-3 vs 4.13e-3 (−5%) | — |
| retained CFD vs full DSMC (n, mean abs) | — | 2.9% (noise floor ≥ 2.2%) |
| hybrid DSMC vs full DSMC (n, mean abs) | — | 8.1% |
| net transverse / axial momentum | 1.5e-3 | 6.5e-4 |

Collisionless full-3-D DSMC (`smoke-collisionless`, D/20 cells, 64 000 cells,
~0.7 M parcels, collisions `none`, 1664 averaged steps) against the
independent reference averaged over the same cells:

| quantity | result |
|---|---|
| centreline n: mean bias / rms / max | +0.08% / 1.3% / 4.6% (max in the exit-plane cell, x/D = 0.025) |
| centreline T rms, Ux rms | 0.45%, 0.30% |
| off-axis rings r/D = 0.25 / 0.5 / 0.75, n mean bias | 0.07% / 0.28% / 2.2% |
| measured vs expected gross injection | +4.4e-7 (counting error 7.0e-4) |
| net transverse / axial momentum | 3.4e-4, 2.7e-4 |
| azimuthal spread of n (median / max over (x, r) bins) | 7.1% / 36% (sampling noise at ~4 parcels/cell/step) |
| T_parallel / T_perp on the axis | 0.94 at the exit rising to 3.75 at x/D = 1.8 |
| cost | 4.1e-6 s per parcel-step per rank (4 ranks), rank imbalance 1.22 |

Mesh sensitivity of the breakdown surface, Kne 0.01 precursor (same domain,
floor and thresholds; D/10 = 8 000 cells, D/20 = 64 000 cells):

| quantity | D/10 | D/20 |
|---|---|---|
| CFD mass residual / source flow | −4.6e-5 | 8.0e-6 |
| exitPlane numerical inflow / source | 4.4% | 2.8% |
| retained volume [m³] (0.5x / 1x / 2x thresholds) | 0.0023 / 0.0044 / 0.0116 | 0.0023 / 0.0041 / 0.0099 |
| retained extent on the axis | 0.25 m | 0.245 m |
| source area covered by the retained region | 75% | 82% |
| interface returning mass fraction | 0.250 | 0.252 |
| subsonic-normal-Mach area fraction (total Mach, area mean) | 0.87 (2.82) | 0.85 (2.78) |
| max Kn_GLL within 2 cm of the lip | 1.0 | 2.5 — mesh-dependent, as expected at the singular lip |
| verdict | unsuitable | unsuitable (and 0 overlap layers) |

The verdict does not depend on the mesh; the threshold factor moves the
retained volume by 4x, the mesh by 7%.

None of these is a converged or validated result. No numerical tolerance was
imposed; the numbers are reported, with their noise floors where measurable.

## 11. Reproduction

Laptop (WSL, 32 GB host, 15 GB VM; never all cores):

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc; export PATH="$FOAM_APPBIN:$PATH"
applications/thermophysicalModels/Allwmake   # libhybridPlumeThermo.so, once per environment
applications/dsmcBoundaryModels/Allwmake     # libplumeDsmcBoundaryModels.so
C=cases/cai2012-hybrid/runs/kn0p01c
cases/cai2012-hybrid/Allrun $C setup --profile smoke-kn0p01-coarse -np 4
cases/cai2012-hybrid/Allrun $C cfd -np 4
cases/cai2012-hybrid/Allrun $C post-cfd [--allow-unsuitable]
cases/cai2012-hybrid/Allrun $C dsmc -np 4
cases/cai2012-hybrid/Allrun $C post-dsmc
cases/cai2012-hybrid/Allrun $C crossings -np 4      # optional, exact exits
python -m plumetools.hybridplume overlap $C --reference cases/cai2012-hybrid/runs/kn0p01c_full
```

HPC (Rocky Linux): compile both libraries there, then
`python -m plumetools.hybridplume hpc <case> --scheduler pbs --np 192` writes
`<case>/hpc/run_<stage>.pbs` sourcing `cases/cai2012-kn-extension/hpc_env.sh`
(override with `HYBRID_HPC_ENV`); `qsub` them in stage order. A walltime kill is
resumed by resubmitting the same script. SLURM scripts carry `__SET_ME__`
partition/account placeholders: this repository has no SLURM settings.

## 12. Limitations

* Smoke meshes are 5–20 mean free paths per cell at the source; nothing here
  is converged. The verification and production profiles are prepared, not run.
* The CFD needs a background floor; its inward leakage at `x = 0` is 1–4% of
  the source flow in the smoke runs. Floor sensitivity must be run per case.
* The Maxwellian interface inlet cannot carry viscous stress or heat flux.
* The returning flux is absorbed at the interface (one-way coupling). Where the
  returning fraction is large the hybrid is unsuitable, and the tool says so.
* Exact crossings need per-step cloud writes; the snapshot estimator is cheaper
  and was 9% low against the exact count in the one comparison made.
* The dsmcFoam restart insertion burst is unexplained.
* The staircase source area error (+0.59% at D/20, +1.86% at D/10) is a direct
  multiplier on injected and continuum fluxes; it is reported, never rescaled.

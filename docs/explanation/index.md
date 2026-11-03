# Explanation

Why the models are what they are, where their equations live in code, what they
assume, and how far the evidence for each goes. Each topic has one canonical
page; tutorials and API pages link here rather than repeating derivations.

## Models, their code and their evidence

| Model | Equations in code | Key assumptions | Tests | Canonical page |
|---|---|---|---|---|
| Legacy plume source flow (3d-inflow, ISS, Fnum) | `plumetools.sourceflow` — E1–E8 | reconstructed from code; **no recoverable citation**; N₂ model gas for an argon solver, θ from +z, ignored pressure — all frozen (SM-01 … SM-10) | `test_sourceflow.py`, regression golden | [The plume source-flow model](../source-flow-model.md) |
| Paper-faithful source flow (Markelov) | `plumetools.markelov1999.sourceflow` — P1 … | as printed in AIAA 99-3455; the 0.41 exponent offset's provenance; inflow hemisphere radius and stagnation temperature are assumptions | `test_markelov_sourceflow.py`, `test_markelov_flux.py` | [The AIAA 99-3455 family](../markelov1999-case.md) |
| Cai's collisionless solution | `plumetools.cai2012.analytical` — C0–C5 | re-derived from the free-molecular integral, not transcribed; uniform drifting Maxwellian on the disk | `test_cai2012_analytical.py` | [The Cai 2012 family §6](../cai2012-case.md#6-the-analytical-comparison-solution) |
| Cai derivation chain (Kn → λ₀ → n₀, S₀ → U₀, weight, Δt) | `plumetools.cai2012.gas` (G1, G2), `cai2012.inflow`, `cai2012.mesh` | T₀ = 300 K, argon VHS constants and the `vhs` λ convention are assumptions; Δt from a Courant target, not Cai's Δt/t₀ = 1 | `test_cai2012_gas.py`, `test_cai2012_inflow.py`, `test_cai2012_mesh.py` | [The Cai 2012 family §2–5](../cai2012-case.md#2-the-derivation-chain) |
| Nozzle-exit state (any gas, optional plate) | `plumetools.nozzleplume.exitstate`, `nozzleplume.gas` | uniform exit profile; one species; diffuse 300 K wall; hard sphere for R-4D | `test_nozzleplume_*`, OpenFOAM conservation test | [r-4d, explained](../r-4d-case.md) |
| Rarefaction diagnostics | `plumetools.rarefaction.core` | averaged moments; λ/Δ ≥ 3; > 20 sampled parcels per step | `test_rarefaction.py` | [Rarefaction diagnostics](../rarefaction-diagnostics.md) |
| Kn-extension feasibility | `plumetools.cai2012ext.feasibility` | the production criteria are the Cai family's, unrelaxed | `test_cai2012ext_policy.py` | [The Kn extension §5](../cai2012-kn-extension.md#5-resolution-and-feasibility-before-meshing) |
| Interface closure (Cai exit → angular source → sphere) | `plumetools.cai2012ext.interface` | Markelov's sonic-orifice angular law imposed on a supersonic exit (measured to under-predict axis density); Bird's B threshold | `test_cai2012ext_interface*.py`, `test_cai2012ext_overlap.py` | [The interface track](../cai2012-kn-interface.md) |
| Empirical field scaling | `plumetools.cai2012ext.scaling`, `scaling_fit` | coordinate-only maps; no conservation; no physical model | `test_cai2012ext_scaling*.py` (synthetic only) | [Empirical whole-field scaling](../cai2012-scaling.md) |
| Reduced virtual collisions | `plumetools.virtualcollision.model`, `kernel`, `source` | single species, elastic VHS; mutually thinned pair events in a finite plume | `test_virtualcollision_*.py` | [Reduced virtual collisions](../cai2012-virtual-collision.md) |
| Hybrid CFD → DSMC | `plumetools.hybridplume.physics`, `halfrange`, `fluxes`, `breakdown`, `interface`, `mapping` | one-way coupling; a Maxwellian inlet (no viscous stress or heat flux); a numerical background floor in the CFD | `test_hybridplume_*.py` | [`plumetools.hybridplume`](../hybridplume.md) |
| Per-face patch-selected DSMC inflow | `applications/dsmcBoundaryModels/plumeFieldInflow` (C++) | `FreeStream`'s Bird eqn 4.22 flux and eqn 12.5 acceptance–rejection, unchanged; patch list and per-face density changed | `tests/openfoam` | [plumeFieldInflow](../plume-field-inflow.md) |

## Cross-cutting topics

- [Rarefaction, collisions and the criteria built on them](rarefaction-and-collisions.md)
  — mean-free-path and Knudsen conventions side by side, resolution and
  time-step criteria, collisionless vs collisional, breakdown parameters,
  conservation checks.
- [What the tests establish](evidence.md) — tiers, the last measured run, and
  what no test establishes.
- [Standard `dsmcFoam` vs `dsmcFoam+`](../solver-compatibility.md) — why the
  repository needs its own inflow model, and the physics compromises of stock
  `FreeStream`.
- [Meshes](../mesh.md) — the 3d-inflow mesh pipelines and why generation waited
  for the centroid fix.
- [Slice imagery](../viz-slices.md) — how the pictures are drawn and the ways
  they mislead.
- [Environment and provenance](../environment.md) — what produced the
  historical results and what cannot be recovered.

## Family explanations

[Cai 2012](../cai2012-case.md) ·
[cai2012-health](../cai2012-health.md) ·
[Kn extension](../cai2012-kn-extension.md) ·
[Interface track](../cai2012-kn-interface.md) ·
[Field scaling](../cai2012-scaling.md) ·
[Virtual collisions](../cai2012-virtual-collision.md) ·
[Hybrid](../hybridplume.md) ·
[Markelov 1999](../markelov1999-case.md) ·
[r-4d](../r-4d-case.md) ·
[Source-flow model](../source-flow-model.md)

Literature: the citations the repository carries are the ones in these pages
and in the module docstrings (Cai & Wang 2012; Lumpkin, Stewart & Markelov 1999;
Stewart & Lumpkin 2012 for the archived validation case; Bird for the DSMC
relations; Boyd, Chen & Candler 1995 and Wang & Boyd 2003 for the gradient-length
Knudsen number). Where a model has no citation, its page says so; none was added
here.

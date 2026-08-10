# Standard `dsmcFoam` vs the MNF `dsmcFoam+` fork

Everything below was verified against **OpenFOAM v2512** by reading its source and
running the case, not inferred.

## The short version

**Standard `dsmcFoam` cannot express this case's boundary conditions.** It will
run, but only after two changes that alter the physics:

| | Standard `dsmcFoam` | MNF `dsmcFoam+` |
|---|---|---|
| Per-face inflow **velocity** | ✅ `0/boundaryU`, read per face | ✅ |
| Per-face inflow **temperature** | ✅ `0/boundaryT` — but **volScalarField** | ✅ — **volVectorField** `(T 0 0)` |
| Per-face inflow **number density** | ❌ one scalar per species | ✅ `dsmcFreeStreamInflowFieldPatch` |
| Choose *which* patches admit inflow | ❌ every `patch`-type boundary | ✅ named per patch |
| Particle-deleting outflow | ✅ any `patch` boundary already does — see [below](#blocker-3-was-not-a-blocker) | ✅ `dsmcDeletionPatch` |

The first two rows are what
[`plumeFieldInflow`](plume-field-inflow.md) adds to standard `dsmcFoam`. This
document describes `cases/3d-inflow`, which does **not** use it — that case is
bound to the regression golden and keeps the workaround below.

The plume's angular *direction* survives. Its angular *density profile* — the
thing the source-flow model exists to compute — does not.

## The three blockers, with sources

### 1. Number density is uniform, not per-face

`FreeStream.C:93`:

```cpp
numberDensities_[i] = numberDensitiesDict.get<scalar>(molecules[i]);
```

One `scalar` per species, read from `constant/dsmcProperties`. There is no field
form. On `cases/3d-inflow` the model's per-face density spans
**2.669e19 → 5.570e21**, a factor of **209**; collapsing that to a constant
discards the entire angular structure.

`runInflow.py` prints the area-weighted mean (**1.915737e+21** for the reference
case), which preserves total particle flux and is the least-wrong single value.
It is a collapse, not a translation.

### 2. Inflow cannot be restricted to the inflow patch

`FreeStream.C:57-62`:

```cpp
forAll(cloud.mesh().boundaryMesh(), p)
{
    const polyPatch& patch = cloud.mesh().boundaryMesh()[p];
    if (isType<polyPatch>(patch))
    {
        patches.append(p);
    }
}
```

Every boundary of geometric type `patch` gets free-stream injection. There is no
selection list. This mesh has `inflow` *and* `vacuum` as `patch`, so the vacuum
boundary would inject plume gas inward.

`isType<>` is an exact type match, so `wall`, `symmetry`, `empty` and `wedge`
patches are excluded — which is the only lever available (see the workaround
below).

Symptom if you do nothing: `boundaryT` is `zeroGradient` on `vacuum`, so it
evaluates to the internal field, 0, and the run aborts at the first timestep with

```
--> FOAM FATAL ERROR: Zero boundary temperature detected, check boundaryT condition.
    From void Foam::FreeStream<CloudType>::inflow()  ... FreeStream.C at line 186
```

### 3. There is no particle-deleting boundary

The complete set of `WallInteractionModel`s in v2512:

```
MaxwellianThermal   MixedDiffuseSpecular   SpecularReflection
```

All three reflect. None deletes. A plume expanding into vacuum needs particles to
leave the domain, which is what the fork's `dsmcDeletionPatch` does — and what
`system/boundariesDict` in this case already configures.

## Making it run anyway

Set `mesh.outer_patch_type: wall` in `case.yaml`. The generator then emits the
outer boundary as a `wall`, which `FreeStream` excludes, and the case runs to
completion — verified end to end in v2512:

```
Particles inserted        = 59251
Collisions                = 518485
Number of dsmc particles  = 296218
ExecutionTime = 2.41 s
End
```

**But the physics is then wrong in two ways**, and neither is subtle:

1. `vacuum` reflects instead of absorbing. The plume is expanding into a closed
   box, not into vacuum. Density will build up rather than reaching steady state.
2. The inflow density is uniform, so the plume has no angular structure.

That may still be useful for a smoke test, a mesh check, or a performance
measurement. It is not a plume-impingement calculation.

It is opt-in rather than the default because reflecting the outflow is a physics
decision. Leaving `outer_patch_type` at `patch` under `dialect: standard` makes
`load_case_config` warn, so the problem is reported when the config is read
rather than as an MPI stack trace at the first timestep:

```
ConfigWarning: output.dialect is 'standard' and mesh.outer_patch_type is 'patch'.
Standard dsmcFoam's FreeStream injects on EVERY patch-type boundary
(FreeStream.C:57-62), so the outer boundary becomes a second inflow and the run
aborts at the first timestep with 'Zero boundary temperature detected'.
```

`Allrun` also extracts the first `FOAM FATAL` block from the log rather than
tailing it — under `mpirun` the last lines are MPI teardown from every rank and
the real error has scrolled past.

## What this repository does about it

`case.yaml` carries `output.dialect`:

```yaml
output:
  dialect: standard    # or: mnf
```

* **`standard`** — writes `0/boundaryT` as a **volScalarField**, and generates the
  nine zeroed measurement fields (`dsmcRhoN fD iDof internalE linearKE momentum q
  rhoM rhoN`) that `dsmcFoam` requires but `dsmcInitialise` does not create.
  `runInflow.py` prints the area-weighted number density for `FreeStreamCoeffs`.
* **`mnf`** — writes `0/boundaryT` as a **volVectorField** holding `(T 0 0)`, the
  pre-refactor behaviour the regression golden pins.

Writing the wrong type is a hard read failure, not a silent mis-run — the two
solvers genuinely disagree about what `boundaryT` is.

`0/boundaryNumberDensity_<species>` is written under both dialects. Standard
`dsmcFoam` ignores it; it is kept because it is the only record of what the model
actually computed, and because a custom `InflowBoundaryModel` would read exactly
that file.

## Dictionary mapping

`constant/dsmcProperties.mnf` and `system/dsmcInitialiseDict.mnf` are the
originals, kept for provenance.

| MNF | Standard | Note |
|---|---|---|
| *(absent — set per patch in `boundariesDict`)* | `WallInteractionModel` | required; `SpecularReflection` chosen as neutral |
| *(absent)* | `InflowBoundaryModel` + `FreeStreamCoeffs` | required |
| `rotationalRelaxationCollisionNumber` | `relaxationCollisionNumber` + `Tref` | renamed, and `Tref` added |
| `rotationalDegreesOfFreedom` | `internalDegreesOfFreedom` | renamed |
| `alpha` (VSS scattering) | *(none)* | dropped |
| `collisionPartnerSelectionModel` | *(none)* | MNF-only |
| `coordinateSystem dsmcAxisymmetric` | *(none)* | MNF-only — and it was wrong here anyway (HA-06) |
| `configurations ( { type dsmcMeshFill; ... } )` | flat `numberDensities`/`temperature`/`velocity` | `dsmcInitialiseDict` schema differs entirely |
| `boundariesDict`, `fieldPropertiesDict`, `controllersDict`, `chemReactDict`, `loadBalanceDict` | *(none)* | MNF-only; standard `dsmcFoam` ignores them |

One incidental benefit: dropping `coordinateSystem dsmcAxisymmetric` resolves
**HA-06** by construction. Standard `dsmcFoam` is Cartesian, so the axisymmetric
setting this case inherited from the archived 2d-wedge case — on a 3D mesh with a
φ-dependent inflow — can no longer apply.

## If you need the real boundary conditions

Option 2 below **has since been written**, and it resolves blockers 1 and 2. It
also establishes that blocker 3 was not a blocker at all.

1. **Use the MNF fork.** `boundariesDict` in this case already configures
   `dsmcFreeStreamInflowFieldPatch` on `inflow` and `dsmcDeletionPatch` on both
   `inflow` and `vacuum`. Nothing needs writing — set `output.dialect: mnf` and
   `application dsmcFoam+` in `controlDict`.
2. **Use `plumeFieldInflow`** —
   [`applications/dsmcBoundaryModels/plumeFieldInflow`](../applications/dsmcBoundaryModels/plumeFieldInflow),
   documented in [`plume-field-inflow.md`](plume-field-inflow.md). It subclasses
   `InflowBoundaryModel<CloudType>`, takes a patch-name list from the dictionary
   instead of `isType<polyPatch>`, and reads the per-face number density from a
   `volScalarField` where `FreeStream` uses `numberDensities_[i]`. Everything else
   is `FreeStream`'s injection algorithm, kept line-for-line.

   `cases/markelov1999` uses it, and therefore keeps its vacuum boundary as an
   open `patch` rather than a reflecting `wall`.

   This case (`3d-inflow`) is **not** switched to it. Its configuration is bound
   to the regression golden, and changing the boundary condition would change
   the results the golden exists to pin. The workaround above remains what this
   case does.
3. **Accept the compromises** for smoke tests and mesh work only, per above.

### Blocker 3 was not a blocker

The "there is no particle-deleting boundary" entry above is **wrong**, and reading
`DSMCParcel.C` and `particleTemplates.C` in v2512 shows why.

`DSMCParcel::hitPatch` returns `false`, so `particle::hitBoundaryFace` falls
through its dispatch chain — wedge, symmetryPlane, symmetry, cyclic, cyclicACMI,
cyclicAMI, processor, wall — and off the end:

```cpp
else
{
    td.keepParticle = false;
}
```

An ordinary `patch` boundary therefore **already deletes** outgoing particles.
The three `WallInteractionModel`s all reflect, but they only ever run on `wall`
patches.

What made the open boundary unusable was `FreeStream` *injecting* on it, not any
failure to delete. Remove that — which a patch-name list does — and no custom
outflow or deleting model is needed. `cases/markelov1999` runs with `vacuum` and
`upstreamVacuum` as plain `patch` boundaries and no deleting model at all;
`tests/openfoam` asserts by measurement that particles leave.

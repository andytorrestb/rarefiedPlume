# `plumeFieldInflow` — a per-face, patch-selected DSMC inflow model

`applications/dsmcBoundaryModels/plumeFieldInflow/`

An `InflowBoundaryModel` for standard OpenFOAM `dsmcFoam` that injects particles
across an **explicitly named list of patches**, using a **per-face number
density** read from a `volScalarField`.

**OpenFOAM version:** written and verified against **v2512**
(`linux64GccDPInt32Opt`). Every source reference below was read in that tree. See
[Version dependence](#version-dependence).

---

## Why stock `FreeStream` is not enough

Two independent reasons, both read out of `FreeStream.C`.

### 1. It injects on every `patch`-type boundary

`FreeStream.C:54-64`:

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

There is no selection list. A plume domain's outer boundary must be an open
`patch` so particles can leave — and that turns it into a second inflow. The
symptom is immediate:

```
--> FOAM FATAL ERROR: Zero boundary temperature detected, check boundaryT condition.
    From void Foam::FreeStream<CloudType>::inflow()  ... FreeStream.C at line 186
```

because `boundaryT` is `zeroGradient` there and evaluates to the zero internal
field.

The usual workaround — declaring the outer boundary a `wall`, which the exact
`isType<>` match excludes — makes it **reflect instead of absorb**. The plume then
expands into a closed box. That is a different physical problem, and it is what
`cases/3d-inflow` had to accept (see
[`solver-compatibility.md`](solver-compatibility.md)).

### 2. The number density is one scalar per species

`FreeStream.C:93`:

```cpp
numberDensities_[i] = numberDensitiesDict.get<scalar>(molecules[i]);
```

There is no field form. A source-flow inflow surface carries a strong angular
density profile, and collapsing it to its area-weighted mean discards exactly the
structure the analytical model exists to produce.

Per-face **velocity** and **temperature** already work: `FreeStream` reads
`cloud.boundaryU()` and `cloud.boundaryT()` face by face. So the plume's angular
*direction* survives stock OpenFOAM; its angular *density* does not.

---

## What this model changes

Exactly two things. Everything else — the Bird eqn 4.22 flux accumulator, the
eqn 12.5 acceptance-rejection on the normal velocity, the triangle-area-weighted
injection position, the equipartition internal energy — is `FreeStream`'s,
deliberately kept line-for-line so a future upstream change can be diffed in.
The two changes are marked `CHANGED` in the source.

```
constant/dsmcProperties

    InflowBoundaryModel   plumeFieldInflow;

    plumeFieldInflowCoeffs
    {
        patches ( inflow );                  // CHANGED 1: explicit, not isType<>

        numberDensityFields
        {
            N2    boundaryNumberDensity_N2;  // CHANGED 2: a field, not a scalar
        }
    }

system/controlDict

    libs ( "libplumeDsmcBoundaryModels.so" );
```

The named field is an ordinary `volScalarField` in the time directory, whose
**boundary values** on the inflow patch carry the number density in 1/m³.
`plumetools.markelov1999.foamfields` writes it.

---

## What is deliberately absent: an outflow model

**None is needed.** Standard OpenFOAM already deletes particles crossing a plain
`patch`.

`DSMCParcel.C`:

```cpp
template<class ParcelType>
template<class TrackCloudType>
bool Foam::DSMCParcel<ParcelType>::hitPatch(TrackCloudType&, trackingData&)
{
    return false;
}
```

Returning `false` sends `particle::hitBoundaryFace` (`particleTemplates.C`) down
its dispatch chain — wedge, symmetryPlane, symmetry, cyclic, cyclicACMI,
cyclicAMI, processor, wall — and off the end:

```cpp
else
{
    td.keepParticle = false;
}
```

So an ordinary `patch` is already a correct absorbing outflow. Reflection happens
only on `wall` patches, through `DSMCParcel::hitWallPatch`, which is also what
records the surface fluxes (`fD`, `q`, `rhoN`, …) and applies the
`WallInteractionModel`.

The only thing that made this unusable before was `FreeStream` *injecting* on
those same open boundaries. Remove that — which the patch list does — and no
custom deleting model is needed. Adding one would duplicate behaviour OpenFOAM
already has.

`tests/openfoam/test_markelov_openfoam.py::test_particles_leave_through_the_open_boundaries`
asserts this by measurement rather than by reading the source: with a reflecting
boundary the particle count would equal everything ever supplied, and it is
lower.

**Consequence for patch types.** They are load-bearing, not labels:

| patch | type | behaviour |
|---|---|---|
| `inflow` | `patch` | injected across (named in the coeffs); particles reaching it are deleted |
| `cylinder`, `plate` | `wall` | `hitWallPatch` runs — records `fD`/`q`, applies the reflection model |
| `vacuum`, `upstreamVacuum` | `patch` | particles deleted |
| `symmetry` | `symmetry` | reflected |

A body declared `patch` would absorb the whole plume and exert no force. The
model refuses a `wall` named as an inflow for the mirror-image reason.

---

## Failure behaviour

Every failure is **fatal and named**. There is no fallback to a uniform density
and no silent skipping — a run that quietly degraded to `FreeStream` behaviour
would produce plausible-looking output for the wrong boundary condition, which is
the failure mode the class exists to remove.

| Condition | Result |
|---|---|
| `patches` empty | `FatalError` — an empty list is always a configuration mistake |
| A named patch is not in the mesh | `FatalError`, listing the mesh's patches |
| A named patch is a `wall` | `FatalError` — a wall is an interaction site, not a source |
| `numberDensityFields` empty | `FatalError` |
| A species is not in `typeIdList` | `FatalError`, listing the cloud's type ids |
| The named field is missing | `FatalError` from `MUST_READ` |
| A negative number density | `FatalError`, with the range |
| Zero boundary temperature on an inflow patch | `FatalError` |
| The library is not loaded | OpenFOAM reports `plumeFieldInflow` as an unknown `InflowBoundaryModel` and lists what it knows |

`cases/markelov1999/baseCase/Allrun` checks for the library *before* running
`dsmcInitialise` and `decomposePar`, so a missing build costs seconds rather than
minutes.

---

## Implementation notes

### Number densities are copied, not scaled in place

`FreeStream` divides its scalars by `nParticle` once at construction, because the
flux accumulator counts **parcels**, not molecules:

```cpp
numberDensities_ /= cloud.nParticle();
```

The obvious translation — dividing the `volScalarField` — is **wrong**, and was a
real bug during development. It scales every patch field in the mesh, including
`zeroGradient` ones whose storage OpenFOAM has not initialised at that point.
Arithmetic on that raised a floating-point exception on the first timestep, deep
inside `inflow()`:

```
#3  Foam::plumeFieldInflow<...>::inflow() in libplumeDsmcBoundaryModels.so
```

So each injection patch's values are **copied out and scaled**, once, at
construction. The field on disk stays in physical 1/m³, which is also what anyone
inspecting `boundaryNumberDensity_N2` expects.

### Restarts

The fields are read `MUST_READ` with `AUTO_WRITE`, exactly as `DSMCCloud`
constructs `boundaryT` and `boundaryU`. `AUTO_WRITE` is what makes a restart
work: the field is rewritten into each output time directory, so a run resumed
from `t > 0` finds it there rather than only in `0/`.

### Parallel

Accumulators are sized per **local** patch face count. In a decomposed run every
rank holds every patch, most with zero faces, so this is correct in serial and in
parallel with no special case: a rank owning no inflow faces accumulates nothing
and its empty patch is skipped. `decomposePar` decomposes the number-density
field along with every other field in the time directory.

### Registration

`makePlumeDsmcBoundaryModels.C` calls `makeInflowBoundaryModelType` only. The
selection **table** is already defined by OpenFOAM's own
`makeDSMCParcelInflowBoundaryModels.C`; calling `makeInflowBoundaryModel` again
would define it twice and fail to link.

---

## Building

```bash
. /path/to/OpenFOAM/etc/bashrc
cd applications/dsmcBoundaryModels
./Allwmake                # -> $FOAM_USER_LIBBIN/libplumeDsmcBoundaryModels.so
./Allclean                # remove objects and the installed library
./Allclean --dry-run      # list what would be removed
```

`Make/options` links `-lDSMC` — OpenFOAM's capitalisation, from
`src/lagrangian/DSMC/Make/files`. `-ldsmc` fails with `cannot find -ldsmc`.

Builds clean with no warnings under `-Wall -Wextra -Wold-style-cast
-Wnon-virtual-dtor`.

---

## Version dependence

Verified against **OpenFOAM v2512**. The model depends on these interfaces:

| Interface | Where |
|---|---|
| `InflowBoundaryModel<CloudType>` and its selection table | `src/lagrangian/DSMC/submodels/InflowBoundaryModel/` |
| `cloud.boundaryT()`, `cloud.boundaryU()`, `cloud.nParticle()`, `cloud.constProps()`, `cloud.typeIdList()`, `cloud.addNewParcel()`, `cloud.maxwellianMostProbableSpeed()`, `cloud.equipartitionInternalEnergy()` | `DSMCCloud` |
| `polyMeshTetDecomposition::faceTetIndices` | `meshTools` |

These have been stable across the v1706–v2512 range the repository spans, but
none of that is verified here. If the build fails on another version, compare
`plumeFieldInflow.C` against that version's `FreeStream.C` — they are kept
line-for-line for exactly that purpose.

---

## Related documents

* [`markelov1999-case.md`](markelov1999-case.md) — the case family that uses this
* [`solver-compatibility.md`](solver-compatibility.md) — the fuller comparison of
  standard `dsmcFoam` against the MNF `dsmcFoam+` fork
* [`cases/markelov1999/README.md`](../cases/markelov1999/README.md) — running it

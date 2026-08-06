/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | rarefiedPlume: plume-impingement DSMC cases
   \\    /   O peration     |
    \\  /    A nd           | https://github.com/andytorrestb/rarefiedPlume
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later, matching OpenFOAM, from which FreeStream is derived.

    This file is part of a derived work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

Class
    Foam::plumeFieldInflow

Group
    grpDsmcInflowBoundaryModels

Description
    An InflowBoundaryModel that injects particles across an explicitly named
    list of patches, using a PER-FACE number density read from a volScalarField.

    Written for a rarefied plume expanding from an analytical source-flow
    solution, where the inflow surface carries a strong angular density profile.
    Stock FreeStream cannot express that case for two independent reasons:

      1. It injects across EVERY boundary of geometric type `patch`, chosen by
         `isType<polyPatch>(patch)` with no selection list
         (FreeStream.C:54-64, OpenFOAM v2512). On a plume domain, whose outer
         boundary must be an open `patch` so particles can leave, that turns the
         vacuum boundary into a second inflow. The usual workaround -- declaring
         the outer boundary a `wall` so the exact type match excludes it -- makes
         it reflect instead of absorb, which is a different physical problem.

      2. It reads ONE scalar per species from constant/dsmcProperties
         (FreeStream.C:93, `numberDensities_[i] = numberDensitiesDict.get<scalar>`).
         There is no field form. Collapsing a source-flow profile to its
         area-weighted mean discards the entire angular structure the model
         exists to compute -- a factor of ~200 across the surface on the
         reference case.

    Per-face velocity and temperature already work in stock FreeStream, which
    reads cloud.boundaryU() and cloud.boundaryT() face by face. Only the number
    density and the patch selection are missing, so this class changes exactly
    those two things and inherits the rest of the injection algorithm unaltered:
    the Bird eqn 4.22 flux accumulator, the eqn 12.5 normal-velocity sampling,
    the triangle-weighted face position, and the equipartition internal energy.

    What is NOT here, deliberately
    ------------------------------
    No outflow or particle-deleting model. Standard OpenFOAM already deletes
    particles that cross a plain `patch`: DSMCParcel::hitPatch returns false, so
    particle::hitBoundaryFace falls through its dispatch chain to
    `td.keepParticle = false` (particleTemplates.C). Reflection happens only on
    `wall` patches, via DSMCParcel::hitWallPatch. Once unwanted injection is
    removed -- which is what the patch list above does -- an ordinary `patch`
    boundary is already a correct absorbing outflow, and adding a custom deleting
    model would only duplicate it. tests/openfoam asserts this behaviour rather
    than assuming it.

    Dictionary
    ----------
    In constant/dsmcProperties:

    \verbatim
    InflowBoundaryModel   plumeFieldInflow;

    plumeFieldInflowCoeffs
    {
        // Patches to inject across. Required, and every name must exist.
        patches ( inflow );

        // species -> volScalarField holding the per-face number density [1/m^3].
        numberDensityFields
        {
            N2    boundaryNumberDensity_N2;
        }
    }
    \endverbatim

    Every failure is fatal and named. There is no fallback to a uniform density
    and no silent skipping of a missing patch: a run that quietly degraded to
    stock FreeStream behaviour would produce plausible-looking output for the
    wrong boundary condition, which is the failure mode this whole class exists
    to remove.

SourceFiles
    plumeFieldInflow.C

\*---------------------------------------------------------------------------*/

#ifndef plumeFieldInflow_H
#define plumeFieldInflow_H

#include "InflowBoundaryModel.H"
#include "polyMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plumeFieldInflow Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class plumeFieldInflow
:
    public InflowBoundaryModel<CloudType>
{
    // Private data

        //- Indices of the patches to inject across, from the coeffs dictionary.
        //  Never derived from the patch type: that is the whole point.
        labelList patches_;

        //- Patch names, kept for error messages after the indices are resolved.
        wordList patchNames_;

        //- The molecule types to be introduced
        List<label> moleculeTypeIds_;

        //- Per-face number density on each injection patch, already divided by
        //  nParticle so it is in PARCELS per m^3 -- which is what the flux
        //  accumulator counts.
        //  + Outer List  - one inner List per patch, matching patches_
        //  + Inner List  - one Field per species
        //  + Field entry - one per patch face
        //
        //  Extracted ONCE at construction rather than read from the field every
        //  timestep, and copied out rather than scaling the volScalarField in
        //  place. Scaling in place would touch the patch fields of boundaries
        //  this model has nothing to do with, including zeroGradient ones whose
        //  storage OpenFOAM has not initialised at that point -- arithmetic on
        //  which raises a floating-point exception. It also leaves the field on
        //  disk in physical 1/m^3, which is what anyone inspecting it expects.
        List<List<Field<scalar>>> numberDensities_;

        //- The number-density fields, kept registered and AUTO_WRITE so they are
        //  rewritten into each output time directory and a restart finds them.
        //  Not otherwise used after construction.
        PtrList<volScalarField> numberDensityFields_;

        //- Names of those fields, for error messages.
        wordList numberDensityFieldNames_;

        //- Carry-over of fractional particle flux from one timestep to the next.
        //  + Outer List  - one inner List per patch
        //  + Inner List  - one Field per species
        //  + Field entry - one per patch face
        List<List<Field<scalar>>> particleFluxAccumulators_;


    // Private Member Functions

        //- Resolve the configured patch names to indices, or abort.
        void readPatches(const CloudType& cloud);

        //- Read the per-species number-density fields, or abort.
        void readNumberDensityFields(CloudType& cloud);

        //- Copy each injection patch's number density out and scale it by
        //  1/nParticle, or abort.
        void extractNumberDensities(CloudType& cloud);


public:

    //- Runtime type information
    TypeName("plumeFieldInflow");


    // Constructors

        //- Construct from dictionary
        plumeFieldInflow
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~plumeFieldInflow() = default;


    // Member Functions

        // Mapping

            //- Remap the accumulators following a mesh change
            virtual void autoMap(const mapPolyMesh&);

        //- Introduce particles
        virtual void inflow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "plumeFieldInflow.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | rarefiedPlume: plume-impingement DSMC cases
   \\    /   O peration     |
    \\  /    A nd           | https://github.com/andytorrestb/rarefiedPlume
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later.

Description
    Runtime-selection registration.

    Instantiating plumeFieldInflow for DSMCCloud<dsmcParcel> and adding it to the
    InflowBoundaryModel selection table is what makes

        InflowBoundaryModel   plumeFieldInflow;

    resolvable in constant/dsmcProperties. The table itself is already defined by
    OpenFOAM's own makeDSMCParcelInflowBoundaryModels.C, so only the entry is
    added here -- calling makeInflowBoundaryModel() again would define the table
    a second time and the link would fail.

    dsmcFoam picks the library up through controlDict:

        libs ( "libplumeDsmcBoundaryModels.so" );

    If that line is missing, OpenFOAM reports plumeFieldInflow as an unknown
    InflowBoundaryModel and lists the models it does know. That is the intended
    behaviour: the run stops rather than falling back to a uniform inflow.

\*---------------------------------------------------------------------------*/

#include "dsmcParcel.H"
#include "DSMCCloud.H"

#include "plumeFieldInflow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef DSMCCloud<dsmcParcel> CloudType;

    makeInflowBoundaryModelType(plumeFieldInflow, CloudType);
}


// ************************************************************************* //

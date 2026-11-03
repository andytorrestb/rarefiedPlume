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
    One extra run-time-selectable psiThermo:

        hePsiThermo<pureMixture<tabulatedTransport<species::thermo<
            hConstThermo<perfectGas<specie>>, sensibleInternalEnergy>>>>

    i.e. a calorically perfect ideal gas whose viscosity and conductivity are
    read from mu(T), kappa(T) tables. OpenFOAM v2512 compiles tabulated
    transport for rhoThermo (heRhoThermo) but not for psiThermo, and
    rhoCentralFoam uses psiThermo. Without this, the only temperature-dependent
    transport rhoCentralFoam can use is Sutherland's law, which cannot represent
    the VHS power law mu = mu_ref (T/T_ref)^omega that the DSMC collision model
    implies over the decade-plus temperature range of a plume expansion.

    Nothing here is new physics: tabulatedTransport, hConstThermo, perfectGas
    and hePsiThermo are OpenFOAM's own; this file only instantiates the
    combination and registers it in psiThermo's selection table, exactly as
    src/thermophysicalModels/basic/psiThermo/psiThermos.C does for the others.

    The tables are written by plumetools.hybridplume.cfd from the same Gas
    object the DSMC dictionaries come from. Select with

        thermoType
        {
            type            hePsiThermo;
            mixture         pureMixture;
            transport       tabulated;
            thermo          hConst;
            equationOfState perfectGas;
            specie          specie;
            energy          sensibleInternalEnergy;
        }

    and load with  libs ("libhybridPlumeThermo.so");  in controlDict.

\*---------------------------------------------------------------------------*/

#include "psiThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "tabulatedTransport.H"

#include "hePsiThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeThermos
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | rarefiedPlume: plume-impingement DSMC cases
   \\    /   O peration     |
    \\  /    A nd           | https://github.com/andytorrestb/rarefiedPlume
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later.

Class
    Foam::cappedTabulatedTransport

Description
    tabulatedTransport with the kinematic viscosity capped at nuMax:

        mu_eff    = min(mu_tab(T), nuMax rho(p, T))
        kappa_eff = kappa_tab(T) mu_eff/mu_tab(T)
        alphah    = kappa_eff/Cp

    where rho(p, T) comes from the equation of state (p/(R T) for a perfect
    gas). Prandtl number is preserved. rhoCentralFoam advances viscous work
    explicitly, so in the near-vacuum cells beside a plume the viscous number
    mu/(rho h^2) dt grows without bound; capping nu = mu/rho bounds it. The
    cap is meant to act only far outside the retained continuum region, and
    plumetools.hybridplume checks that it does not act inside it.

Usage
    \verbatim
    transport
    {
        mu      ( (T0 mu0) (T1 mu1) ... );
        kappa   ( (T0 k0)  (T1 k1)  ... );
        nuMax   1.0e-02;    // [m^2/s]
    }
    \endverbatim

    Selected as  transport cappedTabulated;  in thermoType.

\*---------------------------------------------------------------------------*/

#ifndef cappedTabulatedTransport_H
#define cappedTabulatedTransport_H

#include "tabulatedTransport.H"
#include "IOstreams.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Thermo> class cappedTabulatedTransport;

template<class Thermo>
Ostream& operator<<(Ostream&, const cappedTabulatedTransport<Thermo>&);


/*---------------------------------------------------------------------------*\
                  Class cappedTabulatedTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class cappedTabulatedTransport
:
    public tabulatedTransport<Thermo>
{
    // Private Data

        //- Maximum kinematic viscosity [m^2/s]
        scalar nuMax_;

        //- The input tables, kept only for write(): tabulatedTransport's
        //  own copies are private and it has no accessors
        List<Tuple2<scalar, scalar>> muValues_;
        List<Tuple2<scalar, scalar>> kappaValues_;


    // Private Member Functions

        //- Ratio mu_eff/mu_tab in [0, 1]
        inline scalar capRatio(const scalar p, const scalar T) const
        {
            const scalar muTab = tabulatedTransport<Thermo>::mu(p, T);
            const scalar muCap = nuMax_*this->rho(p, T);
            return (muCap < muTab) ? muCap/muTab : 1.0;
        }


public:

    // Constructors

        //- Construct as named copy
        inline cappedTabulatedTransport
        (
            const word& name,
            const cappedTabulatedTransport& pt
        )
        :
            tabulatedTransport<Thermo>(name, pt),
            nuMax_(pt.nuMax_),
            muValues_(pt.muValues_),
            kappaValues_(pt.kappaValues_)
        {}

        //- Construct from dictionary
        explicit cappedTabulatedTransport(const dictionary& dict)
        :
            tabulatedTransport<Thermo>(dict),
            nuMax_(dict.subDict("transport").get<scalar>("nuMax")),
            muValues_
            (
                dict.subDict("transport")
                    .get<List<Tuple2<scalar, scalar>>>("mu")
            ),
            kappaValues_
            (
                dict.subDict("transport")
                    .get<List<Tuple2<scalar, scalar>>>("kappa")
            )
        {
            if (!(nuMax_ > 0))
            {
                FatalIOErrorInFunction(dict)
                    << "nuMax must be positive, got " << nuMax_
                    << exit(FatalIOError);
            }
        }

        //- Return a clone
        inline autoPtr<cappedTabulatedTransport> clone() const
        {
            return autoPtr<cappedTabulatedTransport<Thermo>>::New(*this);
        }

        //- Selector from dictionary
        inline static autoPtr<cappedTabulatedTransport> New
        (
            const dictionary& dict
        )
        {
            return autoPtr<cappedTabulatedTransport<Thermo>>::New(dict);
        }


    // Member Functions

        //- The instantiated type name
        static word typeName()
        {
            return "cappedTabulated<" + Thermo::typeName() + '>';
        }

        //- Maximum kinematic viscosity [m^2/s]
        scalar nuMax() const
        {
            return nuMax_;
        }

        //- Dynamic viscosity, capped at nuMax rho [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const
        {
            return min
            (
                tabulatedTransport<Thermo>::mu(p, T),
                nuMax_*this->rho(p, T)
            );
        }

        //- Thermal conductivity, scaled with the viscosity cap [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const
        {
            return tabulatedTransport<Thermo>::kappa(p, T)*capRatio(p, T);
        }

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const
        {
            return kappa(p, T)/this->Cp(p, T);
        }

        //- Write to Ostream
        void write(Ostream& os) const
        {
            // tabulatedTransport::write closes its own blocks, so re-emit the
            // whole entry here with nuMax added to the transport block.
            // Thermo::write is species::thermo's, as in tabulatedTransport.
            os.beginBlock(this->name());
            Thermo::write(os);
            {
                os.beginBlock("transport");
                os.writeEntry("mu", muValues_);
                os.writeEntry("kappa", kappaValues_);
                os.writeEntry("nuMax", nuMax_);
                os.endBlock();
            }
            os.endBlock();
        }


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const cappedTabulatedTransport&
        );
};


template<class Thermo>
Ostream& operator<<(Ostream& os, const cappedTabulatedTransport<Thermo>& pt)
{
    pt.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

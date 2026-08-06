/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | rarefiedPlume: plume-impingement DSMC cases
   \\    /   O peration     |
    \\  /    A nd           | https://github.com/andytorrestb/rarefiedPlume
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    GPL-3.0-or-later. Derived from OpenFOAM's FreeStream inflow boundary model,
    Copyright (C) 2011-2017 OpenFOAM Foundation, (C) 2022 OpenCFD Ltd.

    The injection algorithm below -- the Bird eqn 4.22 flux accumulator, the
    eqn 12.5 normal-velocity sampling, the triangle-area-weighted face position
    and the equipartition internal energy -- is FreeStream's, kept deliberately
    line-for-line so that a future OpenFOAM change to it can be diffed in.

    Exactly two things differ, both marked CHANGED below:
      1. the patch list comes from the dictionary, not from isType<polyPatch>;
      2. the number density is read per face from a volScalarField, not as one
         scalar per species.

\*---------------------------------------------------------------------------*/

#include "plumeFieldInflow.H"
#include "constants.H"
#include "tetIndices.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::plumeFieldInflow<CloudType>::readPatches(const CloudType& cloud)
{
    const polyBoundaryMesh& boundary = cloud.mesh().boundaryMesh();

    // CHANGED (1 of 2): an explicit list, where FreeStream walks every patch and
    // keeps the ones whose type is exactly polyPatch.
    this->coeffDict().readEntry("patches", patchNames_);

    if (patchNames_.empty())
    {
        FatalErrorInFunction
            << "plumeFieldInflowCoeffs/patches is empty, so no particles would"
            << " ever be injected." << nl
            << "This model exists to name the inflow patch explicitly; an empty"
            << " list is always a configuration mistake." << nl
            << abort(FatalError);
    }

    patches_.setSize(patchNames_.size());

    forAll(patchNames_, i)
    {
        const label patchi = boundary.findPatchID(patchNames_[i]);

        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Patch " << patchNames_[i] << " named in"
                << " plumeFieldInflowCoeffs/patches is not in the mesh." << nl
                << "Available patches: " << boundary.names() << nl
                << abort(FatalError);
        }

        // A body must not be an inflow. Injecting across a wall would both
        // create particles inside a solid and suppress the surface fluxes
        // hitWallPatch records, and the symptom -- a body with no drag -- is a
        // long way from the cause.
        if (isA<wallPolyPatch>(boundary[patchi]))
        {
            FatalErrorInFunction
                << "Patch " << patchNames_[i] << " has geometric type 'wall',"
                << " but is named as an inflow patch." << nl
                << "A wall is a gas-surface interaction site, not a source."
                << nl << abort(FatalError);
        }

        patches_[i] = patchi;
    }

    Info<< "    plumeFieldInflow: injecting across " << patchNames_
        << " (patch ids " << patches_ << ")" << endl;
}


template<class CloudType>
void Foam::plumeFieldInflow<CloudType>::readNumberDensityFields(CloudType& cloud)
{
    const fvMesh& mesh = cloud.mesh();

    // CHANGED (2 of 2): FreeStream reads
    //     numberDensities_[i] = numberDensitiesDict.get<scalar>(molecules[i]);
    // -- one scalar per species. Here each species names a volScalarField whose
    // BOUNDARY values carry the number density face by face.
    const dictionary& fieldsDict
    (
        this->coeffDict().subDict("numberDensityFields")
    );

    const wordList molecules(fieldsDict.toc());

    if (molecules.empty())
    {
        FatalErrorInFunction
            << "plumeFieldInflowCoeffs/numberDensityFields is empty." << nl
            << abort(FatalError);
    }

    moleculeTypeIds_.setSize(molecules.size());
    numberDensityFieldNames_.setSize(molecules.size());
    numberDensityFields_.setSize(molecules.size());

    forAll(molecules, i)
    {
        moleculeTypeIds_[i] = cloud.typeIdList().find(molecules[i]);

        if (moleculeTypeIds_[i] == -1)
        {
            FatalErrorInFunction
                << "typeId " << molecules[i] << " is not defined in the cloud."
                << nl << "typeIdList: " << cloud.typeIdList() << nl
                << abort(FatalError);
        }

        const word fieldName(fieldsDict.get<word>(molecules[i]));
        numberDensityFieldNames_[i] = fieldName;

        // MUST_READ and AUTO_WRITE, exactly as DSMCCloud constructs boundaryT
        // and boundaryU. AUTO_WRITE is what makes a restart work: the field is
        // rewritten into each output time directory, so a run resumed from
        // t > 0 finds it there rather than only in 0/.
        numberDensityFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                mesh
            )
        );

        Info<< "    plumeFieldInflow: " << molecules[i] << " number density from "
            << fieldName << endl;
    }
}


template<class CloudType>
void Foam::plumeFieldInflow<CloudType>::extractNumberDensities(CloudType& cloud)
{
    const scalar nParticle = cloud.nParticle();

    if (nParticle <= 0)
    {
        FatalErrorInFunction
            << "nEquivalentParticles is " << nParticle
            << "; it must be positive." << nl << abort(FatalError);
    }

    numberDensities_.setSize(patches_.size());

    forAll(patches_, p)
    {
        const label patchi = patches_[p];

        numberDensities_[p].setSize(moleculeTypeIds_.size());

        forAll(moleculeTypeIds_, i)
        {
            // Copy the patch's own values out, then scale. FreeStream applies
            // the same nParticle division once at construction:
            //     numberDensities_ /= cloud.nParticle();
            // but its numberDensities_ is a plain list of scalars. Doing the
            // equivalent to a volScalarField would scale every patch field in
            // the mesh, including ones this model never reads.
            numberDensities_[p][i] =
                numberDensityFields_[i].boundaryField()[patchi];

            numberDensities_[p][i] /= nParticle;

            if (min(numberDensities_[p][i]) < 0)
            {
                FatalErrorInFunction
                    << "Negative number density on inflow patch "
                    << patchNames_[p] << " in field "
                    << numberDensityFieldNames_[i] << "." << nl
                    << "Range: " << min(numberDensities_[p][i])*nParticle
                    << " .. " << max(numberDensities_[p][i])*nParticle
                    << " 1/m^3" << nl
                    << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::plumeFieldInflow<CloudType>::plumeFieldInflow
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InflowBoundaryModel<CloudType>(dict, cloud, typeName),
    patches_(),
    patchNames_(),
    moleculeTypeIds_(),
    numberDensities_(),
    numberDensityFields_(),
    numberDensityFieldNames_(),
    particleFluxAccumulators_()
{
    readPatches(cloud);
    readNumberDensityFields(cloud);
    extractNumberDensities(cloud);

    // One accumulator field per (patch, species), sized to the patch's LOCAL
    // face count. In a decomposed run every rank holds every patch, most with
    // zero faces, so this is correct in serial and in parallel without a special
    // case: a rank owning no inflow faces simply accumulates nothing.
    particleFluxAccumulators_.setSize(patches_.size());

    forAll(patches_, p)
    {
        const polyPatch& patch = cloud.mesh().boundaryMesh()[patches_[p]];

        particleFluxAccumulators_[p] = List<Field<scalar>>
        (
            moleculeTypeIds_.size(),
            Field<scalar>(patch.size(), Zero)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::plumeFieldInflow<CloudType>::autoMap(const mapPolyMesh&)
{
    const CloudType& cloud(this->owner());
    const polyMesh& mesh(cloud.mesh());

    forAll(patches_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patches_[p]];
        List<Field<scalar>>& pFA = particleFluxAccumulators_[p];

        forAll(pFA, i)
        {
            pFA[i].setSize(patch.size(), 0);
        }
    }
}


template<class CloudType>
void Foam::plumeFieldInflow<CloudType>::inflow()
{
    CloudType& cloud(this->owner());

    const polyMesh& mesh(cloud.mesh());

    const scalar deltaT = mesh.time().deltaTValue();

    Random& rndGen = cloud.rndGen();

    const scalar sqrtPi = sqrt(pi);

    label particlesInserted = 0;

    const volScalarField::Boundary& boundaryT
    (
        cloud.boundaryT().boundaryField()
    );

    const volVectorField::Boundary& boundaryU
    (
        cloud.boundaryU().boundaryField()
    );


    forAll(patches_, p)
    {
        const label patchi = patches_[p];

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.empty())
        {
            continue;
        }

        // Add mass to the accumulators. Negative face area dotted with the
        // velocity to point the flux INTO the domain.
        List<Field<scalar>>& pFA = particleFluxAccumulators_[p];

        forAll(pFA, i)
        {
            const label typeId = moleculeTypeIds_[i];

            const scalar mass = cloud.constProps(typeId).mass();

            if (min(boundaryT[patchi]) < SMALL)
            {
                FatalErrorInFunction
                    << "Zero boundary temperature on inflow patch "
                    << patch.name() << ", check the boundaryT condition." << nl
                    << abort(FatalError);
            }

            // CHANGED: per-face, where FreeStream has one scalar for the patch.
            // Already divided by nParticle at construction, so this is in
            // parcels per m^3 -- the same units FreeStream's numberDensities_
            // carries by the time it reaches this expression.
            const scalarField& numberDensity = numberDensities_[p][i];

            const scalarField mostProbableSpeed
            (
                cloud.maxwellianMostProbableSpeed
                (
                    boundaryT[patchi],
                    mass
                )
            );

            // Boundary velocity dotted with the face unit normal (which points
            // out of the domain, hence the negation), over the most probable
            // speed: the molecular speed ratio times cos(theta).
            const scalarField sCosTheta
            (
                (boundaryU[patchi] & -patch.faceAreas()/mag(patch.faceAreas()))
              / mostProbableSpeed
            );

            // Bird eqn 4.22: the half-range flux of a drifting Maxwellian
            // through a surface. Unchanged from FreeStream except that
            // numberDensity is now a field rather than a scalar, so this is an
            // element-wise product instead of a scalar one.
            pFA[i] +=
                mag(patch.faceAreas())*numberDensity*deltaT
               *mostProbableSpeed
               *(
                   exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
               /(2.0*sqrtPi);
        }

        forAll(patch, pFI)
        {
            // Loop over faces as the outer loop to avoid recomputing the
            // geometric properties for every species.

            const face& f = patch[pFI];

            const label globalFaceIndex = pFI + patch.start();

            const label celli = mesh.faceOwner()[globalFaceIndex];

            const vector fC = patch.faceCentres()[pFI];

            const scalar fA = mag(patch.faceAreas()[pFI]);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh,
                globalFaceIndex,
                celli
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), Zero);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh).mag()/fA
                  + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last fraction to 1.0 so rounding on a non-flat face
            // cannot leave it below 1.
            cTriAFracs.last() = 1.0;

            // Normal unit vector NEGATED, so it points into the domain
            const vector n = -normalised(patch.faceAreas()[pFI]);

            // Wall tangential unit vector, from the face centre to its first
            // vertex
            const vector t1 = normalised(fC - (mesh.points()[f[0]]));

            // The other tangential unit vector, rescaled in case the face is not
            // flat and n and t1 are not perfectly orthogonal
            const vector t2 = normalised(n ^ t1);

            const scalar faceTemperature = boundaryT[patchi][pFI];

            const vector& faceVelocity = boundaryU[patchi][pFI];

            forAll(pFA, i)
            {
                scalar& faceAccumulator = pFA[i][pFI];

                // Number of whole particles to insert
                label nI = max(label(faceAccumulator), 0);

                // Add one more with a probability equal to the remainder, so the
                // long-run mean flux is exact rather than truncated.
                if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
                {
                    nI++;
                }

                faceAccumulator -= nI;

                const label typeId = moleculeTypeIds_[i];

                const scalar mass = cloud.constProps(typeId).mass();

                for (label j = 0; j < nI; j++)
                {
                    // Choose a triangle to insert on, weighted by area
                    const scalar triSelection = rndGen.sample01<scalar>();

                    label selectedTriI = -1;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    const tetIndices& faceTetIs = faceTets[selectedTriI];

                    point pos = faceTetIs.faceTri(mesh).randomPoint(rndGen);

                    // Velocity generation

                    const scalar mostProbableSpeed
                    (
                        cloud.maxwellianMostProbableSpeed
                        (
                            faceTemperature,
                            mass
                        )
                    );

                    const scalar sCosTheta =
                        (faceVelocity & n)/mostProbableSpeed;

                    // Coefficients for Bird eqn 12.5
                    const scalar uNormProbCoeffA =
                        sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                    const scalar uNormProbCoeffB =
                        0.5*
                        (
                            1.0
                          + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                        );

                    // Equivalent to the QA value in Bird's DSMC3.FOR
                    scalar randomScaling = 3.0;

                    if (sCosTheta < -3)
                    {
                        randomScaling = mag(sCosTheta) + 1;
                    }

                    scalar P = -1;

                    scalar uNormal;
                    scalar uNormalThermal;

                    // Acceptance-rejection on Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());

                    const vector U =
                        sqrt(physicoChemical::k.value()*faceTemperature/mass)
                       *(
                            rndGen.GaussNormal<scalar>()*t1
                          + rndGen.GaussNormal<scalar>()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    const scalar Ei = cloud.equipartitionInternalEnergy
                    (
                        faceTemperature,
                        cloud.constProps(typeId).internalDegreesOfFreedom()
                    );

                    cloud.addNewParcel(pos, celli, U, Ei, typeId);

                    particlesInserted++;
                }
            }
        }
    }

    Info<< "    Particles inserted              = "
        << returnReduce(particlesInserted, sumOp<label>()) << endl;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Simple central-difference snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "skewCorrectedSnGrad.H"
#include "fv.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "correctedSnGrad.H"
#include "uncorrectedSnGrad.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
tmp<GeometricField<vector, fvsPatchField, surfaceMesh> >
skewCorrectedSnGrad<vector>::correction
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<vector, fvsPatchField, surfaceMesh>
    tmp<GeometricField<vector, fvsPatchField, surfaceMesh> > tssf
    (
        new GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "snGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    GeometricField<vector, fvsPatchField, surfaceMesh>& ssf = tssf();

    ssf = dimensioned<vector>("0", ssf.dimensions(), vector::zero);

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const vectorField& Sf = mesh.Sf().internalField();
//     const scalarField& magSf = mesh.magSf().internalField();

    vectorField nf = Sf/mag(Sf);

    const vectorField& Cf = mesh.Cf().internalField();
    const vectorField& C = mesh.C().internalField();

    const scalarField& deltaCoeffs = 
        mesh.deltaCoeffs().internalField();

    surfaceVectorField kP
    (
        IOobject
        (
            "kP",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    );
    vectorField& kPI = kP.internalField();

    surfaceVectorField kN
    (
        IOobject
        (
            "kN",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, vector::zero)
    );
    vectorField& kNI = kN.internalField();

    kPI = Cf - vectorField(C, owner);
    kPI -= Sf*(Sf&kPI)/magSqr(Sf);
//     kPI -= Sf*(Sf&kPI)/sqr(magSf);

    kNI = Cf - vectorField(C, neighbour);
    kNI -= Sf*(Sf&kNI)/magSqr(Sf);
//     kNI -= Sf*(Sf&kNI)/sqr(magSf);

    // Ununiformity correction
    if (false)
    {
        vectorField delta = 
            Cf 
          - (vectorField(C, neighbour) + kN + vectorField(C, owner) + kP)/2.0;

        kPI += delta;
        kNI += delta;
    }

    forAll(kP.boundaryField(), patchI)
    {
        if (kP.boundaryField()[patchI].coupled())
        {
            kP.boundaryField()[patchI] = 
                mesh.boundary()[patchI].Cf()
              - mesh.boundary()[patchI].Cn();

            kP.boundaryField()[patchI] -= 
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kP.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

            kN.boundaryField()[patchI] = 
                mesh.Cf().boundaryField()[patchI]
              - (
                    mesh.boundary()[patchI].Cn()
                  + mesh.boundary()[patchI].delta()
                );

            kN.boundaryField()[patchI] -= 
                mesh.boundary()[patchI].Sf()
               *(
                    mesh.boundary()[patchI].Sf()
                  & kN.boundaryField()[patchI]
                )
               /sqr(mesh.boundary()[patchI].magSf());

//             Info << mesh.boundary()[patchI].name() << endl;
//             Info << 1.0/mesh.boundary()[patchI].deltaCoeffs() << endl;

//             // Ununiformity corr
//             vectorField delta = 
//                 mesh.boundary()[patchI].Cf()
//               - (
//                     (
//                         mesh.boundary()[patchI].Cn()
//                       + mesh.boundary()[patchI].delta()
//                     )
//                   + kN.boundaryField()[patchI]
//                   + mesh.boundary()[patchI].Cn()
//                   + kP.boundaryField()[patchI]
//                 )/2.0;

//             kP.boundaryField()[patchI] += delta;
//             kN.boundaryField()[patchI] += delta;
        }
    }

//     volTensorField gradVf =
//         gradScheme<vector>::New
//         (
//             mesh,
//             mesh.schemesDict().gradScheme(ssf.name())
//         )()
//        .grad(vf);

    const volTensorField* gradVfPtr(NULL);

    bool calcGradient = false;

    if ( mesh.found("grad(" + vf.name() + ")") )
    {
//         Info << "Use stored gradient " << vf.name() << endl;

        const volTensorField& gradVf_ = 
            mesh.lookupObject<volTensorField>("grad(" + vf.name() + ")");

        gradVfPtr = &gradVf_;
    }
    else
    {
//         Info << "Calc gradient " << vf.name() << endl;
        gradVfPtr =
            new volTensorField
            (
                gradScheme<vector>::New
                (
                    mesh,
                    mesh.schemesDict().gradScheme(ssf.name())
                )()
               .grad(vf)
            );

        calcGradient = true;
    }

    const volTensorField& gradVf = *gradVfPtr;

//     Info << "magGrad, avg: " << gMax(mag(gradVf.internalField()))
//         << "min: " << gMin(mag(gradVf.internalField()))
//         << "avg: " << gAverage(mag(gradVf.internalField())) << endl;

//     const volTensorField& gradVf = 
//         mesh.lookupObject<volTensorField>("grad(" + vf.name() + ")");

    ssf.internalField() = 
        (
            (kNI&Field<tensor>(gradVf, neighbour))
          - (kPI&Field<tensor>(gradVf, owner))
        )
       *deltaCoeffs;

    // Add third order correction
    if (false)
    {
//         Info << "Third order correction" << endl;
        const scalarField& weights = 
            mesh.weights().internalField();
        const vectorField& psiI = vf.internalField();
        const tensorField& gradPsiI = gradVf.internalField();

        forAll(ssf.internalField(), faceI)
        {
            scalar d = 1.0/deltaCoeffs[faceI];
            scalar w = weights[faceI];
            scalar df = (1-w)*d;

            vector psiP = 
                psiI[owner[faceI]] 
              + (kPI[faceI]&gradPsiI[owner[faceI]]);
            vector psiN = 
                psiI[neighbour[faceI]] 
              + (kNI[faceI]&gradPsiI[neighbour[faceI]]);

            vector dPsiP = (nf[faceI]&gradPsiI[owner[faceI]]);
            vector dPsiN = (nf[faceI]&gradPsiI[neighbour[faceI]]);

            vector a = 3*(psiN - psiP)/sqr(d) - (dPsiN + 2*dPsiP)/d;
            vector b = (dPsiN + dPsiP)/sqr(d) - 2*(psiN - psiP)/(d*sqr(d));

            vector dPsif = dPsiP + 2*a*df + 3*b*sqr(df);

            ssf.internalField()[faceI] += dPsif - (psiN - psiP)/d;
        }
    }

    forAll(ssf.boundaryField(), patchI)
    {
        if (ssf.boundaryField()[patchI].coupled())
        {
            ssf.boundaryField()[patchI] =
            (
                (
                    (
                        kN.boundaryField()[patchI]
                      & gradVf.boundaryField()[patchI].patchNeighbourField()
                    )
                  - (
                        kP.boundaryField()[patchI]
                      & gradVf.boundaryField()[patchI].patchInternalField()
                    )
                )
               *mesh.deltaCoeffs().boundaryField()[patchI]
            );                
        }
    }

    surfaceScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag
            (
                uncorrectedSnGrad<vector>::snGrad
                (
                    vf, 
                    this->deltaCoeffs(vf), 
                    "orthSnGrad"
                )
              + ssf
            )
           /(
                (1 - limitCoeff_)*mag(ssf)
              + dimensionedScalar("small", ssf.dimensions(), SMALL)
            ),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    if (fv::debug)
    {
        Info<< "limitedSnGrad :: limiter min: " << min(limiter.internalField())
            << " max: "<< max(limiter.internalField())
            << " avg: " << average(limiter.internalField()) << endl;
    }

    ssf *= limiter;

    if (calcGradient)
    {
        deleteDemandDrivenData(gradVfPtr);
    }

    return tssf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

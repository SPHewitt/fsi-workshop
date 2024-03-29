/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "consistentIcoFluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "adjustPhi.H"

#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

#include "elasticSlipWallVelocityFvPatchVectorField.H"
#include "elasticWallVelocityFvPatchVectorField.H"
#include "elasticWallPressureFvPatchScalarField.H"

#include "findRefCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(consistentIcoFluid, 0);
addToRunTimeSelectionTable(fluidSolver, consistentIcoFluid, dictionary);
// addToRunTimeSelectionTable(icoFluid, consistentIcoFluid, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void consistentIcoFluid::makeSf() const
{
    // Find global face zones
    if (SfPtr_)
    {
        FatalErrorIn
        (
            "void fluidSolver::makeSf() const"
        )
            << "Face surface vectors alrady created"
                << abort(FatalError);
    }

//     const fvMesh& mesh = fluidSolver::mesh();

    IOobject SfHeader
    (
        "Sf",
        runTime().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    SfPtr_ =
        new surfaceVectorField
        (
            SfHeader,
            mesh(),
            dimensionedVector("0", dimArea, vector::zero)
        );
    surfaceVectorField& Sf = *SfPtr_;

    if(!SfHeader.headerOk())
    {
        const vectorField& allFaceAreas = mesh().faceAreas();

        Sf.internalField() =
            vectorField::subField(allFaceAreas, mesh().nInternalFaces());

        const fvPatchList& patches = mesh().boundary();

        forAll (patches, patchI)
        {
            Sf.boundaryField()[patchI] =
                patches[patchI].patchSlice(allFaceAreas);
        }
    }
}


void consistentIcoFluid::updateSf()
{
    Sf().oldTime();

    const vectorField& allFaceAreas = mesh().faceAreas();

    Sf().internalField() =
        vectorField::subField(allFaceAreas, mesh().nInternalFaces());

    const fvPatchList& patches = mesh().boundary();

    forAll (patches, patchI)
    {
        Sf().boundaryField()[patchI] =
            patches[patchI].patchSlice(allFaceAreas);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

consistentIcoFluid::consistentIcoFluid(const fvMesh& mesh)
:
    icoFluid(mesh),
    SfPtr_(NULL)
{
    phi().oldTime();
    updateSf();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void consistentIcoFluid::evolve()
{
    Info << "Evolving fluid solver: " << this->type() << endl;

    const fvMesh& mesh = fluidSolver::mesh();

    updateSf();

    int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

    int nNonOrthCorr =
        readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), fluidProperties(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p().name());
    mesh.schemesDict().setFluxRequired(phi().name());

    if(mesh.moving())
    {
        // Make the fluxes relative
        phi() -= fvc::meshPhi(U());
    }

    // CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;

        if (mesh.nInternalFaces())
        {
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*mag(phi());

            CoNum = max(SfUfbyDelta/mesh.magSf())
                .value()*runTime().deltaT().value();

            meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
                .value()*runTime().deltaT().value();

            velMag = max(mag(phi())/mesh.magSf()).value();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum
            << " velocity magnitude: " << velMag << endl;
    }

    fvVectorMatrix UEqn
    (
        fvm::ddt(U())
      + fvm::div(phi(), U())
      - fvm::laplacian(nu(), U())
    );

    solve(UEqn == -gradp());

    // --- PISO loop

    volScalarField AU = UEqn.A();
    surfaceScalarField rAUf("rAUf", 1.0/fvc::interpolate(AU, "interpolate(U)"));

    for (int corr=0; corr<nCorr; corr++)
    {
        volVectorField HU = UEqn.H();

        U() = HU/AU;

        // U() = UEqn.H()/AU;
        // phi() = (fvc::interpolate(U())&mesh.Sf());

#       include "calcPhi.H"

#       include "updateRobinFsiInterface.H"

        adjustPhi(phi(), U(), p());

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    rAUf, p(), "laplacian((1|A(U)),p)"
                )
             == fvc::div(phi())
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            gradp() = fvc::grad(p());

            if (nonOrth == nNonOrthCorr)
            {
                phi() -= pEqn.flux();
            }
        }

        // Continuity error
        {
            volScalarField contErr = fvc::div(phi());

            scalar sumLocalContErr = runTime().deltaT().value()*
                mag(contErr)().weightedAverage(mesh.V()).value();

            scalar globalContErr = runTime().deltaT().value()*
                contErr.weightedAverage(mesh.V()).value();

            Info<< "time step continuity errors : sum local = "
                << sumLocalContErr << ", global = " << globalContErr << endl;
        }

        U() -= gradp()/AU;
        U().correctBoundaryConditions();

        gradU() = fvc::grad(U());
    }
}

const surfaceVectorField& consistentIcoFluid::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


surfaceVectorField& consistentIcoFluid::Sf()
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidSolvers
} // End namespace Foam

// ************************************************************************* //

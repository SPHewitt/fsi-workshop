/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "dynamicFreeSurfaceFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "tetMotionSolver.H"
#include "laplaceTetMotionSolver.H"
#include "velocityLaplacianFvMotionSolver.H"
#include "fixedValueTetPolyPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "transformField.H"
#include "areaFields.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFreeSurfaceFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh, 
        dynamicFreeSurfaceFvMesh, 
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicFreeSurfaceFvMesh::clearOut()
{
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
}

void Foam::dynamicFreeSurfaceFvMesh::makeControlPoints()
{
    if (debug)
    {
        Info<< "dynamicFreeSurfaceFvMesh::makeControlPoints() : "
            << "making control points"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (controlPointsPtr_)
    {
        FatalErrorIn("Foam::dynamicFreeSurfaceFvMesh::makeControlPoints()")
            << "control points already exist"
            << abort(FatalError);
    }

    IOobject controlPointsHeader
    (
        "controlPoints",
        this->time().timeName(),
        *this,
        IOobject::MUST_READ
    );

    if (controlPointsHeader.headerOk())
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    this->time().timeName(),
                    *this,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        controlPointsPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "controlPoints",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh_.areaCentres().internalField()
            );
    }
}


void Foam::dynamicFreeSurfaceFvMesh::makeMotionPointsMask()
{
    if (debug)
    {
        Info<< "dynamicFreeSurfaceFvMesh::makeMotionPointsMask() : "
            << "making motion points mask"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (motionPointsMaskPtr_)
    {
        FatalErrorIn("dynamicFreeSurfaceFvMesh::motionPointsMask()")
            << "motion points mask already exists"
            << abort(FatalError);
    }

    motionPointsMaskPtr_ = new labelList
    (
        mesh().boundaryMesh()[freeSurfacePatchID_].nPoints(),
        1
    );
}


void Foam::dynamicFreeSurfaceFvMesh::makeDirections()
{
    if (debug)
    {
        Info<< "Foam::dynamicFreeSurfaceFvMesh::makeDirections() : "
            << "making displacement directions for points and "
            << "control points"
            << endl;
    }


    // It is an error to attempt to recalculate
    // if the pointer is already set
    if 
    (
        pointsDisplacementDirPtr_ ||  
        facesDisplacementDirPtr_
    )
    {
        FatalErrorIn("Foam::dynamicFreeSurfaceFvMesh::makeDirections()")
            << "points and control points displacement directions "
            << "already exists"
            << abort(FatalError);
    }


    pointsDisplacementDirPtr_ = 
        new vectorField
        (
            mesh().boundaryMesh()[freeSurfacePatchID_].nPoints(),
            vector::zero
        );

    facesDisplacementDirPtr_ = 
        new vectorField
        (
            mesh().boundaryMesh()[freeSurfacePatchID_].size(),
            vector::zero
        );

//     if(!normalMotionDir())
    {
//         if(mag(motionDir_) < SMALL)
//         {
//             FatalErrorIn("freeSurface::makeDirections()")
//                 << "Zero motion direction"
//                     << abort(FatalError);
//         }

        vector motionDir_(0, 1, 0);

        facesDisplacementDir() = motionDir_;
        pointsDisplacementDir() = motionDir_;
    }

//     updateDisplacementDirections();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicFreeSurfaceFvMesh::dynamicFreeSurfaceFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    motionPtr_(motionSolver::New(*this)),
    freeSurfacePatchName_
    (
        dynamicMeshCoeffs_.lookup("freeSurfacePatchName")
    ),
    freeSurfacePatchID_(-1),
    aMesh_(*this),
    fixedFreeSurfacePatches_
    (
        dynamicMeshCoeffs_.lookup("fixedFreeSurfacePatches")
    ),
    controlPointsPtr_(NULL),
    motionPointsMaskPtr_(NULL),
    pointsDisplacementDirPtr_(NULL),
    facesDisplacementDirPtr_(NULL),
    deltaH_(),
    H_(),
    V_(),
    W_()
{
    freeSurfacePatchID_ = boundaryMesh().findPatchID(freeSurfacePatchName_);

    if(freeSurfacePatchID_<0)
    {
        FatalErrorIn
        (
            "dynamicFreeSurfaceFvMesh::"
            "dynamicFreeSurfaceFvMesh(const IOobject& io)"
        )
            << "Can't find patch: " << freeSurfacePatchName_
                << exit(FatalError);
    }

    deltaH_.setSize(boundaryMesh()[freeSurfacePatchID_].size(), 0);
    H_.setSize(boundaryMesh()[freeSurfacePatchID_].size(), 0);

    // Mark fixed free surface boundary points 
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn
            (
                "dynamicFreeSurfaceFvMesh::dynamicFreeSurfaceFvMesh(...)"
            )
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the dynamic mesh dictionary"
                    << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicFreeSurfaceFvMesh::~dynamicFreeSurfaceFvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicFreeSurfaceFvMesh::update()
{
    scalar deltaT = time().deltaT().value();

    fvMesh& mesh = *this;

    scalar fsDeltaT = 0.001;

    const surfaceScalarField& phi = 
        mesh.lookupObject<surfaceScalarField>("phi");
    
    scalarField sweptVolCorr = 
        phi.boundaryField()[freeSurfacePatchID_]*fsDeltaT;
    
    const scalarField& Sf = aMesh().S();
    const vectorField& Nf = aMesh().faceAreaNormals().internalField();

    deltaH_ = sweptVolCorr/(Sf*(Nf & facesDisplacementDir()));

    // This part of the code is valid only 
    // for quad mesh at the interface
    forAll(fixedFreeSurfacePatches_, patchI)
    {
        label fixedPatchID = 
            aMesh().boundary().findPatchID
            (
                fixedFreeSurfacePatches_[patchI]
            );

        if(fixedPatchID == -1)
        {
            FatalErrorIn("()")
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the dynamicMesh dictionary"
                    << abort(FatalError);
        }
        
        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaH_[eFaces[edgeI]] *= 2.0;
        }
    }
    
//     IOdictionary transportProperties
//     (
//         IOobject
//         (
//             "transportProperties",
//             mesh.time().constant(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );
//     dimensionedScalar nu(transportProperties.lookup("nu"));

//     const volScalarField& p = 
//         mesh.lookupObject<volScalarField>("p");
//     const scalarField& pfs = 
//         p.boundaryField()[freeSurfacePatchID_];

//     const volVectorField& U = 
//         mesh.lookupObject<volVectorField>("U");

//     vectorField n = mesh.boundary()[freeSurfacePatchID_].nf();
//     scalarField dn = 1.0/mesh.boundary()[freeSurfacePatchID_].deltaCoeffs();

//     // Add modes from previous iteration
//     if (time().timeIndex() > 1)
//     {
//         W_.append(H_);

//         scalarField Fn = 
//             2*nu.value()
//            *(U.boundaryField()[freeSurfacePatchID_].snGrad() & n)
//           - pfs;

//         V_.append(Fn);
//     }


//     if (time().timeIndex() < 10)
//     {
//         scalar fsDeltaT = 0.001;

//         scalarField UnFs = 
//             pfs*dn/(2*nu.value()) 
//           + (
//                 n 
//               & U.boundaryField()[freeSurfacePatchID_].patchInternalField()
//             );

//         scalarField deltaN = UnFs*fsDeltaT;

//         const vectorField& Nf = aMesh().faceAreaNormals().internalField();

//         deltaH_ = deltaN/(Nf & facesDisplacementDir());
//     }
//     else
//     {
//         scalarRectangularMatrix v(pfs.size(), V_.size()-1, 0.0);
//         scalarRectangularMatrix w(pfs.size(), W_.size()-1, 0.0);

//         for (label j=0; j<v.m(); j++)
//         {
//             for (label i=0; i<v.n(); i++)
//             {
//                 v[i][j] = V_[j+1][i] - V_[0][i];
//                 w[i][j] = W_[j+1][i] - W_[0][i];
//             }
//         }

//         scalarSquareMatrix lsv(v.m(), 0.0);

//         for (label i=0; i<v.m(); i++)
//         {
//             for (label j=0; j<v.m(); j++)
//             {
//                 for (label k=0; k<v.n(); k++)
//                 {
//                     lsv[i][j] += v[k][i]*v[k][j];
//                 }
//             }
//         }

//         // Calculate inverse
//         scalarSquareMatrix invlsv = lsv.LUinvert();

//         scalarRectangularMatrix invLsVVt(V_.size()-1, pfs.size(), 0.0);

//         for (label i=0; i<v.m(); i++)
//         {
//             for (label j=0; j<v.n(); j++)
//             {
//                 for (label k=0; k<v.m(); k++)
//                 {
//                     invLsVVt[i][j] += invlsv[i][k]*v[j][k];
//                 }
//             }
//         }

//         scalarRectangularMatrix F(pfs.size(), pfs.size(), 0.0);        
//         multiply(F, w, invLsVVt);

//         scalarField newH = W_[0];

//         for (label i=0; i<pfs.size(); i++)
//         {
//             for (label j=0; j<pfs.size(); j++)
//             {
//                 newH[i] -= F[i][j]*V_[0][j]; 
//             }
//         }

//         deltaH_ = newH - H_;
//     }

    if
    (
        motionPtr_->type()
     == velocityLaplacianFvMotionSolver::typeName
    )
    {
        velocityLaplacianFvMotionSolver& mSolver =
            refCast<velocityLaplacianFvMotionSolver>
            (
                motionPtr_()
            );

        pointVectorField& pointMotionU = mSolver.pointMotionU();
        
        if
        (
            pointMotionU.boundaryField()[freeSurfacePatchID_].type()
         == fixedValuePointPatchVectorField::typeName
        )
        {
            fixedValuePointPatchVectorField& pMotionUFsPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    pointMotionU.boundaryField()[freeSurfacePatchID_]
                );

            pointField displacement = pointDisplacement(deltaH_);

            pMotionUFsPatch == displacement/deltaT;

//             const pointVectorField& pointU =
//                 this->objectRegistry::lookupObject<pointVectorField>
//                 (
//                     "pointU"
//                 );

//             vectorField pointDisplacement =
//                 pointU.boundaryField()[freeSurfacePatchID_]
//                .patchInternalField()*fsDeltaT;

//             const vectorField& curPointNormals = aMesh_.pointAreaNormals();

//             pointDisplacement = 
//                 curPointNormals*(curPointNormals & pointDisplacement);

//             pointDisplacement =
//                 (pointDisplacement & curPointNormals)*pointsDisplacementDir()
//                /(curPointNormals & pointsDisplacementDir());

//             pMotionUFsPatch == pointDisplacement/deltaT;
        }
        else
        {
            FatalErrorIn("dynamicFreeSurfaceFvMesh::update()")
                << "Bounary condition on " << pointMotionU.name()
                    <<  " for " << freeSurfacePatchName_ << " patch is "
                    << pointMotionU.boundaryField()[freeSurfacePatchID_].type()
                    << ", instead "
                    << fixedValuePointPatchVectorField::typeName
                    << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("dynamicFreeSurfaceFvMesh::update()")
            << "Selected mesh motion solver is "
                << motionPtr_->type()
                << ", instead "
                << velocityLaplacianFvMotionSolver::typeName
                << exit(FatalError);
    }

    fvMesh::movePoints(motionPtr_->newPoints());

    aMesh().movePoints();

    // Mesh motion only - return false
    return false;
}


Foam::vectorField& Foam::dynamicFreeSurfaceFvMesh::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}


Foam::labelList& Foam::dynamicFreeSurfaceFvMesh::motionPointsMask()
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}


Foam::vectorField& Foam::dynamicFreeSurfaceFvMesh::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}


Foam::vectorField& Foam::dynamicFreeSurfaceFvMesh::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}


Foam::faMesh& Foam::dynamicFreeSurfaceFvMesh::aMesh()
{
//     if (!aMeshPtr_)
//     {
//         makeFaMesh();
//     }
    
    return aMesh_;
}

// ************************************************************************* //

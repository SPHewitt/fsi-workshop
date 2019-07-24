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

#include "leastSquaresVectorsBiLinear.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "mapPolyMesh.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresVectorsBiLinear, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::leastSquaresVectorsBiLinear::leastSquaresVectorsBiLinear(const fvMesh& mesh)
:
    MeshObject<fvMesh, leastSquaresVectorsBiLinear>(mesh),
    invLsMatrices_(0),
    pVectorsPtr_(NULL),
    nVectorsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::leastSquaresVectorsBiLinear::~leastSquaresVectorsBiLinear()
{
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::leastSquaresVectorsBiLinear::makeLeastSquaresVectors() const
{
    if (debug)
    {
        Info<< "leastSquaresVectorsBiLinear::makeLeastSquaresVectors() :"
            << "Constructing least square gradient vectors"
            << endl;
    }

    pVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresP",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsP = *pVectorsPtr_;

    nVectorsPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "LeastSquaresN",
            mesh().pointsInstance(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    );
    surfaceVectorField& lsN = *nVectorsPtr_;

    emptyPatchPVectorsPtr_ = new vectorField(0, vector::zero);
    vectorField& emptyPatchPVectors = *emptyPatchPVectorsPtr_;

    forAll(mesh().boundary(), patchI)
    {
        const fvPatch& p = mesh().boundary()[patchI];

        if (isA<emptyFvPatch>(p))
        {
            emptyPatchPVectors = 
                vectorField(p.patch().size(), vector::zero);
            break;
        } 
    }

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& C = mesh().cellCentres();
//     const surfaceScalarField& w = mesh().weights();
//     const surfaceScalarField& magSf = mesh().magSf();

    // Set up temporary storage for the dd tensor (before inversion)
//     symmTensorField dd(mesh().nCells(), symmTensor::zero);

    // Set up storage for inverse least-squares normal equation matrix
    invLsMatrices_.setSize(mesh().nCells());
    PtrList<scalarRectangularMatrix> M(mesh().nCells());

    label nCoeffs = 6;
    const cellList& cells = mesh().cells();
    forAll(invLsMatrices_, cellI)
    {
        invLsMatrices_.set
        (
            cellI, 
            new scalarRectangularMatrix
            (
                nCoeffs,
                cells[cellI].size(),
                0.0
            )
        );

        M.set
        (
            cellI, 
            new scalarRectangularMatrix
            (
                cells[cellI].size(),
                nCoeffs,
                0.0
            )
        );
    }

    
    labelList curRow(mesh().nCells(), 0);
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];

        M[own][curRow[own]][0] = d.x();
        M[own][curRow[own]][1] = d.y();
        M[own][curRow[own]][2] = d.z();
        M[own][curRow[own]][3] = d.x()*d.y();
        M[own][curRow[own]][4] = d.x()*d.z();
        M[own][curRow[own]][5] = d.y()*d.z();

        curRow[own]++;

        M[nei][curRow[nei]][0] = -d.x();
        M[nei][curRow[nei]][1] = -d.y();
        M[nei][curRow[nei]][2] = -d.z();
        M[nei][curRow[nei]][3] = d.x()*d.y();
        M[nei][curRow[nei]][4] = d.x()*d.z();
        M[nei][curRow[nei]][5] = d.y()*d.z();

        curRow[nei]++;

//         M[own][curRow[own]][0] = d.x();
//         M[own][curRow[own]][1] = d.y();
//         M[own][curRow[own]][2] = d.x()*d.y();
//         M[own][curRow[own]][3] = d.z();
//         M[own][curRow[own]][4] = d.x()*d.z();
//         M[own][curRow[own]][5] = d.y()*d.z();

//         curRow[own]++;

//         M[nei][curRow[nei]][0] = -d.x();
//         M[nei][curRow[nei]][1] = -d.y();
//         M[nei][curRow[nei]][2] = d.x()*d.y();
//         M[nei][curRow[nei]][3] = -d.z();
//         M[nei][curRow[nei]][4] = d.x()*d.z();
//         M[nei][curRow[nei]][5] = d.y()*d.z();

//         curRow[nei]++;
    }

    forAll(lsP.boundaryField(), patchi)
    {
        const fvPatch& p = mesh().boundary()[patchi];

        const unallocLabelList& faceCells = p.patch().faceCells();

        vectorField pd = p.delta();

        if (isA<emptyFvPatch>(p))
        {
            const polyPatch& pp = p.patch();

            vectorField cc(pp.size(), vector::zero);
            const unallocLabelList& faceCells = pp.faceCells();
            forAll (faceCells, faceI)
            {
                cc[faceI] = C[faceCells[faceI]];
            }

            pd = pp.faceCentres() - cc;
        }

        forAll(pd, patchFacei)
        {
            const vector& d = pd[patchFacei];

            label own = faceCells[patchFacei];

            M[own][curRow[own]][0] = d.x();
            M[own][curRow[own]][1] = d.y();
            M[own][curRow[own]][2] = d.z();
            M[own][curRow[own]][3] = d.x()*d.y();
            M[own][curRow[own]][4] = d.x()*d.z();
            M[own][curRow[own]][5] = d.y()*d.z();

//             M[own][curRow[own]][0] = d.x();
//             M[own][curRow[own]][1] = d.y();
//             M[own][curRow[own]][2] = d.x()*d.y();
//             M[own][curRow[own]][3] = d.z();
//             M[own][curRow[own]][4] = d.x()*d.z();
//             M[own][curRow[own]][5] = d.y()*d.z();

            curRow[own]++;
        }

//         if (p.coupled())
//         {
//             forAll(pd, patchFacei)
//             {
//                 const vector& d = pd[patchFacei];

//                 dd[faceCells[patchFacei]] += (1.0/magSqr(d))*sqr(d);
//             }
//         }
//         else
//         {
//             forAll(pd, patchFacei)
//             {
//                 const vector& d = pd[patchFacei];

//                 dd[faceCells[patchFacei]] += (1.0/magSqr(d))*sqr(d);
//             }
//         }
    }


    // Invert the dd tensor
//     symmTensorField invDd = inv(dd);
    // Fix: householder inverse.  HJ, 3/Nov/2009
//     symmTensorField invDd = hinv(dd);

    forAll(invLsMatrices_, cellI)
    {
        scalarRectangularMatrix& curMatrix = invLsMatrices_[cellI];

        scalarSquareMatrix lsM(nCoeffs, 0.0);

        if (M[cellI].n() == nCoeffs)
        {
            for (label i=0; i<nCoeffs; i++)
            {
                for (label j=0; j<nCoeffs; j++)
                {
                    lsM[i][j] = M[cellI][i][j];
                }
            }            

            scalarSquareMatrix invLsM = lsM.LUinvert();

            for (label i=0; i<nCoeffs; i++)
            {
                for (label j=0; j<nCoeffs; j++)
                {
                    curMatrix[i][j] = invLsM[i][j];
                }
            }

            Info << cellI << ", " << invLsM << endl;            
        }
        else
        {
            for (label i=0; i<nCoeffs; i++)
            {
                for (label j=0; j<nCoeffs; j++)
                {
                    for (label k=0; k<M[cellI].n(); k++)
                    {
                        lsM[i][j] += M[cellI][k][i]*M[cellI][k][j];
                    }
                }
            }

            scalarSquareMatrix invLsM = lsM.LUinvert();

            for (label i=0; i<nCoeffs; i++)
            {
                for (label j=0; j<M[cellI].n(); j++)
                {
                    for (label k=0; k<nCoeffs; k++)
                    {
                        curMatrix[i][j] += invLsM[i][k]*M[cellI][j][k];
                    }
                }
            }

            Info << cellI << ", " << curMatrix << endl;
        }

//         Info << cellI << ", " << M[cellI] << endl;
//         Info << cellI << ", " << lsM << endl;

//         scalarSquareMatrix invTmpM = tmpM.LUinvert();

//         Info << cellI << ", " << invTmpM << endl;
    }

    // Revisit all faces and calculate the lsP and lsN vectors
    curRow = 0;
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

//         vector d = C[nei] - C[own];
//         scalar magSfByMagSqrd = 1.0/magSqr(d);
//         lsP[facei] = magSfByMagSqrd*(invDd[own] & d);
//         lsN[facei] = -magSfByMagSqrd*(invDd[nei] & d);

        lsP[facei] = 
            vector
            (
                invLsMatrices_[own][0][curRow[own]],
                invLsMatrices_[own][1][curRow[own]],
                invLsMatrices_[own][2][curRow[own]]
            );

        curRow[own]++;

        lsN[facei] = 
            vector
            (
                invLsMatrices_[nei][0][curRow[nei]],
                invLsMatrices_[nei][1][curRow[nei]],
                invLsMatrices_[nei][2][curRow[nei]]
            );

        curRow[nei]++;
    }

    forAll(lsP.boundaryField(), patchi)
    {
        const fvPatch& p = mesh().boundary()[patchi];
        const unallocLabelList& faceCells = p.patch().faceCells();

        fvsPatchVectorField& patchLsP = lsP.boundaryField()[patchi];

        
        if (!isA<emptyFvPatch>(p))
        {
            forAll(patchLsP, patchFacei)
            {
                label own = faceCells[patchFacei];

                patchLsP[patchFacei] = 
                    vector
                    (
                        invLsMatrices_[own][0][curRow[own]],
                        invLsMatrices_[own][1][curRow[own]],
                        invLsMatrices_[own][2][curRow[own]]
                    );
                
                curRow[own]++;
            }
        }
        else
        {
            forAll(faceCells, patchFacei)
            {
                label own = faceCells[patchFacei];

                emptyPatchPVectors[patchFacei] = 
                    vector
                    (
                        invLsMatrices_[own][0][curRow[own]],
                        invLsMatrices_[own][1][curRow[own]],
                        invLsMatrices_[own][2][curRow[own]]
                    );

                curRow[own]++;
            }
        }

//         if (p.coupled())
//         {
//             forAll(pd, patchFacei)
//             {
//                 const vector& d = pd[patchFacei];

//                 patchLsP[patchFacei] =
//                     (1.0/magSqr(d))
//                    *(invDd[faceCells[patchFacei]] & d);
//             }
//         }
//         else
//         {
//             forAll(pd, patchFacei)
//             {
//                 const vector& d = pd[patchFacei];

//                 patchLsP[patchFacei] =
//                     (1.0/magSqr(d))
//                    *(invDd[faceCells[patchFacei]] & d);
//             }
//         }
    }

//     if (debug)
    {
        Info<< "leastSquaresVectorsBiLinear::makeLeastSquaresVectorsBiLinear() :"
            << "Finished constructing least square gradient vectors"
            << endl;
    }
}


const Foam::PtrList<Foam::scalarRectangularMatrix>& 
Foam::leastSquaresVectorsBiLinear::invLsMatrices() const
{
    label size = invLsMatrices_.size();

    reduce(size, maxOp<label>());

    if (size == 0)
    {
        makeLeastSquaresVectors();
    }

    return invLsMatrices_;
}


const Foam::surfaceVectorField& Foam::leastSquaresVectorsBiLinear::pVectors() const
{
    if (!pVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *pVectorsPtr_;
}


const Foam::surfaceVectorField& Foam::leastSquaresVectorsBiLinear::nVectors() const
{
    if (!nVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *nVectorsPtr_;
}


const Foam::vectorField& Foam::leastSquaresVectorsBiLinear::emptyPatchPVectors() const
{
    if (!emptyPatchPVectorsPtr_)
    {
        makeLeastSquaresVectors();
    }

    return *emptyPatchPVectorsPtr_;
}


bool Foam::leastSquaresVectorsBiLinear::movePoints() const
{
    if (debug)
    {
        InfoIn("bool leastSquaresVectorsBiLinear::movePoints() const")
            << "Clearing least square data" << endl;
    }

    invLsMatrices_.clear();
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
    deleteDemandDrivenData(emptyPatchPVectorsPtr_);

    return true;
}

bool Foam::leastSquaresVectorsBiLinear::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn("bool leastSquaresVectorsBiLinear::updateMesh(const mapPolyMesh&) const")
            << "Clearing least square data" << endl;
    }

    invLsMatrices_.clear();
    deleteDemandDrivenData(pVectorsPtr_);
    deleteDemandDrivenData(nVectorsPtr_);
    deleteDemandDrivenData(emptyPatchPVectorsPtr_);

    return true;
}

// ************************************************************************* //

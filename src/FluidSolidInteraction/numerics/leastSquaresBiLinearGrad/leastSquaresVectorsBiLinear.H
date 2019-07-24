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

Class
    Foam::leastSquaresVectorsBiLinear

Description
    Least-squares gradient scheme vectors

SourceFiles
    leastSquaresVectorsBiLinear.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresVectorsBiLinear_H
#define leastSquaresVectorsBiLinear_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "labelPair.H"
#include "scalarMatrices.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                     Class leastSquaresVectorsBiLinear Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresVectorsBiLinear
:
    public MeshObject<fvMesh, leastSquaresVectorsBiLinear>
{
    // Private data

        //- Inverse of least-squares normal equation matrix
        mutable PtrList<scalarRectangularMatrix> invLsMatrices_;

        //- Least-squares gradient vectors
        mutable surfaceVectorField* pVectorsPtr_;
        mutable surfaceVectorField* nVectorsPtr_;

        //- Least-squares gradient vectors for empty patch
        mutable vectorField* emptyPatchPVectorsPtr_;


    // Private member functions

        //- Construct Least-squares gradient vectors
        void makeLeastSquaresVectors() const;


public:

    // Declare name of the class and its debug switch
    TypeName("leastSquaresVectorsBiLinear");


    // Constructors

        //- Construct given an fvMesh
        explicit leastSquaresVectorsBiLinear(const fvMesh&);


    // Destructor

        virtual ~leastSquaresVectorsBiLinear();


    // Member functions

        //- Get interpolation matrices
        const PtrList<scalarRectangularMatrix>& invLsMatrices() const;

        //- Return reference to owner least square vectors
        const surfaceVectorField& pVectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& nVectors() const;

        //- Return reference to owner least square vectors for empty patch
        const vectorField& emptyPatchPVectors() const;

        //- Update after mesh motion:
        //  Delete the least square vectors when the mesh moves
        virtual bool movePoints() const;

        //- Update after topo change:
        //  Delete the least square vectors when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

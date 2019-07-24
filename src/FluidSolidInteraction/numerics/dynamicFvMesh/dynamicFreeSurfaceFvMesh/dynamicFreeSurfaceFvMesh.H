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
    dynamicFreeSurfaceFvMesh

Description
    Automatic motion of the mesh around a free-surface.

SourceFiles
    dynamicFreeSurfaceFvMesh.C

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef dynamicFreeSurfaceFvMesh_H
#define dynamicFreeSurfaceFvMesh_H

#include "dynamicFvMesh.H"
#include "dimensionedTypes.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class dynamicFreeSurfaceFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicFreeSurfaceFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary dynamicMeshCoeffs_;

        //- Moving mesh solver
        autoPtr<motionSolver> motionPtr_;

        //- Name of the patch which encloses the body
        word freeSurfacePatchName_;

        //- Body patch ID
        label freeSurfacePatchID_;

        //- Finite area mesh which consists
        //  free-surface faces
        faMesh aMesh_;

        //- Free surface patches which do not move
        wordList fixedFreeSurfacePatches_;

        //- Points which are attached to the free-surface faces
        //  and which defines the free-surface shape
        mutable vectorIOField* controlPointsPtr_;

        //- Field which additionally determines 
        //  the motion of free-surface points
        mutable labelList* motionPointsMaskPtr_;

        //- Displacement direction of free-surface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of free-surface control points
        mutable vectorField* facesDisplacementDirPtr_;

        // Control points displacement
        scalarField deltaH_;

        // Total control points displacement
        scalarField H_;

        //- ROM fields
        DynamicList<scalarField> V_;
        DynamicList<scalarField> W_;

    // Private Member Functions

        // Make demand-driven data

            void makeControlPoints();
            void makeMotionPointsMask();
            void makeDirections();

        //- Clear all demand-driven data
        void clearOut();

        //- Move control ponts for deltaH and calculate
        //  free surface points displacement
        //  for given new control points position
        tmp<vectorField> pointDisplacement(const scalarField& deltaH);

        // Calc least sqare plane point and normal
        tmp<vectorField> lsPlanePointAndNormal
        (
            const vectorField& points,
            const vector& origin,
            const vector& axis
        ) const;

        //- Disallow default bitwise copy construct
        dynamicFreeSurfaceFvMesh(const dynamicFreeSurfaceFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicFreeSurfaceFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicFreeSurfaceFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicFreeSurfaceFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicFreeSurfaceFvMesh();


    // Member Functions

        //- Return reference to finite area mesh
        faMesh& aMesh();

        //- Return reference to finite area mesh
        fvMesh& mesh()
        {
            return *this;
        }

        //- Return control points
        vectorField& controlPoints();

        //- Return reference to motion points mask field
        labelList& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

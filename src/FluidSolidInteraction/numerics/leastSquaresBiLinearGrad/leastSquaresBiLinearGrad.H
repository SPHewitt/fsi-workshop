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
    Foam::fv::leastSquaresBiLinearGrad

Description
    Second-order gradient scheme using least-squares.

SourceFiles
    leastSquaresBiLinearGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresBiLinearGrad_H
#define leastSquaresBiLinearGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class leastSquaresBiLinearGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresBiLinearGrad
:
    public fv::gradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        leastSquaresBiLinearGrad(const leastSquaresBiLinearGrad&);

        //- Disallow default bitwise assignment
        void operator=(const leastSquaresBiLinearGrad&);


public:

    //- Runtime type information
    TypeName("leastSquaresBiLinear");


    // Constructors

        //- Construct from mesh
        leastSquaresBiLinearGrad(const fvMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresBiLinearGrad(const fvMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


    // Member Functions

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > grad
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "leastSquaresBiLinearGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

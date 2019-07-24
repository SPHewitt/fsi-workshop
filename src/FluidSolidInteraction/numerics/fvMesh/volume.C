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

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ggiFvPatch.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<scalarField> volume(const fvMesh& mesh)
{
    tmp<scalarField> tV
    (
        new scalarField(mesh.nCells(), 0)
    );

    const vectorField& points = mesh.points();

    const faceList& faces = mesh.faces();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

//     Info << mesh.nInternalFaces() << ", " << owner.size() << ", "
//         << neighbour.size() << endl;

    scalarField& V = tV();

    forAll(owner, faceI)
    {
        const face& curFace = faces[faceI];

        // If the face is a triangle, do a direct calculation
        if (curFace.size() == 3)
        {
            scalar SR = (curFace.normal(points)&curFace.centre(points));

            V[owner[faceI]] += SR;
            V[neighbour[faceI]] -= SR;
        }
        else
        {
            label nPoints = curFace.size();

            point centrePoint = point::zero;

            for (register label pI=0; pI<nPoints; pI++)
            {
                centrePoint += points[curFace[pI]];
            }

            centrePoint /= nPoints;

            for (register label pI=0; pI<nPoints; pI++)
            {
                // Calculate triangle area
                vector St =
                    (
                        (points[curFace[pI]] - centrePoint)
                      ^ (
                            points[curFace[(pI + 1) % nPoints]] 
                          - centrePoint
                        )
                    );
                St /= 2.0;

                // Calculate triangle centre
                vector Ct = 
                    (
                        centrePoint
                      + points[curFace[pI]]
                      + points[curFace[(pI + 1) % nPoints]]
                    )/3;

                V[owner[faceI]] += (St&Ct);
                V[neighbour[faceI]] -= (St&Ct);
            }
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundaryMesh()[patchI].faceCells();

        forAll(mesh.boundaryMesh()[patchI], faceI)
        {
            label globalFaceID = 
                mesh.boundaryMesh()[patchI].start() + faceI;

            const face& curFace = faces[globalFaceID];

            if (isA<wedgeFvPatch>(mesh.boundary()[patchI]))
            {
                V[pFaceCells[faceI]] +=
                    (curFace.normal(points)&curFace.centre(points));          
            }
            else if (curFace.size() == 3)
            {
                V[pFaceCells[faceI]] +=
                    (curFace.normal(points)&curFace.centre(points));
            }
            else
            {
                label nPoints = curFace.size();

                point centrePoint = point::zero;

                for (register label pI=0; pI<nPoints; pI++)
                {
                    centrePoint += points[curFace[pI]];
                }

                centrePoint /= nPoints;

                for (register label pI=0; pI<nPoints; pI++)
                {
                    // Calculate triangle area
                    vector St =
                        (
                            (points[curFace[pI]] - centrePoint)
                          ^ (
                                points[curFace[(pI + 1) % nPoints]] 
                              - centrePoint
                            )
                        );
                    St /= 2.0;

                    // Calculate triangle centre
                    vector Ct = 
                        (
                            centrePoint
                          + points[curFace[pI]]
                          + points[curFace[(pI + 1) % nPoints]]
                        )/3;

                    V[pFaceCells[faceI]] += (St&Ct);
                }
            }
        }
    }

    V /= 3;

    return tV;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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

Description

\*---------------------------------------------------------------------------*/

#include "EulerD2dt2Scheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    makeFvD2dt2Scheme(EulerD2dt2Scheme)


template<>
tmp<fvMatrix<vector> >
EulerD2dt2Scheme<vector>::fvmD2dt2
(
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
//     Info << "EulerD2dt2Scheme<vector>::fvmD2dt2" << endl;

    tmp<fvMatrix<vector> > tfvm
    (
        new fvMatrix<vector>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<vector>& fvm = tfvm();

    scalar deltaT = mesh().time().deltaT().value();
    scalar deltaT0 = mesh().time().deltaT0().value();

    scalar coefft   = 1.0;
    scalar coefft00 = deltaT/deltaT0;
    scalar coefft0  = coefft + coefft00;

    scalar rDeltaT2 = 1.0/sqr(deltaT);

    if (mesh().moving())
    {
        scalarField V0oV = mesh().V0()/mesh().V();

        fvm.diag() = (coefft*rDeltaT2)*mesh().V();

        fvm.source() = rDeltaT2*mesh().V()*
        (
            (1.0 + coefft00*V0oV)*vf.oldTime().internalField()
          - (coefft00*V0oV)*vf.oldTime().oldTime().internalField()
        );

        if 
        (
            mesh().objectRegistry::found("grad(" + vf.name() + ")")
         && mesh().objectRegistry::found("meshU")
        )
        {
            const volTensorField& gradVf = 
                mesh().objectRegistry::lookupObject<volTensorField>
                (
                    "grad(" + vf.name() + ")"
                );

            const volVectorField& meshU = 
                mesh().objectRegistry::lookupObject<volVectorField>
                (
                    "meshU"
                );

            fvm.source() += 
                rDeltaT2*mesh().V().field()*deltaT
               *(
                    (meshU.internalField() & gradVf.oldTime().internalField())
                  - (
                        meshU.oldTime().internalField()
                      & (
                            gradVf.oldTime().oldTime().internalField()*V0oV
                        )
                    )
                );
        }
    }
    else
    {
        fvm.diag() = (coefft*rDeltaT2)*mesh().V();

        fvm.source() = rDeltaT2*mesh().V()*
        (
            coefft0*vf.oldTime().internalField()
          - coefft00*vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}

}
}

// ************************************************************************* //

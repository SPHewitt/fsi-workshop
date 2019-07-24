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

#include "backwardVectorD2dt2Scheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<>
scalar backwardD2dt2Scheme<vector>::deltaT0_
(
    GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
    // Bug fix, Zeljko Tukovic: solver with outer iterations over a time-step
    // HJ, 12/Feb/2010
//     if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())

    if (mesh().found("ddt(" + vf.name() + ")"))
    {
        if 
        (
            vf.oldTime().timeIndex()
         == vf.oldTime().oldTime().timeIndex()
        )
        {
            return GREAT;
        }
        else
        {
            return deltaT0_();
        }
    }
    else
    {
        if 
        (
            vf.oldTime().oldTime().timeIndex()
         == vf.oldTime().oldTime().oldTime().timeIndex()
        )
        {
            return GREAT;
        }
        else
        {
            return deltaT0_();
        }
    }
}


template<>
tmp<fvMatrix<vector> >
backwardD2dt2Scheme<vector>::fvmD2dt2
(
    GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<vector> > tfvm
    (
        new fvMatrix<vector>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<vector>& fvm = tfvm();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        notImplemented
        (
            type()
          + "::fvmD2dt2(GeometricField<vector, fvPatchField, volMesh>& vf)"
        );
    }
    else
    {
        fvm = coefft*dimensionedScalar("rDeltaT", dimless/dimTime, rDeltaT)
           *backwardDdtScheme<vector>(mesh()).fvmDdt(vf);

        if (mesh().found("ddt(" + vf.name() + ")"))
        {
//             Info << "Using stored ddt field" << endl;

            const volVectorField& ddtVf = 
                mesh().lookupObject<volVectorField>("ddt(" + vf.name() + ")");

            fvm.source() += rDeltaT*mesh().V()*
            (
                coefft0*ddtVf.oldTime().internalField()
              - coefft00*ddtVf.oldTime().oldTime().internalField()
            );
        }
        else
        {
//             Info << "Calculating ddt field" << endl;
//             Info << "ddt(" + vf.name() + ")" << endl;

            fvm.source() += rDeltaT*mesh().V()*
            (
                coefft0
               *backwardDdtScheme<vector>(mesh()).fvcDdt(vf.oldTime())
                ().internalField()
              - coefft00
               *backwardDdtScheme<vector>(mesh()).fvcDdt(vf.oldTime().oldTime())
                ().internalField()
            );
        }
    }

    return tfvm;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

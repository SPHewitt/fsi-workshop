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

#include "EulerD2dt2Scheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerD2dt2Scheme<Type>::fvcD2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    IOobject d2dt2IOobject
    (
        "d2dt2("+vf.name()+')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT2 = 1.0/sqr(mesh().time().deltaT());

    scalar deltaT = mesh().time().deltaT().value();
    scalar deltaT0 = mesh().time().deltaT0().value();

    scalar coefft   = 1.0;
    scalar coefft00 = deltaT/deltaT0;
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        scalarField V0oV = mesh().V0()/mesh().V();

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*vf.dimensions(),
                rDeltaT2.value()*
                (
                    coefft*vf.internalField()
                  - (1.0 + coefft00*V0oV)*vf.oldTime().internalField()
                  + coefft00*V0oV*vf.oldTime().oldTime().internalField()
                ),
                rDeltaT2.value()*
                (
                    coefft*vf.boundaryField()
                  - coefft0*vf.oldTime().boundaryField()
                  + coefft00*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                rDeltaT2*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT2 =
        4.0/sqr(mesh().time().deltaT() + mesh().time().deltaT0());

    IOobject d2dt2IOobject
    (
        "d2dt2("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    scalar deltaT = mesh().time().deltaT().value();
    scalar deltaT0 = mesh().time().deltaT0().value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2.value();
        scalar quarterRdeltaT2 = 0.25*rDeltaT2.value();

        scalarField VV0rhoRho0 =
            (mesh().V() + mesh().V0())
           *(rho.internalField() + rho.oldTime().internalField());

        scalarField V0V00rho0Rho00 =
            (mesh().V0() + mesh().V00())
           *(
               rho.oldTime().internalField()
             + rho.oldTime().oldTime().internalField()
            );

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                mesh(),
                rDeltaT2.dimensions()*rho.dimensions()*vf.dimensions(),
                quarterRdeltaT2*
                (
                    coefft*VV0rhoRho0*vf.internalField()

                  - (coefft*VV0rhoRho0 + coefft00*V0V00rho0Rho00)
                   *vf.oldTime().internalField()

                  + (coefft00*V0V00rho0Rho00)
                   *vf.oldTime().oldTime().internalField()
                )/mesh().V(),
                halfRdeltaT2*
                (
                    coefft
                   *(rho.boundaryField() + rho.oldTime().boundaryField())
                   *vf.boundaryField()

                  - (
                        coefft
                       *(
                           rho.boundaryField()
                         + rho.oldTime().boundaryField()
                        )
                      + coefft00
                       *(
                           rho.oldTime().boundaryField()
                         + rho.oldTime().oldTime().boundaryField()
                        )
                    )*vf.oldTime().boundaryField()

                  + coefft00
                   *(
                       rho.oldTime().boundaryField()
                     + rho.oldTime().oldTime().boundaryField()
                    )*vf.oldTime().oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        dimensionedScalar halfRdeltaT2 = 0.5*rDeltaT2;

        volScalarField rhoRho0 = rho + rho.oldTime();
        volScalarField rho0Rho00 = rho.oldTime() +rho.oldTime().oldTime();

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                d2dt2IOobject,
                halfRdeltaT2*
                (
                    coefft*rhoRho0*vf
                  - (coefft*rhoRho0 + coefft00*rho0Rho00)*vf.oldTime()
                  + coefft00*rho0Rho00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
EulerD2dt2Scheme<Type>::fvmD2dt2
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

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


template<class Type>
tmp<fvMatrix<Type> >
EulerD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol
            /dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar deltaT = mesh().time().deltaT().value();
    scalar deltaT0 = mesh().time().deltaT0().value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2;

        scalarField VV0 = mesh().V() + mesh().V0();

        scalarField V0V00 = mesh().V0() + mesh().V00();

        fvm.diag() = rho.value()*(coefft*halfRdeltaT2)*VV0;

        fvm.source() = halfRdeltaT2*rho.value()*
        (
            (coefft*VV0 + coefft00*V0V00)
           *vf.oldTime().internalField()

          - (coefft00*V0V00)*vf.oldTime().oldTime().internalField()
        );
    }
    else
    {
        fvm.diag() = (coefft*rDeltaT2)*mesh().V()*rho.value();

        fvm.source() = rDeltaT2*mesh().V()*rho.value()*
        (
            (coefft + coefft00)*vf.oldTime().internalField()
          - coefft00*vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol
            /dimTime/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar deltaT = mesh().time().deltaT().value();
    scalar deltaT0 = mesh().time().deltaT0().value();

    scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
    scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);

    scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);

    if (mesh().moving())
    {
        scalar quarterRdeltaT2 = 0.25*rDeltaT2;

        scalarField VV0rhoRho0 =
            (mesh().V() + mesh().V0())
           *(rho.internalField() + rho.oldTime().internalField());

        scalarField V0V00rho0Rho00 =
            (mesh().V0() + mesh().V00())
           *(
               rho.oldTime().internalField()
             + rho.oldTime().oldTime().internalField()
            );

        fvm.diag() = (coefft*quarterRdeltaT2)*VV0rhoRho0;

        fvm.source() = quarterRdeltaT2*
        (
            (coefft*VV0rhoRho0 + coefft00*V0V00rho0Rho00)
           *vf.oldTime().internalField()

          - (coefft00*V0V00rho0Rho00)
           *vf.oldTime().oldTime().internalField()
        );
    }
    else
    {
        scalar halfRdeltaT2 = 0.5*rDeltaT2;

        scalarField rhoRho0 =
            (rho.internalField() + rho.oldTime().internalField());

        scalarField rho0Rho00 =
        (
            rho.oldTime().internalField()
          + rho.oldTime().oldTime().internalField()
        );

        fvm.diag() = (coefft*halfRdeltaT2)*mesh().V()*rhoRho0;

        fvm.source() = halfRdeltaT2*mesh().V()*
        (
            (coefft*rhoRho0 + coefft00*rho0Rho00)
           *vf.oldTime().internalField()

          - (coefft00*rho0Rho00)
           *vf.oldTime().oldTime().internalField()
        );
    }

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

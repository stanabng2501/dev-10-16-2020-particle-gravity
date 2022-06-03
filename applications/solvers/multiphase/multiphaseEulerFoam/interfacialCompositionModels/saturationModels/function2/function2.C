/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "function2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(function2, 0);
    addToRunTimeSelectionTable(saturationModel, function2, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::function2::function2
(
    const dictionary& dict,
    const phasePair& pair
)
:
    saturationModel(pair),
    function_
    (
        Function1<scalar>::New("function", dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::function2::~function2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::function2::Tsat
(
    const volScalarField& p
) const
{
    NotImplemented;
    return volScalarField::null();
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::function2::pSatPrime
(
    const volScalarField& T
) const
{
    NotImplemented;
    return volScalarField::null();
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::function2::lnPSat
(
    const volScalarField& T
) const
{
    NotImplemented;
    return volScalarField::null();
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::function2::pSat
(
    const volScalarField& T
) const
{
    tmp<volScalarField> tPsat
    (
        volScalarField::New
        (
            "Psat",
            T.mesh(),
            dimensionedScalar(dimPressure, 0)
        )
    );

    volScalarField& Psat = tPsat.ref();

    forAll(Psat, celli)
    {
        Psat[celli] = function_->value(T[celli]);
    }

    volScalarField::Boundary& PsatBf = Psat.boundaryFieldRef();

    forAll(Psat.boundaryField(), patchi)
    {
        scalarField& Psatp = PsatBf[patchi];
        const scalarField& pp = T.boundaryField()[patchi];

        forAll(Psatp, facei)
        {
            Psatp[facei] = function_->value(pp[facei]);

        }
    }

    return tPsat;
}


// ************************************************************************* //

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

#include "function3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationDensityModels
{
    defineTypeNameAndDebug(function3, 0);
    addToRunTimeSelectionTable(saturationDensityModel, function3, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationDensityModels::function3::function3
(
    const dictionary& dict,
    const phasePair& pair
)
:
    saturationDensityModel(pair),
    function_
    (
        Function1<scalar>::New("function", dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationDensityModels::function3::~function3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

Foam::tmp<Foam::volScalarField>
Foam::saturationDensityModels::function3::rhoSat
(
    const volScalarField& T
) const
{
    tmp<volScalarField> trhoSat
    (
        volScalarField::New
        (
            "rhoSat",
            T.mesh(),
            dimensionedScalar(dimDensity, 0)
        )
    );

    volScalarField& rhoSat = trhoSat.ref();

    forAll(rhoSat, celli)
    {
        rhoSat[celli] = function_->value(T[celli]);
    }

    volScalarField::Boundary& rhoSatBf = rhoSat.boundaryFieldRef();

    forAll(rhoSat.boundaryField(), patchi)
    {
        scalarField& rhoSatp = rhoSatBf[patchi];
        const scalarField& tt = T.boundaryField()[patchi];

        forAll(rhoSatp, facei)
        {
            rhoSatp[facei] = function_->value(tt[facei]);

        }
    }

    return trhoSat;
}


// ************************************************************************* //

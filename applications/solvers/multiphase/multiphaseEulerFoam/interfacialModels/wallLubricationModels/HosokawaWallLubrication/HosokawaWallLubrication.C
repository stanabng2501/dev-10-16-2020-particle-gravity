/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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

#include "HosokawaWallLubrication.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(HosokawaWallLubrication, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        HosokawaWallLubrication,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::HosokawaWallLubrication::HosokawaWallLubrication
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallLubricationModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::HosokawaWallLubrication::~HosokawaWallLubrication()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::wallLubricationModels::HosokawaWallLubrication::Fi() const
{
    volVectorField Ur(pair_.Ur());

    const volVectorField& n(nWall());
    const volScalarField& y(yWall());

    volScalarField Eo(pair_.Eo());

    return zeroGradWalls
    (
        pair_.continuous().rho()
        *0.0217*Eo
        *(2/pair_.dispersed().d())
         *magSqr(Ur)
         *sqr(pair_.dispersed().d()/(2*y))
         *n       
    );
}


// ************************************************************************* //

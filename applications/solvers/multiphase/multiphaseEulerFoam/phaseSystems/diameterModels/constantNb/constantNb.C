/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "constantNb.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(constantNb, 0);
    addToRunTimeSelectionTable(diameterModel, constantNb, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::constantNb::calcD() const
{
    return d_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::constantNb::constantNb
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    spherical(diameterProperties, phase),
    dmin_("dmin", dimLength, diameterProperties),
    nB_("nB", dimensionSet(0, -3, 0, 0, 0), diameterProperties),
    d_(dRef())
{
    d_ = dmin_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::constantNb::~constantNb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::constantNb::correct()
{

    d_ = max(pow((6*phase()/(nB_* constant::mathematical::pi)), 1.0/3.0), dmin_)  ;
}


bool Foam::diameterModels::constantNb::read(const dictionary& phaseProperties)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("dmin") >> dmin_;
    diameterProperties().lookup("nB") >> nB_;

    return true;
}


// ************************************************************************* //

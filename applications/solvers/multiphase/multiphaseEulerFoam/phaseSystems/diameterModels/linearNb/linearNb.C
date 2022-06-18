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

#include "linearNb.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(linearNb, 0);
    addToRunTimeSelectionTable(diameterModel, linearNb, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::linearNb::calcD() const
{
    return d_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::linearNb::linearNb
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    spherical(diameterProperties, phase),
    dmin_("dmin", dimLength, diameterProperties),
    nB_("nB", dimensionSet(0, -3, 0, 0, 0), diameterProperties),
    d_(dRef()),
    dmax_(pow(6/(nB_* constant::mathematical::pi) , 1.0/3.0)) 
{
    d_ = dmin_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::linearNb::~linearNb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::linearNb::correct()
{

    d_ = dmin_ + phase()*(dmax_ - dmin_) ;
    
        
    
}


bool Foam::diameterModels::linearNb::read(const dictionary& phaseProperties)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("dmin") >> dmin_;
    diameterProperties().lookup("nB") >> nB_;
    diameterProperties().lookup("dmax_") >> dmax_;

    return true;
}


// ************************************************************************* //

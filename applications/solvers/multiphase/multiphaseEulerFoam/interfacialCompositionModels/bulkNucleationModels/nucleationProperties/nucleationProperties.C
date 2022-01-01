/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "nucleationProperties.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceTensionModel.H" 
#include "fundamentalConstants.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bulkNucleationModels
{
    defineTypeNameAndDebug(nucleationProperties, 0);
    addToRunTimeSelectionTable
    (
        bulkNucleationModel,
        nucleationProperties,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulkNucleationModels::nucleationProperties::nucleationProperties
(
    const dictionary& dict,
    const phasePair& pair
)
:
    bulkNucleationModel(dict, pair) 
{
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulkNucleationModels::nucleationProperties::~nucleationProperties()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
 
 
Foam::tmp<Foam::volScalarField>
Foam::bulkNucleationModels::nucleationProperties::B
 (
   const phaseModel& phase1 ,
   const phaseModel& phase2 ,
   const volScalarField& pSat,
   const word& transfertype 
 ) const
 
 {
      Info << "nucleation Properties model does not calculate frequency of nucleation"<<endl;
      NotImplemented;
      return volScalarField::null();
            
 } 
 


 


// ************************************************************************* //

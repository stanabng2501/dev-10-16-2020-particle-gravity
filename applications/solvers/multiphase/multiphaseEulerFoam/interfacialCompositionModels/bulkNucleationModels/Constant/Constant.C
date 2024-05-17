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

#include "Constant.H"
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
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable
    (
        bulkNucleationModel,
        Constant,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulkNucleationModels::Constant::Constant
(
    const dictionary& dict,
    const phasePair& pair
)
:
    bulkNucleationModel(dict, pair),
    bFac_("B", dimensionSet(0, 0, -1, 0, 0), dict) 
{
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulkNucleationModels::Constant::~Constant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
 
Foam::tmp<Foam::volScalarField>
Foam::bulkNucleationModels::Constant::B
 (
   const phaseModel& phase1 ,
   const phaseModel& phase2 ,
   const volScalarField& pSat,
   const word& transfertype 
 ) const
 
 {
    
      return volScalarField::New
    (
        "bFac_",
        pSat.mesh(),
        bFac_
    );
     
          
 } 

 

 //    Info<< "explkT = " << min(explkT.primitiveField()) << "  explkT   max = " << max(explkT.primitiveField()) <<  " explkT  dimensions = " << explkT.dimensions() <<endl;      

 


// ************************************************************************* //

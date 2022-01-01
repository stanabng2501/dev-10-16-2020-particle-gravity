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

#include "Doring.H"
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
    defineTypeNameAndDebug(Doring, 0);
    addToRunTimeSelectionTable
    (
        bulkNucleationModel,
        Doring,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulkNucleationModels::Doring::Doring
(
    const dictionary& dict,
    const phasePair& pair
)
:
    bulkNucleationModel(dict, pair) 
{
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulkNucleationModels::Doring::~Doring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
 
Foam::tmp<Foam::volScalarField>
Foam::bulkNucleationModels::Doring::B
 (
   const phaseModel& phase1 ,
   const phaseModel& phase2 ,
   const volScalarField& pSat,
   const word& transfertype 
 ) const
 
 {
      
     const rhoThermo& thermo1 = phase1.thermo();
     const rhoThermo& thermo2 = phase2.thermo();
     const volScalarField& rho1(thermo1.rho());
     const volScalarField& rho2(thermo2.rho());
     const volScalarField& p(thermo1.p());
     const volScalarField& T1(thermo1.T());
     const volScalarField sigma
                     (
                        phase1.fluid().lookupSubModel<surfaceTensionModel>
                          (
                             phasePair(phase1, phase2)
                          ).sigma()
                      );
    const volScalarField pg( p+(pSat-p)*(1-rho1/rho2) ) ;
                   
    volScalarField b(p/pg) ; 
    b.max(small);
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number  
    const volScalarField lkT( ((thermo1.ha() - thermo2.ha())*thermo1.W()) / (k*T1*1000*Av )  );
            
    return ( exp(-lkT) *sqrt( 1000*Av*6*sigma/((3-b)* constant::mathematical::pi*thermo1.W())) );
            
 } 

 


 


// ************************************************************************* //

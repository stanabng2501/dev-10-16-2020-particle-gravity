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

#include "Alesksandrov.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(Alesksandrov, 0);
    addToRunTimeSelectionTable(heatTransferModel, Alesksandrov, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::Alesksandrov::Alesksandrov
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict, pair) ,
    saturationModels_
    (
       pair.continuous().mesh().foundObject<saturationModel>(IOobject::groupName("saturationModel", pair.name()))
     ? pair.continuous().mesh().lookupObject<saturationModel>(IOobject::groupName("saturationModel", pair.name()))
     :
        saturationModel::New
        (
            dict.subDict("SaturationTemperature"),
            pair
        )
    )

{
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::Alesksandrov::~Alesksandrov()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::Alesksandrov::K(const scalar residualAlpha) const
{
     volScalarField D = (pair_.continuous().thermo().kappa() /( pair_.continuous().rho()* pair_.continuous().thermo().Cpv()));

     // Peclet number
     volScalarField Pe( pair_.magUr() *pair_.dispersed().d() / D);  
     

     const volScalarField& Tc = pair_.continuous().thermo().T();
     const volScalarField& p = pair_.dispersed().thermo().p();
     const volScalarField& Td = pair_.dispersed().thermo().T();     
 
    volScalarField Tsat( saturationModels_.Tsat(p));   
    volScalarField L(pair_.continuous().thermo().ha(p,Tc) -  pair_.dispersed().thermo().ha(p,Td));
    
    
    // Jakob Number                                                                 
    volScalarField Ja(pair_.continuous().thermo().Cpv() * pair_.continuous().rho()* (Tc-Tsat)  / (pair_.dispersed().rho()*L) ); 

    const scalar pi =constant::mathematical::pi;
    volScalarField Nu(sqrt(12*sqr(Ja/pi) +  Pe/(3*pi))); 
      
    return
        6
       *max(pair_.dispersed(), residualAlpha)
       *pair_.continuous().thermo().kappa()
       *Nu
       /sqr(pair_.dispersed().d());
}


// ************************************************************************* //

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

#include "RanzMarshallLimited.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(RanzMarshallLimited, 0);
    addToRunTimeSelectionTable(heatTransferModel, RanzMarshallLimited, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::RanzMarshallLimited::RanzMarshallLimited
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict, pair),
    dNuc_("dNuc", dimLength, dict),
    alphaPhaseTypeName_(dict.lookup("alphaPhaseTypeName")),
    alphaPhaseType_(alphaPhaseTypeName_ == "continuous"? pair.continuous() : pair.dispersed()),
    ReMin_("ReMin", dimless, dict),
    ReMax_("ReMax", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::RanzMarshallLimited::~RanzMarshallLimited()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::RanzMarshallLimited::K(const scalar residualAlpha) const
{
    volScalarField Re(pair_.magUr()*dNuc_
                       /pair_.continuous().thermo().nu()); 
    
    
 //   Info << "Reynolds Number before limiting for pair = " << pair_.name() 
 //        << ",  min " << min(Re.primitiveField()) 
 //        << ",  max " << max(Re.primitiveField()) 
 //        << endl;
    
    Re =   max( min(Re, ReMax_), ReMin_);  
    
    Info << "RanzMarshall Re Limit   min = " << min(Re.primitiveField())
         << ",    max = " << max(Re.primitiveField())
          << endl; 
                
                             
    volScalarField Pr(pair_.Pr()); 
    volScalarField Nu(2 + 0.6*sqrt(Re)*cbrt(Pr));
    volScalarField alphaNu(alphaPhaseType_ * Nu);

    Info << "Alpha Nusselt number for Ranz Marshall for pair = " << pair_.name() 
         << ",  mean " << average(alphaNu.primitiveField()) 
         << ",  min " << min(alphaNu.primitiveField()) 
         << ",  max " << max(alphaNu.primitiveField()) 
         << endl;
         
         
    return
        6
       *max(alphaPhaseType_, residualAlpha)
       *pair_.continuous().thermo().kappa()
       *Nu
       /sqr(dNuc_);
}


// ************************************************************************* //

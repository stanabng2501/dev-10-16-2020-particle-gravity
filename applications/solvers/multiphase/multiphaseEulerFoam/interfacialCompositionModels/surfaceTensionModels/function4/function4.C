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

#include "function4.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(function4, 0);
    addToRunTimeSelectionTable(surfaceTensionModel, function4, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::function4::function4
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    surfaceTensionModel(dict, pair, registerObject),
    phaseName_(dict.lookup("phaseName")),
    function_
    (
        Function1<scalar>::New("function", dict)
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::function4::~function4()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
 


 
Foam::tmp<Foam::scalarField>
Foam::surfaceTensionModels::function4::sigma
(
    label patchi
) const
{

   // Used in KocamustafaogullariIshiiDepartureFrequency.C
    NotImplemented;
    return scalarField::null();
}



Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::function4::sigma
(
    
) const
{

     const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );
        
     const phaseModel& phase(fluid.phases()[phaseName_]);   
     
     const volScalarField& T =phase.thermo().T();
     
      
    tmp<volScalarField> tSigma
    (
        volScalarField::New
        (
            "sigma",
            T.mesh(),
            dimensionedScalar(dimSigma, 0)
        )
    );

    volScalarField& sigma = tSigma.ref();

    forAll(sigma, celli)
    {
        sigma[celli] = function_->value(T[celli]);
    }

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();

    forAll(sigma.boundaryField(), patchi)
    {
        scalarField& sigmaP = sigmaBf[patchi];
        const scalarField& ss = T.boundaryField()[patchi];

        forAll(sigmaP, facei)
        {
            sigmaP[facei] = function_->value(ss[facei]);

        }
    }

    return tSigma;
}


// ************************************************************************* //

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

Class
    Foam::heatTransferModels::RanzMarshall

Description
    Ranz-Marshall correlation for turbulent heat transfer from the surface of a
    sphere to the surrounding fluid.

SourceFiles
    RanzMarshall.C

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshallLimited_H
#define RanzMarshallLimited_H

#include "heatTransferModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                        Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/

class RanzMarshallLimited
:
    public heatTransferModel
{

    //- Nucleation Sauter Diameter
        dimensionedScalar dNuc_;
        
    // Name of phase
       const word alphaPhaseTypeName_;     
        
    // phase type : dispersed or continuous
       const phaseModel& alphaPhaseType_;
       
             
      //- Min Reynolds number
      const dimensionedScalar ReMin_;
      
      //- Max Reynolds number
      const dimensionedScalar ReMax_;
public:

    //- Runtime type information
    TypeName("RanzMarshallLimited");


    // Constructors

        //- Construct from components
        RanzMarshallLimited
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~RanzMarshallLimited();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::interfaceCompositionModels::Doring

Description
    Volmer Doring Nucleation Model.

SourceFiles
    Doring.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationProperties_H
#define nucleationProperties_H

#include "bulkNucleationModel.H"
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace bulkNucleationModels
{

/*---------------------------------------------------------------------------*\
                          Class saturated Declaration
\*---------------------------------------------------------------------------*/

class nucleationProperties
:
    public bulkNucleationModel
{
protected:

    // Private Data

        //- saturated species name
        

        //- saturated species index
        

        


    // Private Member Functions

        //- Constant of proportionality between partial pressure and mass
        //  fraction
      


public:

    //- Runtime type information
    TypeName("nucleationProperties");

    // Constructors

        //- Construct from components
        nucleationProperties
        (
            const dictionary& dict,
            const phasePair& pair
        );

  
       virtual tmp<volScalarField> B
            (
                const phaseModel& phase1 ,
                 const phaseModel& phase2 ,
                  const volScalarField& pSat,
                const word& transfertype 
               
            ) const;
   

    //- Destructor
    virtual ~nucleationProperties();


    // Member Functions

  

        //- The interface species fraction
      

        //- The interface species fraction derivative w.r.t. temperature
         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bulkNucleationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

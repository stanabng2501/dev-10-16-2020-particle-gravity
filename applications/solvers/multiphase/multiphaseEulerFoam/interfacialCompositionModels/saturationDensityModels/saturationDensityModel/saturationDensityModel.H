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
    Foam::saturationModel

Description

SourceFiles
    saturationDensityModel.C
    saturationDensityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef saturationDensityModel_H
#define saturationDensityModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class saturationModel Declaration
\*---------------------------------------------------------------------------*/

class saturationDensityModel
:
    public IOdictionary
{
public:

    //- Runtime type information
    TypeName("saturationDensityModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationDensityModel,
        dictionary,
        (
            const dictionary& dict, const phasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct null
        saturationDensityModel(const phasePair& pair);

        //- Disallow default bitwise copy construction
        saturationDensityModel(const saturationDensityModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<saturationDensityModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~saturationDensityModel();


    // Member Functions

        //- Saturation density
        virtual tmp<volScalarField> rhoSat
        (
            const volScalarField& T
        ) const = 0;

  


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const saturationDensityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

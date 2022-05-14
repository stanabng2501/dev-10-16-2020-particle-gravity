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
    Foam::liftModels::wallDamped

Description

SourceFiles
    wallDamped.C

\*---------------------------------------------------------------------------*/

#ifndef wallDampedDrag_H
#define wallDampedDrag_H

#include "dragModel.H"
#include "wallDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                         Class wallDamped Declaration
\*---------------------------------------------------------------------------*/

class wallDamped
:
    public dragModel
{
    // Private Data

        //- The drag model to damp
        autoPtr<dragModel> dragModel_;

        //- The wall-damping model
        autoPtr<wallDampingModel> wallDampingModel_;


public:

    //- Runtime type information
    TypeName("wallDamped");


    // Constructors

        //- Construct from a dictionary and a phase pair
        wallDamped
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~wallDamped();


    // Member Functions

 

         //- Return the reference to phase-intensive drag coefficient Ki
        virtual tmp<volScalarField> Ki() const;

        //- Return the drag coefficient K
        virtual tmp<volScalarField> K() const;

        //- Return the drag coefficient Kf
        virtual tmp<surfaceScalarField> Kf() const;
        
        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

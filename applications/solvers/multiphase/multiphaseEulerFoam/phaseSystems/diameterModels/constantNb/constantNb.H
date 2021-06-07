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
    Foam::diameterModels::constantNb

Description
    Isothermal dispersed-phase particle diameter model.

SourceFiles
    constantNb.C

\*---------------------------------------------------------------------------*/

#ifndef constantNb_H
#define constantNb_H

#include "sphericalDiameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                         Class isothermal Declaration
\*---------------------------------------------------------------------------*/

class constantNb
:
    public spherical
{
    // Private Data

        //- Minimum  diameter 
        dimensionedScalar dmin_;

        //- Reference bubble number density
        dimensionedScalar nB_;

        //- Diameter field
        volScalarField& d_;


protected:

    // Protected Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> calcD() const;


public:

    //- Runtime type information
    TypeName("constantNb");


    // Constructors

        //- Construct from components
        constantNb
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~constantNb();


    // Member Functions

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

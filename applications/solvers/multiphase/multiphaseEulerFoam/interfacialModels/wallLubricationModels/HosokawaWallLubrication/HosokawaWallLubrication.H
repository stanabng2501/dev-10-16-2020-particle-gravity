/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
    Foam::wallLubricationModels::HosokawaWallLubrication

Description
    Wall lubrication model of Hosokawa.

    References:
v

    \verbatim
        Hosokawa, S., Tomiyama, A., Misaki, S. and Hamada, T., 2002, January. 
        Lateral migration  of single bubbles due to the presence of wall. 
        In Fluids Engineering Division Summer Meeting (Vol. 36150, pp. 855-860).
    \endverbatim

SourceFiles
    HosokawaWallLubrication.C

\*---------------------------------------------------------------------------*/

#ifndef HosokawaWallLubrication_H
#define HosokawaWallLubrication_H

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                   Class TomiyamaWallLubrication Declaration
\*---------------------------------------------------------------------------*/

class HosokawaWallLubrication
:
    public wallLubricationModel
{
    // Private Data


public:

    //- Runtime type information
    TypeName("Hosokawa");


    // Constructors

        //- Construct from components
        HosokawaWallLubrication
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~HosokawaWallLubrication();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

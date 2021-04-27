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
    Foam::heatTransferModels::Alesksandrov

Description
    Alesksandrov correlation for turbulent heat transfer from the surface of a
    sphere to the surrounding fluid.
    
        Based on:

    \verbatim
        Alexandrov, Y. A., Delone, N. B., Gorbunkov, V. M., & Voronov, G. S. (1967).
         Bubble chambers.
    \endverbatim
    
    Validated in:
 
    \verbatim
        Liao, Y. and Lucas, D., 2015. 3D CFD simulation of flashing flows in a 
        converging-diverging nozzle. Nuclear Engineering and Design, 292, pp.149-163.
    \endverbatim

SourceFiles
    Alesksandrov.C

\*---------------------------------------------------------------------------*/

#ifndef Alesksandrov_H
#define Alesksandrov_H

#include "heatTransferModel.H"
 #include "saturationModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                        Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/

class Alesksandrov
:
    public heatTransferModel
{
 


        //- Saturation model
       const saturationModel& saturationModels_;
public:

    //- Runtime type information
    TypeName("Alesksandrov");


    // Constructors

        //- Construct from components
        Alesksandrov
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Alesksandrov();


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

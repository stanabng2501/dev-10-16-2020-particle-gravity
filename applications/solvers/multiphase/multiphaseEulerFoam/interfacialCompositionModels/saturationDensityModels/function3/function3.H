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

Class
    Foam::saturationModels::function1

Description
    Saturation vapour temperature in terms of
    the vapour pressure (in Pa). The saturation temperature in Kelvins is
    specified as a Foam::Function1 type, to enable use of, e.g. constant,
    polynomial, table values.

    Currently this class only provides \f$T_sat\f$, the inverse function to
    return the vapour pressure for a given temperature are not implemented.

    Examples:

    \verbatim
        type function1;
        function polynomial
        (
            (308.0422       0)
            (0.0015096      1)
            (-1.61589e-8    2)
            (1.114106e-13   3)
            (-4.52216e-19   4)
            (1.05192e-24    5)
            (-1.2953e-30    6)
            (6.5365e-37     7)
        )
    \endverbatim

    \verbatim
        type function1;
        function tableFile;
        functionCoeffs
        {
            file                "filename.csv";
            format              csv;
            nHeaderLine         1;
            refColumn           0;
            componentColumns    (1);
            separator           ",";
            mergeSeparators     no;
            outOfBounds         clamp;
            interpolationScheme linear;
        };
    \endverbatim

SourceFiles
    function3.C

\*---------------------------------------------------------------------------*/

#ifndef function3_H
#define function3_H

#include "saturationDensityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationDensityModels
{

/*---------------------------------------------------------------------------*\
                          Class function1 Declaration
\*---------------------------------------------------------------------------*/

class function3
:
    public saturationDensityModel
{
    // Private Data

        //- Saturation temperature as a function of pressure
        autoPtr<Function1<scalar>> function_;


public:

    //- Runtime type information
    TypeName("function3");

    // Constructors

        //- Construct from a dictionary
        function3(const dictionary& dict, const phasePair& pair);


    //- Destructor
    virtual ~function3();


    // Member Functions

        //- Saturation pressure
        virtual tmp<volScalarField> rhoSat(const volScalarField& T) const;

 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

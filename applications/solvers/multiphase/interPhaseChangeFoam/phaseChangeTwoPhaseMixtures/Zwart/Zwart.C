/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "Zwart.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Zwart, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeTwoPhaseMixture,
        Zwart,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Zwart::Zwart
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    Cc_("Cc", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    Cv_("Cv", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    rNuc_("rNuc", dimless, phaseChangeTwoPhaseMixtureCoeffs_),
    Rb_("Rb", dimLength, phaseChangeTwoPhaseMixtureCoeffs_),

    p0_(pSat().dimensions(), Zero)
{
    correct();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Zwart::pCoeff
    (
        const volScalarField& p
    ) const
{

    return
        (3*rho2())*sqrt(2/(3*rho1()))
        /(Rb_*sqrt(mag(p - pSat()) + 0.01*pSat()));
}

    Foam::Pair<Foam::tmp<Foam::volScalarField> >
    Foam::phaseChangeTwoPhaseMixtures::Zwart::mDotAlphal() const
    {
        const volScalarField& p = alpha1().db().lookupObject<volScalarField>("p");
        volScalarField limitedAlpha1(min(max(alpha1(), scalar(0)), scalar(1)));

        volScalarField pCoeff(this->pCoeff(p));

        return Pair<tmp<volScalarField> >
        (
            Cc_*pCoeff*max(p - pSat(), p0_),

            Cv_*rNuc_*pCoeff*min(p - pSat(), p0_) 
        
        );
    }

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::Zwart::mDotP() const
{
    const volScalarField& p = alpha1().db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1(), scalar(0)), scalar(1)));

    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat())*pCoeff,

        (-Cv_)*rNuc_*limitedAlpha1*neg(p - pSat())*pCoeff
    );
}


void Foam::phaseChangeTwoPhaseMixtures::Zwart::correct()
{}


bool Foam::phaseChangeTwoPhaseMixtures::Zwart::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");
	
	phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("rNuc") >> rNuc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Rb") >> Rb_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

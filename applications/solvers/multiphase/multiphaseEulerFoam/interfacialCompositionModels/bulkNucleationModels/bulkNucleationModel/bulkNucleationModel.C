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

\*---------------------------------------------------------------------------*/

#include "bulkNucleationModel.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "phaseSystem.H"
#include "rhoReactionThermo.H"
#include "surfaceTensionModel.H" 
#include "fundamentalConstants.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bulkNucleationModel, 0);
    defineRunTimeSelectionTable(bulkNucleationModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bulkNucleationModel::bulkNucleationModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair),
    WcrkTNmin_(dict.lookup<scalar>("WcrkTNmin")),
    WcrkTDelta_(dict.lookup<scalar>("WcrkTDelta")),    
    WcrFrac_(dict.lookup<scalar>("WcrFrac")),
    D32_(dict.lookupOrDefault<scalar>("D32",1e-4)),   
    pSat_("pSat", dimPressure, dict),  
    rhoVSat_("rhoVSat", dimDensity, dict),     
    rhoLSat_("rhoLSat", dimDensity, dict),     
    saturationModel_
    (
        saturationModel::New
        (
            dict.subDict("saturationPressure"),
            pair
        )
    )
    
{
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulkNucleationModel::~bulkNucleationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// This has many errors
Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::dmdts2to1(volScalarField& dmdtf ) const
{

     const phaseModel& phase1 = pair_.phase1();
     const phaseModel& phase2 = pair_.phase2();
     const rhoThermo& thermo1 = phase1.thermo();
     const rhoThermo& thermo2 = phase2.thermo();
     const volScalarField& T1(thermo1.T());
     const volScalarField& rho1(thermo1.rho());
     const volScalarField& rho2(thermo2.rho());
     const volScalarField& p(thermo1.p());
     const volScalarField sigma
                     (
                       pair_.phase1().fluid().lookupSubModel<surfaceTensionModel>
                          (
                             phasePair(pair_.phase1(), pair_.phase2())
                          ).sigma()
                     );
            
     volScalarField meshVol  
           (
             IOobject
              (
                 "meshVol",
                  pair_.phase1().mesh().time().timeName(),
                  pair_.phase1().mesh(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
               ),
              pair_.phase1().mesh(),
              dimensionedScalar( dimVolume, small)  
    
           );
    meshVol.ref() =   pair_.phase1().mesh().V(); 
    
    const volScalarField pSat(saturationModel_->pSat(T1));
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number 


    
    volScalarField L(thermo1.ha()-thermo2.ha());
    
      Info<< " Latent heat mean  = " << average(L.primitiveField()) <<endl;
    
    volScalarField rc ((2*sigma)/ ( (1-(rho1/rho2))* (pSat-p))); 
 
    Info<< "rc   min = " << min(rc.primitiveField()) << "  rc   max = " << max(rc.primitiveField()) <<  "  rc  dimensions = " << rc.dimensions() <<endl;

    Info<< "p   min = " << min(p.primitiveField()) << "  p   max = " << max(p.primitiveField()) <<  "  p  dimensions = " << p.dimensions() <<endl;   
     const dimensionedScalar rctemp(dimLength,0.01); // setting all negative values to 0.01
   
     forAll(rc, celli)
    { 
        if(rc[celli] < 0.0) 
        
        {
           rc[celli] =  0.01;      
        }    
    }
      
    Info<< "rc   min = " << min(rc.primitiveField()) << "  rc   max = " << max(rc.primitiveField()) <<  "  rc  dimensions = " << rc.dimensions() <<endl;
 
    volScalarField Wcr ( 4*constant::mathematical::pi*sqr(rc)*sigma/3);    
         
    volScalarField kT (k*T1);
                                                              
    volScalarField WcrkTN ( 
                              (4*constant::mathematical::pi*sqr(rc)*sigma/3)/   
                              (k*T1*(4/3)*constant::mathematical::pi*pow(rc,3)*rho1*(1000/thermo1.W())*Av )                              
                            );   
     
    volScalarField Ja ( 
                           (1000/thermo1.W())*Av* rho2 *
                           B(phase1, phase2,pSat, "dmdts2to1") * //B function
                           exp (-WcrkTN) * 
                           (1-exp(-( sigma*thermo1.W()/ (kT *rho1*1000*Av))/rc))                    // distribution percentage                            
                         );     
 //     Info<< "Ja   min = " << min(Ja.primitiveField()) << "  Ja   max = " << max(Ja.primitiveField()) <<  "  Ja  dimensions = " << Ja.dimensions() <<endl;
 
  
    volScalarField D32vol   
           (
             IOobject
              (
                 "meshVol",
                  pair_.phase1().mesh().time().timeName(),
                  pair_.phase1().mesh(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
               ),
              pair_.phase1().mesh(),
              dimensionedScalar( dimVolume, D32_)  
    
           );
           
           
 
     label  ncvs  = 0;
      scalar nucvol = 0.0;
     forAll(dmdtf, celli)
    {   
    
       if(WcrkTN[celli] > WcrkTNmin_  ) 
         {
           dmdtf[celli] = Ja[celli]*meshVol[celli]*rho1[celli]*D32vol[celli];   
            ++ ncvs ;
           nucvol += meshVol[celli];     
         }
    }
    
    Info << "Total Nucleating CV's = " << ncvs 
          << ",  Total nucleating volume = " << nucvol 
          << endl;
 
    Info  << "WcrkTN values"
           << ": min = " << min(WcrkTN.primitiveField())
           << ", mean = " << average(WcrkTN.primitiveField())
           << ", max = " << max(WcrkTN.primitiveField())
           << endl;
    
     return dmdtf;
}



Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::pSat
(
    const volScalarField& T 
) const
{
    return saturationModel_->pSat(T);
}

Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::Tsat
(
    const volScalarField& p
) const
{
    return saturationModel_->Tsat(p);
}


// This is fixed
 Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::calcWcrkTN( ) const
{
     const phaseModel& phase1 = pair_.phase1();
 //    const phaseModel& phase2 = pair_.phase2();
     const rhoThermo& thermo1 = phase1.thermo();
 //    const rhoThermo& thermo2 = phase2.thermo();
     const volScalarField& T1(thermo1.T());
 //    const volScalarField& rho1(thermo1.rho());
 //    const volScalarField& rho2(thermo2.rho());
     const volScalarField& p(thermo1.p());
     const volScalarField sigma
                     (
                       pair_.phase1().fluid().lookupSubModel<surfaceTensionModel>
                          (
                             phasePair(pair_.phase1(), pair_.phase2())
                          ).sigma()
                     );
 
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number 
      
           
    volScalarField rc((2*sigma)/ ( (1-(rhoVSat_/rhoLSat_))* (pSat_-p))); 
 
  
//    Info<< "gas volume min = " << min(rho1.primitiveField()) << ",  gas volume max = " << max(rho1.primitiveField())  <<endl;     
     forAll(rc, celli)
    { 
        if(rc[celli] < 0) 
        
        {
           rc[celli] =  0.01;      
        }    
    }
    
//   Info<< "rc   min = " << min(rc.primitiveField()) << ",  rc   max = " << max(rc.primitiveField()) <<endl;
   volScalarField Wcr (4*constant::mathematical::pi*sqr(rc)*sigma/3);      
   
     volScalarField Nc(
                        4
                        *constant::mathematical::pi
                        *pow(rc,3.0)
                        *rhoVSat_
                        *1000
                        *Av
                        /(3*thermo1.W()) 
                      );
                                 
    return (Wcr/(Nc*k*T1));
}


Foam::scalar Foam::bulkNucleationModel::WcrkTNmin( ) const
{
    return WcrkTNmin_;
}


Foam::scalar Foam::bulkNucleationModel::WcrFrac( ) const
{
    return WcrFrac_;
}

Foam::scalar Foam::bulkNucleationModel::WcrkTDelta( ) const
{
    return WcrkTDelta_;
}
// ************************************************************************* //

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
    WcrkTNmin2_(dict.lookup<scalar>("WcrkTNmin2",0.03)),
    WcrkTDelta_(dict.lookup<scalar>("WcrkTDelta",0.001)),    
    WcrFrac_(dict.lookupOrDefault<scalar>("WcrFrac",1)),
    dNuc_(dict.lookup<scalar>("dNuc",1e-4)),  
    n_(dict.lookup<scalar>("n",1)), 
    saturationTmodel_
    (
        saturationModel::New
        (
            dict.subDict("saturationTemperature"),
            pair
        )
    ),
    saturationPmodel_
    (
        saturationModel::New
        (
            dict.subDict("saturationPressure"),
            pair
        )
    ),
    saturationRho1model_
    (
        saturationDensityModel::New
        (
            dict.subDict("saturationRho1Density"),
            pair
        )
    ), 
    saturationRho2model_
    (
        saturationDensityModel::New
        (
            dict.subDict("saturationRho2Density"),
            pair
        )
    )             
{
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bulkNucleationModel::~bulkNucleationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 
Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::dmdts2to1(
 const volScalarField& phi
 ) const
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
    
    const volScalarField pSat(saturationPmodel_->pSat(T1));
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number 


    
    volScalarField L(thermo1.ha()-thermo2.ha());
    
//      Info<< " Latent heat mean  = " << average(L.primitiveField()) <<endl;
    
    volScalarField rc ((2*sigma)/ ( (1-(rho1/rho2))* (pSat-p))); 
 
//    Info<< "rc   min = " << min(rc.primitiveField()) << "  rc   max = " << max(rc.primitiveField()) <<  "  rc  dimensions = " << rc.dimensions() <<endl;


   
   
     forAll(rc, celli)
    { 
        if(rc[celli] < 0.0) 
        
        {
           rc[celli] =  0.01;      
        }    
    }
      
    Info<< "rc   min = " << min(rc).value()<< "  rc   max = " << max(rc).value() <<  "  rc  dimensions = " << rc.dimensions() <<endl;
     const volScalarField WcrkTN =     4*constant::mathematical::pi*sqr(rc)*sigma/(3*k*T1*n_);  
     scalar   WcrLimit = WcrkTNmin2_ + WcrkTDelta_;
     

     Info << "WcrkTN min " << min(WcrkTN).value() 
          << ", WcrkTN max " << max(WcrkTN).value()
          <<endl;  
                     
     // blending factor 
     volScalarField factor =  ( WcrLimit - WcrkTN)/ (2.0* WcrkTDelta_) ;  
     factor = min( max (factor,0.0), 1.0);     
     Info<< "factor min = " << min(factor .primitiveField()) << ", factor max = " << max(factor .primitiveField()) <<endl;  

     // calculate the p corressponding to the limit
     
     scalar pNucLimit = -1.0* min( -1.0 *pos0(WcrLimit - WcrkTN ) *p).value() ;     
     Info << "p min " << min(p).value() 
          << ", p max " << max(p).value()
          << ", p NucStart " << pNucLimit 
          <<endl;            
         
         
     const volScalarField WcrkTNPhi =    phi*WcrkTN ;   
       
  // Info<< "WcrkTN    min = " << min(WcrkTN .primitiveField()) << "  WcrkTN    max = " << max(WcrkTN .primitiveField()) <<  "  WcrkTN   dimensions = " << WcrkTN .dimensions() <<endl;         
 //   const volScalarField Ba (B(phase1, phase2,pSat, "dmdts2to1") ); 
    
 //    const volScalarField expWktn= exp(-WcrkTN);
//     Info<< "Ba   min = " << min(Ba.primitiveField()) << "  Ba   max = " << max(Ba.primitiveField()) <<  "  Ba  dimensions = " << Ba.dimensions() <<endl;         
 //      Info<< "expWktn   min = " <<endl;     
       
               
    const volScalarField Ja =    (1000/thermo1.W())*Av* rho2 * B(phase1, phase2,pSat, "dmdts2to1") * exp(-WcrkTNPhi);  
                         
                     
 //    Info<< "Ja   min = " << min(Ja.primitiveField()) << "  Ja   max = " << max(Ja.primitiveField()) <<  "  Ja  dimensions = " << Ja.dimensions() <<endl;
 
    dimensionedScalar dNucVol = (4/3) *constant::mathematical::pi * pow((dNuc_/2),3.0) ;
          
           
//     Info << "line 201" << endl;

     return factor*phase2*Ja*meshVol*rho1*dNucVol ;
}



Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::pSat
(
    const volScalarField& T 
) const
{
    return saturationPmodel_->pSat(T);
}

Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::Tsat
(
    const volScalarField& p
) const
{
    return saturationTmodel_->Tsat(p);
}


// This is fixed
 Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::calcWcrkTN2( ) const
{
     const phaseModel& phase1 = pair_.phase1();
     const phaseModel& phase2 = pair_.phase2();
     const rhoThermo& thermo2 = phase2.thermo();
     const rhoThermo& thermo1 = phase1.thermo();
 //    const rhoThermo& thermo2 = phase2.thermo();
     const volScalarField& T2(thermo2.T());
 //    const volScalarField& rho1(thermo1.rho());
 //    const volScalarField& rho2(thermo2.rho());
     const volScalarField& p(thermo1.p());
     const volScalarField pSat(saturationPmodel_->pSat(T2));
     
     const volScalarField rho1Sat(saturationRho1model_->rhoSat(T2)); 
            
     const volScalarField rho2Sat(saturationRho2model_->rhoSat(T2));  
 
     const volScalarField sigma
                     (
                       pair_.phase1().fluid().lookupSubModel<surfaceTensionModel>
                          (
                             phasePair(pair_.phase1(), pair_.phase2())
                          ).sigma()
                     );
 
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number 
      
           
    volScalarField rc((2*sigma)/ ( (1-(rho1Sat/rho2Sat))*(pSat-p)));  //pSat was pSat_, rho1Sat  was rhoVSat, rho2Sat was rhoLSat_, 
 
  
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
                        *rho1Sat
                        *1000
                        *Av
                        /(3*thermo1.W()) 
                      );
                                 
    return (Wcr/(Nc*k*T2));
}


Foam::scalar Foam::bulkNucleationModel::WcrkTNmin2( ) const
{
    return WcrkTNmin2_;
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

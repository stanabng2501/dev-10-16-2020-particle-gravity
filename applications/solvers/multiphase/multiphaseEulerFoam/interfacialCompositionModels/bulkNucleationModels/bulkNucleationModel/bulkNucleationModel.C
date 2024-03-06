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
    WcrkTNmin2_(dict.lookupOrDefault<scalar>("WcrkTNmin2",0.03)),
    WcrkTDelta_(dict.lookupOrDefault<scalar>("WcrkTDelta",0.001)),    
    WcrFrac_(dict.lookupOrDefault<scalar>("WcrFrac",1)),
    dNuc_(dict.lookupOrDefault<scalar>("dNuc",1e-4)),  
    n_(dict.lookupOrDefault<scalar>("n",1)),
    bFactor_(dict.lookupOrDefault<bool>("bubbleFactor",false)),  
    alpha1min_(dict.lookupOrDefault<scalar>("minBlendedAlpha1",0.0)),
    alpha1max_(dict.lookupOrDefault<scalar>("maxBlendedAlpha1",1.0)),   
    alpha2min_(dict.lookupOrDefault<scalar>("minBlendedAlpha2",0.0)),
    alpha2max_(dict.lookupOrDefault<scalar>("maxBlendedAlpha2",1.0)), 
    residualAlpha_(dict.lookupOrDefault<scalar>("residualAlpha",1e-6)), 
    dAlpha_(dict.lookupOrDefault<scalar>("dAlpha",1e-5)),    
    dWcr_(dict.lookupOrDefault<scalar>("dWcr",1e-5)),    
    blend_(dict.lookupOrDefault<bool>("blending",false)),
    axi_(dict.lookupOrDefault<bool>("axisymmetric",false)), 
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



Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::CalcWcrKTN2to1() const
 
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
              dimensionedScalar( dimVolume, 0.0)  
    
           );
    meshVol.ref() =   pair_.phase1().mesh().V(); 
/*    
    Info<< "meshVol   min = " << min(meshVol).value() 
        << "  meshVol  max = " << max(meshVol).value() 
        << endl;

    Info<< "meshVol   min = " << min(meshVol.ref() ) 
        << "  meshVol  max = " << max(meshVol.ref() )  
        << endl; 
*/           
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
      
 //   Info<< "rc   min = " << min(rc).value()<< "  rc   max = " << max(rc).value() <<  "  rc  dimensions = " << rc.dimensions() <<endl;
     const volScalarField WcrkTN =     4*constant::mathematical::pi*sqr(rc)*sigma/(3*k*T1*n_);   
      scalar   WcrLimit = WcrkTNmin2_ + WcrkTDelta_;
      
      return min(WcrkTN,WcrLimit);
 
 }


// this is for gas to liquid phase change
Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::dmdts1to2(
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
              dimensionedScalar( dimVolume, 0.0)  
    
           );
    meshVol.ref() =   pair_.phase1().mesh().V(); 
    
    const volScalarField pSat(saturationPmodel_->pSat(T1));
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number 
    
    volScalarField rc ((2*sigma)/ ( ((rho2/rho1) -1)* (p - pSat)));                      
              
    const dimensionedScalar rcMax( dimLength,0.01);
    rc =neg0(rc)*rcMax + pos(rc)*rc;    

    const volScalarField WcrkTN = 4*constant::mathematical::pi*sqr(rc)*sigma/(3*k*T1*n_);  


    volScalarField dropletFactor = phase2;  
    dropletFactor.max(residualAlpha_); 
    
   dropletFactor = pos0(dropletFactor - alpha1min_)* neg0(dropletFactor - alpha1max_)*dropletFactor;     
   
   const volScalarField WcrkTNPhi =    phi*WcrkTN ;              
   const volScalarField Ja =    (1000/thermo1.W())*Av* rho1 * B(phase1, phase2,pSat, "dmdts2to1") * exp(-WcrkTNPhi);  
   dimensionedScalar dNucVol = (4/3) *constant::mathematical::pi * pow((dNuc_/2),3.0) ;      
   Info << "dropletFactor 1 to 2 : min = " << min(dropletFactor).value() 
        <<  "     max = " << max(dropletFactor).value()  << endl;
        
        
   return dropletFactor * phase1*Ja*meshVol*rho2*dNucVol ; // last three terms number of liquid droplets X mass of each droplets

}

// this is for liquid to gas phase change 
Foam::tmp<Foam::volScalarField> Foam::bulkNucleationModel::dmdts2to1(
 const volScalarField& phi
 ) const
{

     const phaseModel& phase1 = pair_.phase1();
     const phaseModel& phase2 = pair_.phase2();
     const rhoThermo& thermo1 = phase1.thermo();
     const rhoThermo& thermo2 = phase2.thermo();
     const volScalarField& T2(thermo2.T());
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
              dimensionedScalar( dimVolume, 0.0)  
    
           );
    meshVol.ref() =   pair_.phase1().mesh().V(); 


// Linear blending
   
       
    const volScalarField pSat(saturationPmodel_->pSat(T2));
    const dimensionedScalar Av(Foam::constant::physicoChemical::NA) ; // avagadro number
    const dimensionedScalar k(Foam::constant::physicoChemical::k) ; // Boltzmann constant number 


    
    volScalarField L(thermo1.ha()-thermo2.ha());
    
//      Info<< " Latent heat mean  = " << average(L.primitiveField()) <<endl;
    
    volScalarField rc ((2*sigma)/ ( (1-(rho1/rho2))* (pSat-p))); 
 
//    Info<< "rc   min = " << min(rc.primitiveField()) << "  rc   max = " << max(rc.primitiveField()) <<  "  rc  dimensions = " << rc.dimensions() <<endl;

     const dimensionedScalar rcMax( dimLength,0.01);
    rc =neg0(rc)*rcMax + pos(rc)*rc;
   
/*   
     forAll(rc, celli)
    { 
        if(rc[celli] < 0.0) 
        
        {
           rc[celli] =  0.01;      
        }    
    }
 */ 
      
 //   Info<< "rc   min = " << min(rc).value()<< "  rc   max = " << max(rc).value() <<  "  rc  dimensions = " << rc.dimensions() <<endl;
     const volScalarField WcrkTN =     4*constant::mathematical::pi*sqr(rc)*sigma/(3*k*T2*n_);  


    volScalarField bubbleFactor = phase1;


    if( max(bubbleFactor).value() > alpha1min_  && axi_  ) 
    
    {
        const scalar alpha1dmax = alpha1min_ +dAlpha_;
        
        volScalarField WcrkTNTemp = neg(bubbleFactor - alpha1min_)*1e6 + pos0(bubbleFactor - alpha1min_)*WcrkTN;
 
            
   
        const scalar WcrkTNTempMin =  min(WcrkTNTemp).value();
        
//        WcrkTNTemp = neg(bubbleFactor - alpha1dmax)*1e6 + pos0(bubbleFactor - alpha1dmax)*WcrkTN;
        const scalar WcrkTNTempMax =  WcrkTNTempMin +dWcr_; 
        
        Info << "WcrkTNTempMin, WcrkTNTempMax  = "<< WcrkTNTempMin << " , " <<  WcrkTNTempMax <<endl; 


 //        Info << "bubbleFactor min, max = "<<  min(bubbleFactor .primitiveField()) << " , " <<         max(bubbleFactor .primitiveField())<<endl;        
 
        bubbleFactor = neg(WcrkTN - WcrkTNTempMin)*alpha1min_ +
                       pos0(WcrkTN - WcrkTNTempMin) * neg(WcrkTN - WcrkTNTempMax)* alpha1dmax + 
                       pos0(WcrkTN - WcrkTNTempMax)*bubbleFactor; 

//         Info << "bubbleFactor min, max = "<<  min(bubbleFactor .primitiveField()) << " , " <<         max(bubbleFactor .primitiveField())<<endl; 
        
                
/*        bubbleFactor =   neg0(WcrkTN - WcrkTNTempMax) * alpha1dmax;
        Info << "bubbleFactor min, max = "<< min(bubbleFactor).value() << " , " <<  max(bubbleFactor).value() <<endl;  
        bubbleFactor =  pos(WcrkTN - WcrkTNTempMin)*bubbleFactor;
*/
        
            
    
    }
    
    if(blend_)
    
    {
   
      bubbleFactor = neg0(bubbleFactor - alpha1min_)* alpha1min_ +  pos(bubbleFactor - alpha1min_)* bubbleFactor; 
      bubbleFactor = pos0(bubbleFactor - alpha1max_)* alpha1max_ + neg(bubbleFactor - alpha1max_)* bubbleFactor ; 
      bubbleFactor =  (bubbleFactor -  alpha1min_)/  (alpha1max_ -  alpha1min_);

        Info << "blended bubbleFactor min, max = "<< min(bubbleFactor).value() << " , " <<  max(bubbleFactor).value() <<endl;  
    }

    else

      {
	bubbleFactor = pos0(bubbleFactor - alpha1min_)* neg0(bubbleFactor - alpha1max_)*bubbleFactor;
       
       Info << "Un-blended bubbleFactor min, max = "<< min(bubbleFactor).value() << " , " <<  max(bubbleFactor).value() <<endl;
      }	
    
    bubbleFactor.max(residualAlpha_); 
     scalar   WcrLimit = WcrkTNmin2_ + WcrkTDelta_;
     

     Info << "WcrkTN min " << min(WcrkTN).value() 
          << ", WcrkTN max " << max(WcrkTN).value()
          <<endl;  
                     
     // blending factor 
     volScalarField factor =  ( WcrLimit - WcrkTN)/ (2.0* WcrkTDelta_) ;  
     factor = min( max (factor,0.0), 1.0);     
//     Info<< "factor min = " << min(factor).value() << ", factor max = " <<  max(factor).value() <<endl;  

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
                         
                     
//   Info<< "Ja   min = " << min(Ja.primitiveField()) << "  Ja   max = " << max(Ja.primitiveField()) <<  "  Ja  dimensions = " << Ja.dimensions() <<endl;
 
    dimensionedScalar dNucVol = (4/3) *constant::mathematical::pi * pow((dNuc_/2),3.0) ;
          
           
//     Info << "line 201" << endl;
      if (bFactor_) 
      {
         Info << "bubbleFactor 2 to 1 : min = " << min(bubbleFactor).value() 
              <<  "     max = " << max(bubbleFactor).value()  << endl;
        return bubbleFactor * factor*phase2*meshVol*Ja*rho1*dNucVol ; // last three terms number of gas bubbles X mass of each bubble
      }
      else
      {     
        return factor*phase2*Ja*meshVol*rho1*dNucVol ;
      }
     
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
     const dimensionedScalar rcMax( dimLength,0.01);
    rc =neg0(rc)*rcMax + pos(rc)*rc;
   
/*   
     forAll(rc, celli)
    { 
        if(rc[celli] < 0.0) 
        
        {
           rc[celli] =  0.01;      
        }    
    }
 */ 
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

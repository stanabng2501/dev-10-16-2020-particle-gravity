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

#include "LimitedEvaporationPhaseChangePhaseSystem.H"
#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fvcVolumeIntegrate.H"
#include "fvmSup.H"
#include "rhoReactionThermo.H"
#include "heThermo.H"
// #include "fvcDdt.H"
// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::addDmdts
(
    PtrList<volScalarField>& dmdts
) const
{
 
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {

        const phasePair& pair = this->phasePairs_[dmdtfIter.key()];
        const volScalarField& dmdtf = *dmdtfIter();

        addField(pair.phase1(), "dmdt", dmdtf, dmdts);
        addField(pair.phase2(), "dmdt", - dmdtf, dmdts);
    }


}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::
LimitedEvaporationPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "saturation",
        saturationModels_
    );

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[saturationModelIter.key()];

        if
        (
            !this->heatTransferModels_.found(pair)
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for the " << pair
                 << "pair is not specified. This is required by the "
                 << "corresponding saturation model"
                 << exit(FatalError);
        }
    }

    // Generate interfacial mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[saturationModelIter.key()];

        dmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "limitedEvaporationPhaseChange:dmdtf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );
        
        H1_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "limitedEvaporationPhaseChange:H1",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimensionSet(1, -1, -3, -1, 0), 0)
            )
        );
        
        H2_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "limitedEvaporationPhaseChange:H2",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimensionSet(1, -1, -3, -1, 0), 0)
            )
        );  

    }
    
    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::
~LimitedEvaporationPhaseChangePhaseSystem()
{}


// * * * * * * *Energy/dimMass * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
const Foam::bulkNucleationModel&
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::saturation
(
    const phasePairKey& key
) const
{
    return saturationModels_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    const scalar dmdtfSign =
        Pair<word>::compare(this->phasePairs_[key], key);

    if (dmdtfs_.found(key))
    {
        tDmdtf.ref() += dmdtfSign**dmdtfs_[key];
    }


    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{

   /// calculating dmdtfs_
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[saturationModelIter.key()];
        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();
        const volScalarField& T1(thermo1.T());
        const volScalarField& T2(thermo2.T());
        const volScalarField& p(thermo1.p());
  //      const volScalarField  dpdt =  fvc::ddt(phase1.thermo().p()); 
  //     const volScalarField&  dpdt(phase1.fluid().dpdt());

        // Interfacial mass transfer update
        {
            Info << "---------------------------------------------------------------------------------------------" <<endl;            
            volScalarField& dmdtf(*this->dmdtfs_[pair]);
            volScalarField L(this->L(pair, dmdtf, T2, latentHeatScheme::symmetric));

            const volScalarField Tsat(saturationModelIter()->Tsat(p));                 
            const volScalarField WcrkTN2(saturationModelIter()->calcWcrkTN2( ));   
            volScalarField H1(this->heatTransferModels_[pair].first()->K());
            volScalarField H2(this->heatTransferModels_[pair].second()->K());

            *H1_[pair] = H1;
            *H2_[pair] = H2;
            // Limit the H[12] to avoid /0
            H1.max(small);
            H2.max(small);
/*
            volScalarField dmdtfNew 
                             (
                                IOobject
                                (
                                   IOobject::groupName("dmdtfNew", pair.name()),
                                    this->mesh().time().timeName(),
                                    this->mesh(),
                                    IOobject::NO_READ,
                                    IOobject::NO_WRITE
                                 ),
                                 this->mesh(),
                                 dimensionedScalar( dmdtf.dimensions(), 0.0)     
                              );
 */                             
            volScalarField meshVol  
             (
               IOobject
                (
                 "meshVol",
                 this->mesh().time().timeName(),
                 this->mesh(),
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
                ),
              this->mesh(),
              dimensionedScalar( dimVolume, small)     
             );
            meshVol.ref() =    phase1().mesh().V(); 
 
            scalar   WcrLimit;         
            const scalar dmdtfRelaxRem =
                this->mesh().fieldRelaxationFactor(dmdtf.member()+ "Removal"); 
            const scalar dmdtfRelaxAdd =
                this->mesh().fieldRelaxationFactor(dmdtf.member()+ "Addition");   
            const  scalar& WcrkTDelta(saturationModelIter()->WcrkTDelta());            
            const  scalar& WcrkTNmin2(saturationModelIter()->WcrkTNmin2());
            
            WcrLimit = WcrkTNmin2 - WcrkTDelta;
           
            // blending factor 
            volScalarField factor =  (WcrkTN2 - WcrLimit)/ (2.0* WcrkTDelta) ;  
            factor = max( min (factor,1.0), 0.0);
 
            Info << "Factor   min = " << min(factor.primitiveField())
                << ",     max = " << max(factor.primitiveField())
                 << ",    mean = " << average(factor.primitiveField()) 
                << endl;           
 
/*            
            const  scalar& WcrFrac(saturationModelIter()->WcrFrac());                        
 //         const scalar dpdtvalue =  (1-dpdtLim ) * min(dpdt.primitiveField()) ;
                
         // setting Wcr Limits
           if(WcrFrac <1)
             {
                WcrLimit  =max( WcrkTNmin2 , (1-WcrFrac)*max(WcrkTN.primitiveField()) ) ;
             }
             
           else
             {
               WcrLimit = WcrkTNmin2;
             }                     
 */         
 
   
            label  ncvs  = 0;
            scalar nucvol = 0.0; 
            scalar dmdtfNewVal = 0.0;                            
            forAll(dmdtf, celli)
             {     
               if(WcrkTN2[celli] > WcrLimit ) 
               {
               
                 dmdtfNewVal =  factor[celli]*(H1[celli]*(Tsat[celli] - T1[celli]) + H2[celli]*(Tsat[celli] - T2[celli]))/L[celli]; 
                 dmdtf[celli] = (1 - dmdtfRelaxAdd)*dmdtf[celli] + dmdtfRelaxAdd*dmdtfNewVal;  
                 ++ ncvs ;
                 nucvol += meshVol[celli];     
                }
                
                else
                {
                  dmdtf[celli]  = dmdtfRelaxRem*dmdtf[celli];
                }
              }
                                                             
            
   //          dmdtf = (1 - dmdtfRelax)*dmdtf + dmdtfRelax*dmdtfNew;
  
           Info << "Total Nucleating CV's = " << ncvs 
                << ",  Total nucleating volume = " << nucvol 
                << endl;

           Info << "p min = " << min(p.primitiveField()) 
                << ",  p max = " << max(p.primitiveField())
                <<endl;
                
           Info << "WcrkTN 2   min = " << min(WcrkTN2.primitiveField())
                << ", WcrkTN 2  max = " << max(WcrkTN2.primitiveField())
                << ", Wcr Limit = " << WcrLimit
                << endl;
                
              Info<< "Latent Heat  min = " << min(L.primitiveField())
                 << ", mean = " << average(L.primitiveField())
                  << ", max = " << max(L.primitiveField())
                 << endl;
                
                 
            Info<< dmdtf.name()
                << ": min = " << min(dmdtf.primitiveField())
                << ", mean = " << average(dmdtf.primitiveField())
                << ", max = " << max(dmdtf.primitiveField())
                << ", sum = " << sum (dmdtf.primitiveField()) //fvc::domainIntegrate(dmdtf).value()
                << endl;
                
            Info << "---------------------------------------------------------------------------------------------" <<endl;    
        }
     }

    // adding dmdtfs_ to dmdts()
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());
    addDmdts(dmdts);
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}




template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::specieTransfer() const
{

    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();
 
        this->addDmdtYfs(dmdtfs_, eqns);
   

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::
correctContinuityError()
{
   // no continuity error yet bcoz no heat transfer
}

template<class BasePhaseSystem>
void Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
   // no interface thermo to correct
}

 



template<class BasePhaseSystem>
bool Foam::LimitedEvaporationPhaseChangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

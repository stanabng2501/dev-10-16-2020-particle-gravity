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

#include "NucleationPhaseChangePhaseSystem.H"
#include "fvcVolumeIntegrate.H"
#include "fvmSup.H"
#include "rhoReactionThermo.H"
#include "heThermo.H"
// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::addDmdts
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
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
NucleationPhaseChangePhaseSystem
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
                        "evaporationPhaseChange:dmdtf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

 
       
 
//      Info << "max pressure value is " <<max(pair.phase2().thermo().p().primitiveField()) << endl; 
 //     Info << "max temperature  value is " <<max(pair.phase2().thermo().T().primitiveField()) << endl;
    }
    
    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
~NucleationPhaseChangePhaseSystem()
{}


// * * * * * * *Energy/dimMass * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
const Foam::saturationModel&
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::saturation
(
    const phasePairKey& key
) const
{
    return saturationModels_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
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
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
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

        // Interfacial mass transfer update
        {
        
            volScalarField& dmdtf(*this->dmdtfs_[pair]);
            volScalarField L(this->L(pair, dmdtf, T2, latentHeatScheme::symmetric));

            const volScalarField Tsat(saturationModelIter()->Tsat(thermo1.p()));  
            volScalarField H1(this->heatTransferModels_[pair]->K());
            volScalarField H2(this->heatTransferModels_[pair]->K());
            volScalarField dmdtfNew( (H1*(Tsat-T1)+H2*(Tsat-sT2))/L);
   
            // Limit the H[12] to avoid /0
            H1.max(small);
            H2.max(small);


            const scalar dmdtfRelax =
                this->mesh().fieldRelaxationFactor(dmdtf.member());

            dmdtf = (1 - dmdtfRelax)*dmdtf + dmdtfRelax*dmdtfNew;


            Info<< "Latent Heat  min = " << min(L.primitiveField())
                << ", mean = " << average(L.primitiveField())
                << ", max = " << max(L.primitiveField())
                << endl;
                
                
            Info<< dmdtf.name()
                << ": min = " << min(dmdtf.primitiveField())
                << ", mean = " << average(dmdtf.primitiveField())
                << ", max = " << max(dmdtf.primitiveField())
                << ", integral = " << fvc::domainIntegrate(dmdtf).value()
                << endl;
        }
     }

    // adding dmdtfs_ to dmdts()
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());
    addDmdts(dmdts);
    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
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
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
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
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::specieTransfer() const
{

    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();
 
        this->addDmdtYfs(dmdtfs_, eqns);
   

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
correctContinuityError()
{
   // no continuity error yet bcoz no heat transfer
}

template<class BasePhaseSystem>
void Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
   // no interface thermo to correct
}

 



template<class BasePhaseSystem>
bool Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::read()
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

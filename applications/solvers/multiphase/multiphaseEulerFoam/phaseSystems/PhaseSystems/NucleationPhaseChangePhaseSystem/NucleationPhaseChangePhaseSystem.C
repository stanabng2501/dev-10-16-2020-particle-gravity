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
        "nucleation",
        nucleationModels_
    );

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        nucleationModelTable,
        nucleationModels_,
        nucleationModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[nucleationModelIter.key()];

        if
        (
            !this->heatTransferModels_.found(pair)
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for the " << pair
                 << "pair is not specified. This is required by the "
                 << "corresponding nucleation model"
                 << exit(FatalError);
        }
    }

    // Generate interfacial mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        nucleationModelTable,
        nucleationModels_,
        nucleationModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[nucleationModelIter.key()];

        dmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "nucleationPhaseChange:dmdtf",
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

 
    }
    
    
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::
~NucleationPhaseChangePhaseSystem()
{}


// * * * * * * *Energy/dimMass * * * * * * * Member Functions  * * * * * * * * * * * * * * //

/*
template<class BasePhaseSystem>
const Foam::bulkNucleationModel&
Foam::NucleationPhaseChangePhaseSystem<BasePhaseSystem>::nucleation
(
    const phasePairKey& key
) const
{
    return nucleationModels_[key];
}
*/

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
        nucleationModelTable,
        nucleationModels_,
        nucleationModelIter
    )
    {
          const phasePair& pair = this->phasePairs_[nucleationModelIter.key()];

          volScalarField& dmdtf(*this->dmdtfs_[pair]); 
          volScalarField dmdtfNew 
            (
             IOobject
              (
                 "dmdtfNew ",
                  pair.phase1().mesh().time().timeName(),
                  pair.phase1().mesh(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
               ),
              pair.phase1().mesh(),
              dimensionedScalar( dmdtf.dimensions(), Zero)  
    
           );
        
         const scalar dmdtfRelax =
                this->mesh().fieldRelaxationFactor(dmdtf.member());
         const dimensionedScalar dmdtfMax ( dmdtf.dimensions(),
                                            this->mesh().fieldRelaxationFactor(dmdtf.member()+ "Max") 
                                          );         
        forAllConstIter(phasePair, pair, pairIter)
        {
            const label sign = pairIter.index() == 0 ? 1 : -1;
            const autoPtr<bulkNucleationModel>& bulkNucleationModelPtr =
                nucleationModelIter()[pairIter.index()];
                
           if (!bulkNucleationModelPtr.valid()) continue;
           
           if (sign ==1)
            {
            const volScalarField dmdtfNew1(nucleationModels_[pair][pairIter.index()]->dmdts2to1());
            dmdtfNew +=sign*dmdtfNew1;
            }
            
            else
            {
        //    replace this by  dmdts1to2         
            const volScalarField dmdtfNew2(nucleationModels_[pair][pairIter.index()]->dmdts2to1());
            dmdtfNew +=sign*dmdtfNew2;
            }
                      

        }

            dmdtf = (1 - dmdtfRelax)*dmdtf + dmdtfRelax*dmdtfNew;
            dmdtf = min( dmdtf, dmdtfMax);
                
                
            Info<< dmdtf.name()
                << ": min = " << min(dmdtf.primitiveField())
                << ", mean = " << average(dmdtf.primitiveField())
                << ", max = " << max(dmdtf.primitiveField())
                << ", sum = " << sum(dmdtf).value()
                << endl;
                
            Info << " ------------------------------------------------" << endl;            
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

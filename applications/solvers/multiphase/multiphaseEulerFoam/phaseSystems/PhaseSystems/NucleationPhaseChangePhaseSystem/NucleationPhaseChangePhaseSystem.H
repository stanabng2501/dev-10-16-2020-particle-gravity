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

Class
    Foam::NucleationChangePhaseSystem

Description
    Class to provide  mass transfer between a number of
    phases according the nucleation theory




SourceFiles
    NucleationPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef NucleationPhaseChangePhaseSystem_H
#define NucleationPhaseChangePhaseSystem_H

#include "phaseSystem.H"
#include "bulkNucleationModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ThermalPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class NucleationPhaseChangePhaseSystem
:
    public BasePhaseSystem
{
    // Private typedefs

        typedef HashTable
        <
            Pair<autoPtr<bulkNucleationModel>>,
            phasePairKey,
            phasePairKey::hash
        > nucleationModelTable;

        using latentHeatScheme = typename BasePhaseSystem::latentHeatScheme;



    // Private data

        //- Name of the volatile specie
        word MassMomentum_;

        //- nucleation models used to evaluate Tsat = Tf
        nucleationModelTable nucleationModels_;

        //- Actual Mass transfer rates
        phaseSystem::dmdtfTable dmdtfs_;

        // Calculated mass transfer rates
         phaseSystem::dmdtfTable dmdtfsC_;

        // Calculated mass transfer rates
         phaseSystem::dmdtfTable WcrKTN2to1_;
         
                  
        //- Heteregenity factor
        phaseSystem::dmidtfTable hetFac_; 

      // Private Member Functions

        //- Sum the mass transfer rates for each phase into the given list
        void addDmdts(PtrList<volScalarField>&) const;


public:

    // Constructors

        //- Construct from fvMesh
        NucleationPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~NucleationPhaseChangePhaseSystem();


    // Member Functions

        //- Return the nucleationModel
 //       const bulkNucleationModel& nucleation(const phasePairKey& key) const;

        //- Return the mass transfer rate for an interface
        virtual tmp<volScalarField> dmdtf(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the momentum transfer matrices for the cell-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Return the momentum transfer matrices for the face-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransferf();


        //- Return the specie transfer matrices
        virtual autoPtr<phaseSystem::specieTransferTable>
            specieTransfer() const;

        //- Store phase dmdts at the during the continuity error update
        virtual void correctContinuityError();

        //- Correct the interface thermodynamics
        virtual void correctInterfaceThermo();
        
        
        //- Read base phaseProperties dictionary
        virtual bool read();
        
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NucleationPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

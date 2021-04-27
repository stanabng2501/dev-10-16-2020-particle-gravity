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
    Foam::EvaporationChangePhaseSystem

Description
    Class to provide interfacial heat and mass transfer between a number of
    phases according the interfacial temperature approximated by the saturation
    temperature.




SourceFiles
    EvaporationPhaseChangePhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef EvaporationPhaseChangePhaseSystem_H
#define EvaporationPhaseChangePhaseSystem_H

#include "phaseSystem.H"
#include "saturationModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ThermalPhaseChangePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class EvaporationPhaseChangePhaseSystem
:
    public BasePhaseSystem
{
    // Private typedefs

        typedef HashTable
        <
            autoPtr<saturationModel>,
            phasePairKey,
            phasePairKey::hash
        > saturationModelTable;

        using latentHeatScheme = typename BasePhaseSystem::latentHeatScheme;



    // Private data



        //- Saturation models used to evaluate Tsat = Tf
        saturationModelTable saturationModels_;

        //- Mass transfer rates
        phaseSystem::dmdtfTable dmdtfs_;

        //- Heat transfer rates 
        phaseSystem::dmdtfTable H_;
        
        //- Latent heat
        phaseSystem::dmdtfTable L_;

        //- Nucleate thermal energy transfer rates
        phaseSystem::dmdtfTable nDmdtLfs_;


      // Private Member Functions

        //- Sum the mass transfer rates for each phase into the given list
        void addDmdts(PtrList<volScalarField>&) const;


public:

    // Constructors

        //- Construct from fvMesh
        EvaporationPhaseChangePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~EvaporationPhaseChangePhaseSystem();


    // Member Functions

        //- Return the saturationModel
        const saturationModel& saturation(const phasePairKey& key) const;

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


        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EvaporationPhaseChangePhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

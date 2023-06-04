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
    Foam::bulkNucleationModel

Description
    Generic base class for bulk nucleation kinetic models.  

SourceFiles
    interfaceCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef bulkNucleationModel_H
#define bulkNucleationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "saturationModel.H"
#include "saturationDensityModel.H"

namespace Foam
{

class phaseModel;
class phasePair;

/*---------------------------------------------------------------------------*\
                  Class interfaceCompositionModel Declaration
\*---------------------------------------------------------------------------*/

class bulkNucleationModel
{
    // Private data

        //- Phase pair
        const phasePair& pair_;
        
        
        //- Critical work limit per molecule min max
        const scalar WcrkTNmin2_;

        //-Range for blending
        const scalar WcrkTDelta_;
        //- Depressurization Limit below which nucleation will happen
        const scalar WcrFrac_;
        
        //- Nucleating  Diameter  
         const scalar dNuc_;
         
         //- Cluster size in the bubble
         const scalar n_;
 
         //- bubbleFactor in nucleation, bubbles promoting nucleaiton
         const bool bFactor_; 

         //- bubbleFactor in nucleation, bubbles promoting nucleaiton
         const bool alpha1min_; 
         
         //- bubbleFactor in nucleation, bubbles promoting nucleaiton
         const bool alpha1max_;  
         
         //- Residual Alpha
         const scalar residualAlpha_; 
    
         // Hetereogenity factor
   //      const scalar phi_;
   
        // mesh volume
 //       volScalarField meshVol_;
        
         
                                
protected:

        //- Saturation Temperature model
        autoPtr<saturationModel> saturationTmodel_; 
        
         //- Saturation Pressuree model
        autoPtr<saturationModel> saturationPmodel_;       
               
        //- Saturation gas density as function of T
        autoPtr<saturationDensityModel> saturationRho1model_; 
        
        //- Saturation liquid density as function of T
        autoPtr<saturationDensityModel> saturationRho2model_; 
        
        
   //       label ncvs_;

public:

    //- Runtime type information
    TypeName("bulkNucleationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            bulkNucleationModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Constructors

        //- Construct from a dictionary and a phase pair
        bulkNucleationModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~bulkNucleationModel();


    // Selectors

        static autoPtr<bulkNucleationModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        // Access

 

         //- Saturation pressure
            tmp<volScalarField> pSat
             (
               const volScalarField& T
             ) const;

    

        //- Saturation temperature
             tmp<volScalarField> Tsat
              (
                const volScalarField& p
              ) const;
 
 
         //- Minimum Critical work
             tmp<volScalarField> calcWcrkTN2
              (
              
              ) const;
              
 
          //-Depressurisation Limit
            scalar WcrFrac () const;             
                
         // return mass transfer from liquid too gas
         tmp<Foam::volScalarField> dmdts2to1 ( const volScalarField& phi) const;
         
        tmp<Foam::volScalarField> CalcWcrKTN2to1() const;
               
         //- B factor
            virtual tmp<volScalarField> B
            (
                const phaseModel& phase1 ,
                 const phaseModel& phase2 ,
                 const volScalarField& pSat,
                const word& transfertype 
                
            ) const = 0; 
                 
        //- Return the minimum critical work
        scalar WcrkTNmin2() const;  
        
       //- Return the range for blending
        scalar WcrkTDelta() const; 
               
        // Evaluation

            //- Interface mass fraction
           

            //- The interface mass fraction derivative w.r.t. temperature
             

            //- Mass fraction difference between the interface and the field
            
            
            
            //- Mass fraction difference between the interface and the field
            //  derivative w.r.t. temperature
            

            //- Mass diffusivity
             


        //- Update the composition
         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

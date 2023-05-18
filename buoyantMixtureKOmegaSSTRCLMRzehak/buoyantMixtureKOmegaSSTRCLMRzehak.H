/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
    Foam::RASModels::buoyantMixtureKOmegaSSTRCLMRzehak

Description
    Mixture-based Langtry-Menter transition k-omega SST model with Alahmadi's rotation/curvature correction,
    Devolder's buoyancy production, and Rzehak's bubble induced turbulence.
    

    Langtry-Menter transition SST model described in:

    \verbatim
        Langtry, R. B., & Menter, F. R. (2009).
        Correlation-based transition modeling for unstructured parallelized
        computational fluid dynamics codes.
        AIAA journal, 47(12), 2894-2906.

        Menter, F. R., Langtry, R., & Volker, S. (2006).
        Transition modelling for general purpose CFD codes.
        Flow, turbulence and combustion, 77(1-4), 277-303.

        Langtry, R. B. (2006).
        A correlation-based transition model using local variables for
        unstructured parallelized CFD codes.
        Phd. Thesis, Universität Stuttgart.
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegaSSTCoeffs
        {
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;

            // Langtry-Menter coefficients
            ca1         2;
            ca2         0.06;
            ce1         1;
            ce2         50;
            cThetat     0.03;
            sigmaThetat 2;

            lambdaErr   1e-6;
            maxLambdaIter 10;

            // Rotation/curvature correction
            cr1         1.0;
            cr2         2.0;
            cr3         1.0;

            //Buoyancy production term
            sigmaT      0.85;

            // Bubble generated turbulence
            Cp          1.0; 
        }
    \endverbatim

SourceFiles
    buoyantMixtureKOmegaSSTRCLMRzehak.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantMixtureKOmegaSSTRCLMRzehak_H
#define buoyantMixtureKOmegaSSTRCLMRzehak_H

#include "buoyantMixtureKOmegaSSTRCRzehak.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                Class buoyantMixtureKOmegaSSTRCLMRzehak Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class buoyantMixtureKOmegaSSTRCLMRzehak
:
public buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>
{
   // Private Data
        mutable buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>
            *liquidTurbulencePtr_;
        
    // Private Member Functions
        //- Return the turbulence model for the other phase
        buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>& liquidTurbulence() const;

protected:
    // Protected data

        // Langtry-Menter coefficients
            dimensionedScalar ca1_;
            dimensionedScalar ca2_;

            dimensionedScalar ce1_;
            dimensionedScalar ce2_;

            dimensionedScalar cThetat_;
            dimensionedScalar sigmaThetat_;

            //- Convergence criterion for the lambda/thetat loop
            scalar lambdaErr_;

            //- Maximum number of iterations to converge the lambda/thetat loop
            label maxLambdaIter_;

            //- Stabilization for division by the magnitude of the velocity
            const dimensionedScalar deltaU_;
            

        // Fields

            //- Transition onset momentum-thickness Reynolds number
            volScalarField ReThetat_;

            //- Intermittency
            volScalarField gammaInt_;

            // Mixture fields
            autoPtr<volScalarField> ReThetatm_;
            autoPtr<volScalarField> gammaIntm_;
            autoPtr<volScalarField> gammaIntEffm_;

    // Protected Member Functions

        //- Modified form of the buoyantMixtureKOmegaSSTRCRzehak F1 function
        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega, const volScalarField& kp, const volScalarField& omegap, const volScalarField& nup) const;

        //- Modified form of the buoyantMixtureKOmegaSSTRCRzehak k production rate
        virtual tmp<volScalarField> Pk
        (
            const volScalarField& G
        ) const;

        //- Modified form of the buoyantMixtureKOmegaSSTRCRzehak epsilon/k
        virtual tmp<volScalarField> epsilonByk
        (
            const volScalarField& F1,
            const volScalarField& F2
        ) const;

        virtual void initMixtureFields();

        //- Freestream blending-function
        tmp<volScalarField> Fthetat
        (
            const volScalarField& Usp,
            const volScalarField& Omegap,
            const volScalarField& omegap,
            const volScalarField& nu,
            const volScalarField& gammaIntp,
            const volScalarField& ReThetatp
        ) const;

        //- Empirical correlation for critical Reynolds number where the
        //  intermittency first starts to increase in the boundary layer
        tmp<volScalarField> ReThetac() const;

        //- Empirical correlation that controls the length of the
        //  transition region
        tmp<volScalarField> Flength
        (
            const volScalarField& nup
        ) const;
        
        //- Transition onset location control function
        tmp<volScalarField> Fonset
        (
            const volScalarField& Revp,
            const volScalarField& RTp
        ) const;

        //- Return the transition onset momentum-thickness Reynolds number
        // (based on freestream conditions)
        tmp<volScalarField> ReThetat0
        (
            const volScalarField& kp,
            const volScalarField& Usp,
            const volScalarField& dUsdsp,
            const volScalarField& nup
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correctReThetatGammaInt();

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("buoyantMixtureKOmegaSSTRCLMRzehak");


    // Constructors

        //- Construct from components
        buoyantMixtureKOmegaSSTRCLMRzehak
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        buoyantMixtureKOmegaSSTRCLMRzehak(const buoyantMixtureKOmegaSSTRCLMRzehak&) = delete;


    //- Destructor
    virtual ~buoyantMixtureKOmegaSSTRCLMRzehak()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();
    
    // Langtry-Menter Member Functions

        //- Return the effective diffusivity for transition onset
        //  momentum-thickness Reynolds number
        tmp<volScalarField> DReThetatEff(const volScalarField& nut, const volScalarField& nu) const
        {
            return volScalarField::New
            (
                "DReThetatEff",
                sigmaThetat_*(nut + nu)
            );
        }

        //- Return the effective diffusivity for intermittency
        tmp<volScalarField> DgammaIntEff(const volScalarField& nut, const volScalarField& nu) const
        {
            return volScalarField::New
            (
                "DgammaIntEff",
                nut + nu
            );
        }

    //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyantMixtureKOmegaSSTRCLMRzehak&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "buoyantMixtureKOmegaSSTRCLMRzehak.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

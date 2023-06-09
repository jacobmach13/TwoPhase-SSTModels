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
    Foam::RASModels::buoyantMixtureKOmegaSSTRCRzehak

Description
    Mixture based K-Omega SST model for two phase flows with 
    Alahmadi's simplified rotation and curvature correction,
    Devolder's buoyancy production, and Rzehak's bubble induced turbulence. 

    Turbulence model described in:
    \verbatim
        Menter, F. R. & Esch, T. (2001).
        Elements of Industrial Heat Transfer Prediction.
        16th Brazilian Congress of Mechanical Engineering (COBEM).
    \endverbatim

    with updated coefficients from
    \verbatim
        Menter, F. R., Kuntz, M., and Langtry, R. (2003).
        Ten Years of Industrial Experience with the SST Turbulence Model.
        Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano,
        & M. Tummers, Begell House, Inc., 625 - 632.
    \endverbatim

    but with the consistent production terms from the 2001 paper as form in the
    2003 paper is a typo, see
    \verbatim
        http://turbmodels.larc.nasa.gov/sst.html
    \endverbatim

    and the addition of the optional F3 term for rough walls from
    \verbatim
        Hellsten, A. (1998).
        "Some Improvements in Menter’s k-omega-SST turbulence model"
        29th AIAA Fluid Dynamics Conference, AIAA-98-2554.
    \endverbatim

    Mixture properties are adapted from:
    \verbatim
        Behzadi, A., Issa, R. I., & Rusche, H. (2004).
        Modelling of dispersed bubble and droplet flow at high phase fractions.
        Chemical Engineering Science, 59(4), 759-770.
    \endverbatim

    Rzehak and Krepper's bubble induced turbulence is described in:
    \verbatim
        Rzehak, R., Krepper, E., 2013. Cfd modeling of bubble-induced turbulence.
        Int. J. Multiph. Flow 55, 138–155.
    \endverbatim

    Alahmadi's simplified rotation/curvature correction is described in:
    \verbatim
        Yaser H. Alahmadi, Sawsan A. Awadh, and Andrzej F. Nowakowski. 2021.
        Simulation of Swirling Flow with a Vortex Breakdown Using Modified Shear Stress Transport Model.
        Industrial Engineering & Chemistry Research 60(16), pp. 6016 - 6026
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficuients in a consistent
    manner.  The paper suggests that sigma is blended but this would not be
    consistent with the blending of the k-epsilon and k-omega models.

    Also note that the error in the last term of equation (2) relating to
    sigma has been corrected.

    Wall-functions are applied in this implementation by using equations (14)
    to specify the near-wall omega as appropriate.

    The blending functions (15) and (16) are not currently used because of the
    uncertainty in their origin, range of applicability and that if y+ becomes
    sufficiently small blending u_tau in this manner clearly becomes nonsense.

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
            // Rotation and curvature correction
            crc         1.4;
            //Buoyancy production term
            sigmaT      0.85;
            // Bubble generated turbulence
            Cp          1.0; 
        }
    \endverbatim

SourceFiles
    buoyantMixtureKOmegaSSTRCRzehak.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantMixtureKOmegaSSTRCRzehak_H
#define buoyantMixtureKOmegaSSTRCRzehak_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseSystem.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                 Class buoyantMixtureKOmegaSSTRCRzehak Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class buoyantMixtureKOmegaSSTRCRzehak
:
public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
    // Private Data
        mutable buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>
            *liquidTurbulencePtr_;
        
    // Private Member Functions
        //- Return the turbulence model for the other phase
        buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>& liquidTurbulence() const;
protected:
    // Protected data
    
        // Default kOmegaSST coefficients

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;

            Switch F3_;

        // Rotation and Curvature correction coefficients
            dimensionedScalar cr1_;
            dimensionedScalar cr2_;
            dimensionedScalar cr3_;

        // Devolder's buoyancy production coefficient
            dimensionedScalar sigmaT_;
        
        // Rzehak and Krepper(2013) bubble induced turbulence
            dimensionedScalar Cp_;
            

        // Fields

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            volScalarField k_;
            volScalarField omega_;

            // Mixture fields
            autoPtr<volScalarField> Ct2_;
            autoPtr<volScalarField> rhom_;
            autoPtr<volScalarField> km_;
            autoPtr<volScalarField> omegam_;
        

    // Protected Member Functions
        
        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega, const volScalarField& kp, const volScalarField& omegap, const volScalarField& nup) const;
        virtual tmp<volScalarField> F2(const volScalarField& kp, const volScalarField& omegap, const volScalarField& nup) const;
        virtual tmp<volScalarField> F3(const volScalarField& omegap, const volScalarField& nup) const;
        virtual tmp<volScalarField> F23(const volScalarField& kp, const volScalarField& omegap, const volScalarField& nup) const;
        
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

       tmp<volScalarField> alphaK(const volScalarField& F1) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega(const volScalarField& F1) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField> beta
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, beta1_, beta2_);
        }

        tmp<volScalarField> gamma
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, gamma1_, gamma2_);
        }

        virtual void correctNut
        (
            const volScalarField& S2,
            const volScalarField& F2
        );
        
        virtual void correctNut();

        //- Return k production rate
        virtual tmp<volScalarField> Pk
        (
            const volScalarField& G
        ) const;

        //- Return epsilon/k which for standard RAS is betaStar*omega
        virtual tmp<volScalarField> epsilonByk
        (
            const volScalarField& F1,
            const volScalarField& F2
        ) const;

        // Mixture 

        wordList omegaBoundaryTypes(const volScalarField& omega) const;

        void correctInletOutlet
        (
            volScalarField& vsf,
            const volScalarField& refVsf
        ) const;

        void initMixtureFields();

        tmp<volScalarField> Ct2() const;

        tmp<volScalarField> rholEff() const;
        tmp<volScalarField> rhogEff() const;
        tmp<volScalarField> rhom() const;

        tmp<volScalarField> mix
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volScalarField> mixU
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<surfaceScalarField> mixFlux
        (
            const surfaceScalarField& fc,
            const surfaceScalarField& fd
        ) const;

        tmp<volScalarField> bubbleG() const;

        virtual tmp<fvScalarMatrix> kSource() const;

        virtual tmp<fvScalarMatrix> omegaSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("buoyantMixtureKOmegaSSTRCRzehak");


    // Constructors

        //- Construct from components
        buoyantMixtureKOmegaSSTRCRzehak
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
        buoyantMixtureKOmegaSSTRCRzehak(const buoyantMixtureKOmegaSSTRCRzehak&) = delete;


    //- Destructor
    virtual ~buoyantMixtureKOmegaSSTRCRzehak()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1, const volScalarField& nut, const volScalarField& nu) const
        {
            return volScalarField::New
            (
                "DkEff",
                alphaK(F1)*nut + nu
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1, const volScalarField& nut, const volScalarField& nu) const
        {
            return volScalarField::New
            (
                "DomegaEff",
                alphaOmega(F1)*nut + nu
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return volScalarField::New
            (
                "epsilon",
                betaStar_*k_*omega_,
                omega_.boundaryField().types()
            );
        }

        //- Return the specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyantMixtureKOmegaSSTRCRzehak&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "buoyantMixtureKOmegaSSTRCRzehak.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

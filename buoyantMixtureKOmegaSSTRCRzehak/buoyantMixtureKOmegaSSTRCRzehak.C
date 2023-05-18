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

\*---------------------------------------------------------------------------*/

#include "buoyantMixtureKOmegaSSTRCRzehak.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
#include "phaseSystem.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "fixedValueFvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "fvmSup.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCRzehak::F1
(
    const volScalarField& CDkOmega,
    const volScalarField& kp,
    const volScalarField& omegap,
    const volScalarField& nup
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(kp)/(omegap*y_),
                scalar(500)*(nup)/(sqr(y_)*omegap)
            ),
            (4*alphaOmega2_)*kp/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCRzehak::F2
(
    const volScalarField& kp,
    const volScalarField& omegap,
    const volScalarField& nup
) const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(kp)/(omegap*y_),
            scalar(500)*(nup)/(sqr(y_)*omegap)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCRzehak::F3
(
    const volScalarField& omegap,
    const volScalarField& nup
) const
{
    tmp<volScalarField> arg3 = min
    (
        150*(nup)/(omegap*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCRzehak::F23
(
    const volScalarField& kp, 
    const volScalarField& omegap,
    const volScalarField& nup
) const
{

    tmp<volScalarField> f23(F2(kp,omegap,nup));
    
    if (F3_)
    {
        f23.ref() *= F3(omegap,nup);
    }


    return f23;
}

template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}

template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23(this->k_,this->omega_,this->mu()/this->rho_));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::Pk
(
    const volScalarField& G
) const
{
    return min(G, (c1_*betaStar_)*this->km_()*this->omegam_());
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField& F1,
    const volScalarField& F2
) const
{
    return betaStar_*this->omegam_();
}

// Mixture

template<class BasicMomentumTransportModel>
wordList buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::omegaBoundaryTypes
(
    const volScalarField& omega
) const
{
    const volScalarField::Boundary& obf = omega.boundaryField();

    wordList obt = obf.types();

    forAll(obf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(obf[patchi]))
        {
            obt[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    return obt;
}


template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::correctInletOutlet
(
    volScalarField& vsf,
    const volScalarField& refVsf
) const
{
    volScalarField::Boundary& bf = vsf.boundaryFieldRef();
    const volScalarField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchScalarField>(bf[patchi])
         && isA<inletOutletFvPatchScalarField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchScalarField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchScalarField>
            (refBf[patchi]).refValue();
        }
    }
}

template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::initMixtureFields()
{
    if (rhom_.valid()) return;

    // Local references to gas-phase properties
    const volScalarField& kg = this->k_;
    const volScalarField& omegag = this->omega_;
    

    // Local references to liquid-phase properties
    buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>& turbc = this->liquidTurbulence();
    const volScalarField& kl = turbc.k_;
    const volScalarField& omegal = turbc.omega_;

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    Ct2_.set
    (
        new volScalarField
        (
            IOobject
            (
                "Ct2",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Ct2()
        )
    );

    rhom_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rhom",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rhom()
        )
    );

    km_.set
    (
        new volScalarField
        (
            IOobject
            (
                "km",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(kl, kg),
            kl.boundaryField().types()
        )
    );
    correctInletOutlet(km_(), kl);

    omegam_.set
    (
        new volScalarField
        (
            IOobject
            (
                "omegam",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(omegal, omegag),
            omegaBoundaryTypes(omegal)
        )
    );
    correctInletOutlet(omegam_(), omegal);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCRzehak
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    liquidTurbulencePtr_(nullptr),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
    cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr1",
            this->coeffDict_,
            1.0
        )
    ),
    cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr2",
            this->coeffDict_,
            2.0
        )
    ),
    cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr3",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaT",
            this->coeffDict_,
            0.85
        )
    ),
    Cp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            1.0
        )
    ),
    y_(wallDist::New(this->mesh_).y()),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
    
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        cr1_.readIfPresent(this->coeffDict());
        cr2_.readIfPresent(this->coeffDict());
        cr3_.readIfPresent(this->coeffDict());
        sigmaT_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>&
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const phaseSystem& fluid = gas.fluid();
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &const_cast<buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>&>
            (
                U.db().lookupObject
                <
                    buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>
                >
                (
                    IOobject::groupName
                    (
                        momentumTransportModel::typeName,
                        liquid.name()
                    )
                )
            );
    }

    return *liquidTurbulencePtr_;
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::Ct2() const
{
    const buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>& liquidTurbulence = this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    const transportModel& liquid = fluid.otherPhase(gas);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    const volScalarField& alphag = this->alpha_;

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    volScalarField beta
    (
        (6*betaStar_/(4*sqrt(3.0/2.0)))
       *drag.K()/(liquid.rho()*liquidTurbulence.omega_)
    );
    volScalarField Ct0((3 + beta)/(1 + beta + 2*gas.rho()/liquid.rho()));
    volScalarField fAlphad((180 + (-4.71e3 + 4.26e4*alphag)*alphag)*alphag);

    return sqr(1 + (Ct0 - 1)*exp(-fAlphad));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::rholEff() const
{
    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    return fluid.otherPhase(gas).rho();
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::rhogEff() const
{
    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    const virtualMassModel& virtualMass = fluid.lookupSubModel<virtualMassModel>(gas, fluid.otherPhase(gas));
    return gas.rho() + virtualMass.Cvm()*fluid.otherPhase(gas).rho();
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::rhom() const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return alphal*rholEff() + alphag*rhogEff();
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::mix
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::mixU
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}

template<class BasicMomentumTransportModel>
tmp<surfaceScalarField> buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::mixFlux
(
    const surfaceScalarField& fc,
    const surfaceScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    surfaceScalarField alphalf(fvc::interpolate(alphal));
    surfaceScalarField alphagf(fvc::interpolate(alphag));

    surfaceScalarField rholEfff(fvc::interpolate(rholEff()));
    surfaceScalarField rhogEfff(fvc::interpolate(rhogEff()));

    return
       (alphalf*rholEfff*fc + alphagf*rhogEfff*fvc::interpolate(Ct2_())*fd)
      /(alphalf*rholEfff + alphagf*rhogEfff*fvc::interpolate(Ct2_()));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::bubbleG() const
{
    const buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>& liquidTurbulence = this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    const transportModel& liquid = fluid.otherPhase(gas);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    // Rzehak and Krepper(2013)
    tmp<volScalarField> bubbleG
    (
        Cp_*liquid*drag.K()*sqr(magUr)
    );

    return bubbleG;
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::kSource() const
{
    return fvm::Su(bubbleG(), km_());
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::omegaSource() const
{
    const buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>& liquidTurbulence = this->liquidTurbulence();
    const volScalarField& kl = liquidTurbulence.k_;
    const volScalarField& omegal = liquidTurbulence.omega_;

    // time scale
    const transportModel& gas = this->transport();
    volScalarField  tau = gas.d()/sqrt(kl);
    
    //Rzehak and Krepper(2013) model
    tmp<volScalarField> SLOmega
    (
        scalar(1)/(betaStar_*kl) * (bubbleG()/tau) - (omegal/kl)*bubbleG()
    );

    return fvm::Su(SLOmega, omegam_());
}

template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::correct()
{
    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();

    // Only solve the mixture turbulence for the gas-phase
    if (&gas != &fluid.phases()[0])
    {
        // This is the liquid phase but check the model for the gas-phase
        // is consistent
        this->liquidTurbulence();

        return;
    }

    if (!this->turbulence_)
    {
        return;
    }

    // Initialise the mixture fields if they have not yet been constructed
    initMixtureFields();

    // Local references to gas-phase properties
    tmp<surfaceScalarField> phig = this->phi();
    const volVectorField& Ug = this->U_;
    const volScalarField& alphag = this->alpha_;
    volScalarField& kg = this->k_;
    volScalarField& omegag = this->omega_;
    volScalarField& nutg = this->nut_;

    // Local references to liquid-phase properties
    buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>& liquidTurbulence = this->liquidTurbulence();
    tmp<surfaceScalarField> phil = liquidTurbulence.phi();
    const volVectorField& Ul = liquidTurbulence.U_;
    const volScalarField& alphal = liquidTurbulence.alpha_;
    volScalarField& kl = liquidTurbulence.k_;
    volScalarField& omegal = liquidTurbulence.omega_;
    volScalarField& nutl = liquidTurbulence.nut_;

    // Local references to mixture properties
    volScalarField& rhom = rhom_();
    volScalarField& km = km_();
    volScalarField& omegam = omegam_();

    fv::options& fvOptions(fv::options::New(this->mesh_));
    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    // Update the effective mixture density
    rhom = this->rhom();

    // Mixture volumetric flux
    surfaceScalarField phim("phim", mixFlux(phil, phig));

    // Mixture mass flux
    surfaceScalarField rhoPhim = fvc::interpolate(rhom) * phim;

    // Mixture velocity divergence
    volScalarField divUm
    (
        mixU
        (
            fvc::div(fvc::absolute(phil, Ul)),
            fvc::div(fvc::absolute(phig, Ug))
        )
    );

    // Individual phase strain rates
    tmp<volTensorField> tgradUl = fvc::grad(Ul);
    volScalarField S2c(2*magSqr(symm(tgradUl())));

    tmp<volTensorField> tgradUg = fvc::grad(Ug);
    volScalarField S2d(2*magSqr(symm(tgradUg())));

    // Mixture strain rate
    volScalarField S2m(mixU(S2c,S2d));

    // Liquid phase GbyNu
    tmp<volScalarField> GbyNuc;
    {
        tmp<volTensorField> tgradUl = fvc::grad(Ul);
        GbyNuc = tmp<volScalarField>
        (
            new volScalarField
            (
                tgradUl() && dev(twoSymm(tgradUl()))
            )
        );
        tgradUl.clear();
    }

    // Gas phase GbyNu
    tmp<volScalarField> GbyNud;
    {
        tmp<volTensorField> tgradUg = fvc::grad(Ug);
        GbyNud = tmp<volScalarField>
        (
            new volScalarField
            (
                tgradUg() && dev(twoSymm(tgradUg()))
            )
        );
        tgradUg.clear();
    }

    // Mixture GbyNu
    volScalarField GbyNum(mix(GbyNuc, GbyNud));

    // Liquid phase turbulence production
    tmp<volScalarField> Gc;
    {
        tmp<volTensorField> tgradUl = fvc::grad(Ul);
        Gc = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutl*(tgradUl() && dev(twoSymm(tgradUl())))
            )
        );
        
        
        tgradUl.clear();

        // Update k, omega and G at the wall
        kl.boundaryFieldRef().updateCoeffs();
        omegal.boundaryFieldRef().updateCoeffs();

        Gc.ref().checkOut();
    }

    // Gas phase turbulence production
    tmp<volScalarField> Gd;
    {
        tmp<volTensorField> tgradUg = fvc::grad(Ug);
        Gd = tmp<volScalarField>
        (
            new volScalarField
            (
                this->GName(),
                nutg*(tgradUg() && dev(twoSymm(tgradUg())))
            )
        );
        
        tgradUg.clear();

        // Update k, omega and G at the wall
        kg.boundaryFieldRef().updateCoeffs();
        omegag.boundaryFieldRef().updateCoeffs();

        Gd.ref().checkOut();
    }

    // Mixture turbulence production
    volScalarField Gm(mix(Gc, Gd));

    // Mixture turbulence viscosity
    volScalarField nutm(mixU(nutl, nutg));

    // Update the mixture k and omega boundary conditions
    km == mix(kl, kg);
    bound(km, this->kMin_);

    omegam == mix(omegal, omegag);
    bound(omegam, this->omegaMin_);

    // Simplified rotation and curvature correction from Alahmadi et al. (2021)
    volScalarField Sl = sqrt(2*magSqr(symm(fvc::grad(Ul)))); // strain rate for the liquid phase
    volScalarField Sg = sqrt(2*magSqr(symm(fvc::grad(Ug)))); // strain rate for the gas phase
    volScalarField Sm(mix(Sl,Sg));                           // mixture strain rate

    volScalarField Wl = sqrt(2*magSqr(skew(fvc::grad(Ul)))); // liquid phase vorticity
    volScalarField Wg = sqrt(2*magSqr(skew(fvc::grad(Ug)))); // gas phase vorticity
    volScalarField Wm(mix(Wl,Wg));                           // mixture vorticity

    volScalarField Ri = Wm/Sm * (Wm/Sm - 1.0);              // Mixture Richardson number
    volScalarField rStar = Sm/max(Wm, dimensionedScalar("minWm", Wm.dimensions(), 1e-6));   // Mixture r*               
    volScalarField rTilda = Ri*rStar;                       // Mixture rTilda

    volScalarField Fr1 // to be multiplied with the k production term
    (
	max
    	(
		min
		(
        		(scalar(1.0)+cr1_)*2.0*rStar/(scalar(1.0)+rStar)*(scalar(1.0) -cr3_*atan(cr2_*rTilda))-cr1_,
        		scalar(1.25)
		),
		scalar(0.0)
	)
    );
    
    //Devolder et al. (2018) buoyancy production term
    const uniformDimensionedVectorField& g = this->mesh_.objectRegistry::template lookupObject<uniformDimensionedVectorField>("g"); // lookup gravitional acceleration

    volScalarField Gb("Gb", -nutm/sigmaT_*(g & fvc::grad(rhom))); 

    // Gas and liquid phase molecular viscosity
    const transportModel& liquid = fluid.otherPhase(gas);
    volScalarField  nug = gas.thermo().nu();
    volScalarField nul = liquid.thermo().nu();

    // mixture molecular viscosity
    volScalarField num(mixU(nul,nug));

    // Gas phase CDKOmega
    volScalarField CDkOmegag
    (
        (2*alphaOmega2_)*(fvc::grad(kg) & fvc::grad(omegag))/omegag
    );

    // Liquid phase CDKOmega
    volScalarField CDkOmegal
    (
        (2*alphaOmega2_)*(fvc::grad(kl) & fvc::grad(omegal))/omegal
    );

    // Mixture CDkOmega
    volScalarField CDkOmegam(mix(CDkOmegal, CDkOmegag));

    // Gas and liquid phase F1
    volScalarField F1g(this->F1(CDkOmegag, kg, omegag, nug));
    volScalarField F1l(this->F1(CDkOmegal, kl, omegal, nul));

    // Mixture F1
    volScalarField F1m(mix(F1l, F1g));

    // Gas and liquid phase F23
    volScalarField F23g(this->F23(kg, omegag, nug));
    volScalarField F23l(this->F23(kl, omegal, nul));

    // Mixture F23
    volScalarField F23(mix(F23l,F23g));

    {
        // Gas and liquid phase gamma
        volScalarField gammag(this->gamma(F1g));
        volScalarField gammal(this->gamma(F1l));

        // Mixture gamma
        volScalarField gammam(mix(gammal, gammag));

        // Gas and liquid phase beta
        volScalarField betag(this->beta(F1g));
        volScalarField betal(this->beta(F1l));

        // Mixture beta
        volScalarField betam(mix(betal, betag));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(rhom,omegam)
          + fvm::div(rhoPhim, omegam)
          - fvm::laplacian(rhom*DomegaEff(F1m, nutm, num), omegam)
         ==
            rhom*gammam
           *min
            (
                GbyNum,
                (c1_/a1_)*betaStar_*omegam
               *max(a1_*omegam, b1_*F23*sqrt(S2m))
            )
          - fvm::SuSp((2.0/3.0)*rhom*gammam*divUm, omegam)
          - fvm::Sp(rhom*betam*omegam, omegam)
          - fvm::SuSp
            (
                rhom*(F1m - scalar(1))*CDkOmegam/omegam,
                omegam
            )
          + omegaSource()	// Rzehak & Krepper bubble induced turbulence
          + fvOptions(rhom,omegam)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omegam.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omegam);
        bound(omegam, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rhom, km)
      + fvm::div(rhoPhim, km)
      - fvm::laplacian(rhom*DkEff(F1m, nutm, num), km)
     ==
        rhom*Fr1*this->Pk(Gm)		     // Fr1 is rotation/curvature correction
      + fvm::Sp(Gb/max(km, this->kMin_), km) //buoyancy correction
      - fvm::SuSp((2.0/3.0)*rhom*divUm, km)
      - fvm::Sp(rhom*this->epsilonByk(F1m,F23), km)
      + kSource()			     // Rzehak & Krepper bubble induced turbulence (k production)
      + fvOptions(rhom,km)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(km);
    bound(km, this->kMin_);
    km.correctBoundaryConditions();

    volScalarField Cc2(rhom/(alphal*rholEff() + alphag*rhogEff()*Ct2_()));

    kl = Cc2*km;
    kl.correctBoundaryConditions();

    omegal = Cc2*omegam;
    omegal.correctBoundaryConditions();

    liquidTurbulence.correctNut();

    Ct2_() = Ct2();

    kg = Ct2_()*kl;
    kg.correctBoundaryConditions();

    omegag = Ct2_()*omegal;
    omegag.correctBoundaryConditions();

    nutg = Ct2_()*(liquidTurbulence.nu()/this->nu())*nutl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

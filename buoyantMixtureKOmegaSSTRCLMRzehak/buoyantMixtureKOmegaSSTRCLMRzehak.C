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

#include "buoyantMixtureKOmegaSSTRCLMRzehak.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCLMRzehak::F1
(
    const volScalarField& CDkOmega,
    const volScalarField& kp,
    const volScalarField& omegap,
    const volScalarField& nup
) const
{
    volScalarField y = wallDist::New(this->mesh_).y();

    volScalarField Ry(y*sqrt(kp)/nup);
    volScalarField F3(exp(-pow(Ry/120.0, 8)));

    return max(buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::F1(CDkOmega, kp, omegap, nup), F3);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::Pk
(
    const volScalarField& G
) const
{

    return this->gammaIntEffm_()*buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::Pk(G);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField& F1,
    const volScalarField& F2
) const
{
    
    return
       min(max(this->gammaIntEffm_(), scalar(0.1)), scalar(1))
       *buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::epsilonByk(F1, F2);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::Fthetat
(
    const volScalarField& Usp,
    const volScalarField& Omegap,
    const volScalarField& omegap,
    const volScalarField& nup,
    const volScalarField& gammaIntp,
    const volScalarField& ReThetatp
) const
{
    volScalarField y = wallDist::New(this->mesh_).y();

    volScalarField delta = scalar(375)*Omegap*nup*ReThetatp*y/sqr(Usp);

    volScalarField ReOmega = (sqr(y)*omegap/nup);

    volScalarField Fwake = exp(-sqr(ReOmega/1e5));

    volScalarField dummyE = exp(-pow4((y/delta)));

    return volScalarField::New
    (
        IOobject::groupName("Fthetat", this->alphaRhoPhi_.group()),
        min
        (
            max
            (
                Fwake*dummyE,
                (1 - sqr((gammaIntp - 1.0/ce2_)/(1 - 1.0/ce2_)))
            ),
            scalar(1)
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::ReThetac() const
{

    tmp<volScalarField> tReThetac
    (
        volScalarField::New
        (
            IOobject::groupName("ReThetac", this->alphaRhoPhi_.group()),
            this->mesh_,
            dimless
        )
    );
    volScalarField ReThetac = tReThetac.ref();

    volScalarField ReThetatp = this->ReThetatm_();

    forAll(ReThetac, celli)
    {
        scalar ReThetat = ReThetatp[celli];

        ReThetac[celli] =
            ReThetat <= 1870
          ?
            ReThetat
          - 396.035e-2
          + 120.656e-4*ReThetat
          - 868.230e-6*sqr(ReThetat)
          + 696.506e-9*pow3(ReThetat)
          - 174.105e-12*pow4(ReThetat)
          :
            ReThetat - 593.11 - 0.482*(ReThetat - 1870);
    }

    return tReThetac;
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::Flength
(
    const volScalarField& nup
) const
{

    tmp<volScalarField> tFlength
    (
        volScalarField::New
        (
            IOobject::groupName("Flength", this->alphaRhoPhi_.group()),
            this->mesh_,
            dimless
        )
    );
    volScalarField& Flength = tFlength.ref();

    volScalarField y = wallDist::New(this->mesh_).y();
    volScalarField ReThetatp = this->ReThetatm_();
    volScalarField omegap = this->omegam_();

    forAll(ReThetatp, celli)
    {
        scalar ReThetat = ReThetatp[celli];

        if (ReThetat < 400)
        {
            Flength[celli] =
                398.189e-1
              - 119.270e-4*ReThetat
              - 132.567e-6*sqr(ReThetat);
        }
        else if (ReThetat < 596)
        {
            Flength[celli] =
                263.404
              - 123.939e-2*ReThetat
              + 194.548e-5*sqr(ReThetat)
              - 101.695e-8*pow3(ReThetat);
        }
        else if (ReThetat < 1200)
        {
            Flength[celli] = 0.5 - 3e-4*(ReThetat - 596);
        }
        else
        {
            Flength[celli] = 0.3188;
        }

        scalar Fsublayer = exp(-sqr(sqr(y[celli])*omegap[celli]/(200*nup[celli])));

        Flength[celli] = Flength[celli]*(1 - Fsublayer) + 40*Fsublayer;
    }

    return tFlength;
}

template<class BasicMomentumTransportModel>
tmp<volScalarField>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::ReThetat0
(
    const volScalarField& kp,
    const volScalarField& Usp,
    const volScalarField& dUsdsp,
    const volScalarField& nup
) const
{
    tmp<volScalarField> tReThetat0
    (
        volScalarField::New
        (
            IOobject::groupName("ReThetat0", this->alphaRhoPhi_.group()),
            this->mesh_,
            dimless
        )
    );
    volScalarField ReThetat0 = tReThetat0.ref();


    label maxIter = 0;

    forAll(ReThetat0, celli)
    {
        const scalar Tu
        (
            max(100*sqrt((2.0/3.0)*kp[celli])/Usp[celli], scalar(0.027))
        );

        // Initialize lambda to zero.
        // If lambda were cached between time-steps convergence would be faster
        // starting from the previous time-step value.
        scalar lambda = 0;

        scalar lambdaErr;
        scalar thetat;
        label iter = 0;

        do
        {
            // Previous iteration lambda for convergence test
            const scalar lambda0 = lambda;

            if (Tu <= 1.3)
            {
                const scalar Flambda =
                    dUsdsp[celli] <= 0
                  ?
                    1
                  - (
                     - 12.986*lambda
                     - 123.66*sqr(lambda)
                     - 405.689*pow3(lambda)
                    )*exp(-pow(Tu/1.5, 1.5))
                  :
                    1
                  + 0.275*(1 - exp(-35*lambda))
                   *exp(-Tu/0.5);

                thetat =
                    (1173.51 - 589.428*Tu + 0.2196/sqr(Tu))
                   *Flambda*nup[celli]
                   /Usp[celli];
            }
            else
            {
                const scalar Flambda =
                    dUsdsp[celli] <= 0
                  ?
                    1
                  - (
                      -12.986*lambda
                      -123.66*sqr(lambda)
                      -405.689*pow3(lambda)
                    )*exp(-pow(Tu/1.5, 1.5))
                  :
                    1
                  + 0.275*(1 - exp(-35*lambda))
                   *exp(-2*Tu);

                thetat =
                    331.50*pow((Tu - 0.5658), -0.671)
                   *Flambda*nup[celli]/Usp[celli];
            }

            lambda = sqr(thetat)/nup[celli]*dUsdsp[celli];
            lambda = max(min(lambda, 0.1), -0.1);

            lambdaErr = mag(lambda - lambda0);

            maxIter = max(maxIter, ++iter);

        } while (lambdaErr > lambdaErr_);

        ReThetat0[celli] = max(thetat*Usp[celli]/nup[celli], scalar(20));
    }

    if (maxIter > maxLambdaIter_)
    {
        WarningInFunction
            << "Number of lambda iterations exceeds maxLambdaIter("
            << maxLambdaIter_ << ')'<< endl;
    }
    return tReThetat0;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::Fonset
(
    const volScalarField& Revp,
    const volScalarField& RTp
) const
{

    volScalarField Fonset1(Revp/(2.193*max(this->ReThetac(), scalar(1e-6))));

    volScalarField Fonset2
    (
        min(max(Fonset1, pow4(Fonset1)), scalar(2))
    );

    volScalarField Fonset3(max(1 - pow3(RTp/2.5), scalar(0)));

    return volScalarField::New
    (
        IOobject::groupName("Fonset", this->alphaRhoPhi_.group()),
        max(Fonset2 - Fonset3, scalar(0))
    );
}




template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::initMixtureFields()
{
    if (ReThetatm_.valid()) return;

    // initialize fields specified in parent class
    buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::initMixtureFields();

    // Local references to gas-phase properties
    const volScalarField& ReThetatg = this->ReThetat_;
    const volScalarField& gammaIntg = this->gammaInt_;

    // Local references to liquid-phase properties
    buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>& turbc = this->liquidTurbulence();
    const volScalarField& ReThetatl = turbc.ReThetat_;
    const volScalarField& gammaIntl = turbc.gammaInt_;

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    ReThetatm_.set
    (
        new volScalarField
        (
            IOobject
            (
                "ReThetatm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mix(ReThetatl, ReThetatg),
            ReThetatl.boundaryField().types()
        )
    );
    this->correctInletOutlet(ReThetatm_(), ReThetatl);

    gammaIntm_.set
    (
        new volScalarField
        (
            IOobject
            (
                "gammaIntm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mix(gammaIntl, gammaIntg),
            gammaIntl.boundaryField().types()
        )
    );
    this->correctInletOutlet(gammaIntm_(), gammaIntl);

    gammaIntEffm_.set
    (
        new volScalarField
        (
            IOobject
            (
                "gammaIntm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh_,
            dimensionedScalar(dimless, 0)
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::buoyantMixtureKOmegaSSTRCLMRzehak
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
    buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    liquidTurbulencePtr_(nullptr),
    ca1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ca1",
            this->coeffDict_,
            2
        )
    ),
    ca2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ce1",
            this->coeffDict_,
            1
        )
    ),
    ce2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50
        )
    ),
    cThetat_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cThetat",
            this->coeffDict_,
            0.03
        )
    ),
    sigmaThetat_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaThetat",
            this->coeffDict_,
            2
        )
    ),
    lambdaErr_
    (
        this->coeffDict_.lookupOrDefault("lambdaErr", 1e-6)
    ),
    maxLambdaIter_
    (
        this->coeffDict_.lookupOrDefault("maxLambdaIter", 10)
    ),
    deltaU_("deltaU", dimVelocity, small),
    ReThetat_
    (
        IOobject
        (
            IOobject::groupName("ReThetat", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    gammaInt_
    (
        IOobject
        (
            IOobject::groupName("gammaInt", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::read()
{
    if (buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::read())
    {
        ca1_.readIfPresent(this->coeffDict());
        ca2_.readIfPresent(this->coeffDict());
        ce1_.readIfPresent(this->coeffDict());
        ce2_.readIfPresent(this->coeffDict());
        sigmaThetat_.readIfPresent(this->coeffDict());
        cThetat_.readIfPresent(this->coeffDict());
        this->coeffDict().readIfPresent("lambdaErr", lambdaErr_);
        this->coeffDict().readIfPresent("maxLambdaIter", maxLambdaIter_);

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>&
buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const phaseSystem& fluid = gas.fluid();
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &const_cast<buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>&>
            (
                U.db().lookupObject
                <
                    buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>
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
void buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::correctReThetatGammaInt()
{
    // Initialise the mixture fields if they have not yet been constructed
    initMixtureFields();

    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();

    // Local references to gas-phase properties
    tmp<surfaceScalarField> phig = this->phi();
    const volVectorField& Ug = this->U_;
    const volScalarField& alphag = this->alpha_;
    volScalarField& kg = this->k_;
    volScalarField& omegag = this->omega_;
    volScalarField& nutg = this->nut_;
    volScalarField& gammaIntg = this->gammaInt_;
    volScalarField& ReThetatg = this->ReThetat_;

    tmp<volScalarField> tnug = this->nu();
    volScalarField nuNg = tnug.ref();

    // Local references to liquid-phase properties
    buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>& liquidTurbulence = this->liquidTurbulence();
    tmp<surfaceScalarField> phil = liquidTurbulence.phi();
    const volVectorField& Ul = liquidTurbulence.U_;
    const volScalarField& alphal = liquidTurbulence.alpha_;
    volScalarField& kl = liquidTurbulence.k_;
    volScalarField& omegal = liquidTurbulence.omega_;
    volScalarField& nutl = liquidTurbulence.nut_;
    volScalarField& gammaIntl = liquidTurbulence.gammaInt_;
    volScalarField& ReThetatl = liquidTurbulence.ReThetat_;

    tmp<volScalarField> tnul = liquidTurbulence.nu();
    volScalarField nuNl = tnul.ref();
    
    // Local references to mixture properties
    volScalarField& rhom = this->rhom_();
    volScalarField& km = this->km_();
    volScalarField& omegam = this->omegam_();

    volScalarField& ReThetatm = ReThetatm_();
    volScalarField& gammaIntm = gammaIntm_();
    volScalarField& gammaIntEffm = gammaIntEffm_();

    // Distance to the nearest wall
    volScalarField y = wallDist::New(this->mesh_).y();

    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Update the effective mixture density
    rhom = this->rhom();

    // Mixture volumetric flux
    surfaceScalarField phim("phim", this->mixFlux(phil, phig));

    // Mixture mass flux
    surfaceScalarField rhoPhim = fvc::interpolate(rhom) * phim;

    //Mixture nuN
    volScalarField nuNm = this->mix(nuNl,nuNg);

    // Gas and liquid phase molecular viscosity
    const transportModel& liquid = fluid.otherPhase(gas);
    volScalarField  nug = gas.thermo().nu();
    volScalarField  nul = liquid.thermo().nu();

    // mixture molecular viscosity
    volScalarField num(this->mixU(nul,nug));

    // Mixture turbulence viscosity
    volScalarField nutm(this->mixU(nutl, nutg));

    // Update the mixture k and omega boundary conditions
    km == this->mix(kl, kg);
    bound(km, this->kMin_);

    omegam == this->mix(omegal, omegag);
    bound(omegam,this->omegaMin_);

    // Update the mixture gammaInt and ReThetat boundary conditions
    gammaIntm == this->mix(gammaIntl, gammaIntg);
    bound(gammaIntm, 0);

    ReThetatm == this->mix(ReThetatl, ReThetatg);
    bound(ReThetatm, 1e-6);

    // Langtry-Menter Fields derived from the velocity gradient (gas phase)
    tmp<volTensorField> tgradUg = fvc::grad(Ug);
    volScalarField S2d(2*magSqr(symm(tgradUg())));
    volScalarField Omegag = sqrt(2*magSqr(skew(tgradUg())));
    volScalarField Usg = max(mag(Ug), deltaU_);
    volScalarField dUsdsg = (Ug & (Ug & tgradUg()))/sqr(Usg);
    tgradUg.clear();

    // Langtry-Menter Fields derived from the velocity gradient (liquid phase)
    tmp<volTensorField> tgradUl = fvc::grad(Ul);
    volScalarField S2c(2*magSqr(symm(tgradUl())));
    volScalarField Omegal = sqrt(2*magSqr(skew(tgradUl())));
    volScalarField Usl = max(mag(Ul), deltaU_);
    volScalarField dUsdsl = (Ul & (Ul & tgradUl()))/sqr(Usl);
    tgradUl.clear();


    // Langtry-Menter Fields derived from velocity gradient (Mixture)
    volScalarField Omegam = this->mix(Omegal,Omegag);
    volScalarField Usm = this->mixU(Usl,Usg);
    volScalarField dUsdsm = this->mix(dUsdsl, dUsdsg);
    volScalarField S2m = this->mix(S2c,S2d);

    // Mixture Fthetat
    volScalarField Fthetatm = this->Fthetat(Usm, Omegam, omegam, nuNm, gammaIntm, ReThetatm);

    {
        // Mixture t
        volScalarField tm = 500*nuNm/sqr(Usm);

        // Mixture Pthetat
        volScalarField Pthetatm
        (
            rhom*(cThetat_/tm)*(1 - Fthetatm)
        );

        // Transition onset momentum-thickness Reynolds number equation
        tmp<fvScalarMatrix> ReThetatEqn
        (
            fvm::ddt(rhom, ReThetatm)
          + fvm::div(rhoPhim, ReThetatm)
          - fvm::laplacian(rhom*DReThetatEff(nutm,num), ReThetatm)
         ==
            Pthetatm*ReThetat0(km, Usm, dUsdsm, nuNm) - fvm::Sp(Pthetatm, ReThetatm)
          + fvOptions(rhom, ReThetatm)
        );

        ReThetatEqn.ref().relax();
        fvOptions.constrain(ReThetatEqn.ref());
        solve(ReThetatEqn);
        fvOptions.correct(ReThetatm);
        bound(ReThetatm, 1e-6);
    }
    
    //Mixture Rev
    volScalarField Revm = sqr(y)*sqrt(S2m)/nuNm;

    // Mixture RT
    volScalarField RTm = km/(nuNm*omegam);

    {

        volScalarField Pgamma1(sqrt(gammaIntm*this->Fonset(Revm, RTm)));

        volScalarField Pgamma2(sqrt(S2m));

        volScalarField Pgamma3(this->Flength(nuNm));

        volScalarField Pgammam
        (
            rhom*ca1_*Pgamma3*Pgamma2*Pgamma1
        );

        volScalarField Fturbm(exp(-pow4(0.25*RTm)));

        volScalarField Egammam
        (
            rhom*ca2_*Omegam*Fturbm*gammaIntm
        );

        // Intermittency equation
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(rhom, gammaIntm)
          + fvm::div(rhoPhim, gammaIntm)
          - fvm::laplacian(rhom*DgammaIntEff(nutm,num), gammaIntm)
        ==
            Pgammam - fvm::Sp(ce1_*Pgammam, gammaIntm)
          + Egammam - fvm::Sp(ce2_*Egammam, gammaIntm)
          + fvOptions(rhom, gammaIntm)
        );

        gammaIntEqn.ref().relax();
        fvOptions.constrain(gammaIntEqn.ref());
        solve(gammaIntEqn);
        fvOptions.correct(gammaIntm);
        bound(gammaIntm, 0);
    }

    volScalarField Freattach = exp(-pow4(RTm/20.0));

    volScalarField gammaSep = min(2*max(Revm/(3.235*max(this->ReThetac(),scalar(1e-6))) - 1, scalar(0))*Freattach, scalar(2))*Fthetatm;


    // Update mixture effective gammaInt
    gammaIntEffm = max(gammaIntm, gammaSep);

    
    volScalarField Cc2(rhom/(alphal*this->rholEff() + alphag*this->rhogEff()*this->Ct2()));

    gammaIntl = Cc2*gammaIntm;
    gammaIntl.correctBoundaryConditions();

    ReThetatl = Cc2*ReThetatm;
    ReThetatl.correctBoundaryConditions();


    gammaIntg = this->Ct2()*gammaIntl;
    gammaIntg.correctBoundaryConditions();

    ReThetatg = this->Ct2()*ReThetatl;
    ReThetatg.correctBoundaryConditions();

}

template<class BasicMomentumTransportModel>
void buoyantMixtureKOmegaSSTRCLMRzehak<BasicMomentumTransportModel>::correct()
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
    
    // Correct gammaInt and ReThetat
    correctReThetatGammaInt();

    // Correct k and omega
    buoyantMixtureKOmegaSSTRCRzehak<BasicMomentumTransportModel>::correct();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

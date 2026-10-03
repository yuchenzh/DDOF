/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "coefficientMultiComponentMixture.H"
#include "coefficientWilkeMultiComponentMixture.H"
#include "singleComponentMixture.H"

#include "psiThermo.H"
#include "psiReactionThermo.H"
#include "hePsiThermo.H"

#include "forGases.H"
#include "makeReactionThermo.H"
#include "ddpolynomialTransport.H"

#include "ode.H"
#include "chemistryModel.H"
 
#include "makeChemistrySolver.H"

#include "ode.H"
#include "chemistryModel.H"

#include "forLiquids.H"
#include "makeChemistryReductionMethod.H"
#include "makeReaction.H"
#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "FallOffReactionRate.H"

#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "LangmuirHinshelwoodReactionRate.H"

#include "MichaelisMentenReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePsiReactionThermos(Mixture, ThermoPhysics)                         \
    makeReactionThermos                                                        \
    (                                                                          \
        psiThermo,                                                             \
        psiReactionThermo,                                                     \
        hePsiThermo,                                                           \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

#define makePsiReactionThermo(Mixture, ThermoPhysics)                          \
    makeReactionThermo                                                         \
    (                                                                          \
        psiReactionThermo,                                                     \
        hePsiThermo,                                                           \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forGasEnergiesAndThermos(ddpolynomialTransport, makePsiReactionThermos, singleComponentMixture);

    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makePsiReactionThermos, coefficientMultiComponentMixture);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, defineChemistrySolvers, nullArg);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeChemistrySolvers, ode);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, defineChemistryReductionMethod, nullArg);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeChemistryReductionMethod, none);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, defineReaction, nullArg);


    // Irreversible/reversible/non-equilibrium-reversible reactions

    //forCoeffGases(makeIRNReactions, ArrheniusReactionRate);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeIRNReactions, ArrheniusReactionRate);

    //forCoeffGases(makeIRNReactions, LandauTellerReactionRate);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeIRNReactions, LandauTellerReactionRate);

    
    //forCoeffGases(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeIRNReactions, thirdBodyArrheniusReactionRate);

    // Irreversible/reversible reactions

    //forCoeffGases(makeIRReactions, JanevReactionRate);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeIRReactions, JanevReactionRate);

    //forCoeffGases(makeIRReactions, powerSeriesReactionRate);
    forCoeffGasEnergiesAndThermos(ddpolynomialTransport, makeIRReactions, powerSeriesReactionRate);


    // Pressure dependent reactions

    // forCoeffGases
    // (
    //     makeIRRPressureDependentReactions,
    //     FallOffReactionRate,
    //     ArrheniusReactionRate,
    //     LindemannFallOffFunction
    // );
    forCoeffGasEnergiesAndThermos
    (
        ddpolynomialTransport,
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );



    // forCoeffGases
    // (
    //     makeIRRPressureDependentReactions,
    //     FallOffReactionRate,
    //     ArrheniusReactionRate,
    //     TroeFallOffFunction
    // );
    forCoeffGasEnergiesAndThermos
    (
        ddpolynomialTransport,
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    // forCoeffGases
    // (
    //     makeIRRPressureDependentReactions,
    //     FallOffReactionRate,
    //     ArrheniusReactionRate,
    //     SRIFallOffFunction
    // );

    forCoeffGasEnergiesAndThermos
    (
        ddpolynomialTransport,
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );

    // forCoeffGases
    // (
    //     makeIRRPressureDependentReactions,
    //     ChemicallyActivatedReactionRate,
    //     ArrheniusReactionRate,
    //     LindemannFallOffFunction
    // );

    forCoeffGasEnergiesAndThermos
    (
        ddpolynomialTransport,
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    // forCoeffGases
    // (
    //     makeIRRPressureDependentReactions,
    //     ChemicallyActivatedReactionRate,
    //     ArrheniusReactionRate,
    //     TroeFallOffFunction
    // );
    forCoeffGasEnergiesAndThermos
    (
        ddpolynomialTransport,
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    // forCoeffGases
    // (
    //     makeIRRPressureDependentReactions,
    //     ChemicallyActivatedReactionRate,
    //     ArrheniusReactionRate,
    //     SRIFallOffFunction
    // );
    forCoeffGasEnergiesAndThermos
    (
        ddpolynomialTransport,
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );


}

// ************************************************************************* //

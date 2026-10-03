/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::turbulenceThermophysicalTransportModels::TabulatedFickian

Description
    Base class for multi-component TabulatedFickian based temperature gradient heat
    flux models with optional Soret thermal diffusion of species.

    The mixture diffusion coefficients are specified as a fitting polynomials relevant 
    to pressure and temperature but independent of composition.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

SourceFiles
    TabulatedFickian.C

\*---------------------------------------------------------------------------*/

#include "Function2.H"

#ifndef TabulatedFickian_H
#define TabulatedFickian_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class TabulatedFickian Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class TabulatedFickian
:
    public BasicThermophysicalTransportModel
{
    // Private data

        // Model coefficients

        //- List of specie Soret thermal diffusion coefficient
        //  functions [kg/m/s]
        PtrList<Function2<scalar>> DTFuncs_;

        //- Mass diffusion coefficients in the mixture
        PtrList<volScalarField> Dm_;

        //- Mass diffusion polynomials coefficients table
        List<List<FixedList<scalar, 5>>>  DpolyCoeffs_;


public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        TabulatedFickian
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~TabulatedFickian()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Update the diffusion coefficients
        virtual void correct();

        //- Return the binary diffusion coefficients given pressure, temperature and fitting coeffs
        inline scalar getDij(const FixedList<scalar, 5>&, const scalar&, const scalar&);

       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TabulatedFickian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

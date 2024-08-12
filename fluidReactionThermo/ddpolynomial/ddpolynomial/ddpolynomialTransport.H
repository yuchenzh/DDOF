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

Class
    Foam::ddpolynomialTransport

Description
    Transport package using ddpolynomial functions of temperature for the
    dynamic viscosity and thermal conductivity:

    \verbatim
        mu = muCoeffs[0] + muCoeffs[1]*T + muCoeffs[2]*sqr(T)
           + muCoeffs[3]*pow3(T) + muCoeffs[4]*pow4(T)
           + muCoeffs[5]*pow(T, 5) + muCoeffs[6]*pow(T, 6)
           + muCoeffs[7]*pow(T, 7)

        kappa = kappaCoeffs[0] + kappaCoeffs[1]*T + kappaCoeffs[2]*sqr(T)
           + kappaCoeffs[3]*pow3(T) + kappaCoeffs[4]*pow4(T)
           + kappaCoeffs[5]*pow(T, 5) + kappaCoeffs[6]*pow(T, 6)
           + kappaCoeffs[7]*pow(T, 7)
    \endverbatim

    The ddpolynomial function is templated on the order of the ddpolynomial which
    defaults to 8.

Usage
    \table
        Property        | Description
        muCoeffs<8>     | Dynamic viscosity ddpolynomial coefficients
        kappaCoeffs<8>  | Thermal conductivity ddpolynomial coefficients
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        muCoeffs<8>     (1000 -0.05 0.003 0 0 0 0 0);
        kappaCoeffs<8>  (2000 -0.15 0.023 0 0 0 0 0);
    }
    \endverbatim

SourceFiles
    ddpolynomialTransportI.H
    ddpolynomialTransport.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef ddpolynomialTransport_H
#define ddpolynomialTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo, int PolySize> class ddpolynomialTransport;

template<class Thermo, int PolySize>
inline ddpolynomialTransport<Thermo, PolySize> operator+
(
    const ddpolynomialTransport<Thermo, PolySize>&,
    const ddpolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline ddpolynomialTransport<Thermo, PolySize> operator*
(
    const scalar,
    const ddpolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const ddpolynomialTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class ddpolynomialTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class ddpolynomialTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity ddpolynomial coefficients [Pa.s/K^i]
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity ddpolynomial coefficients [W/m/K/K^i]
        Polynomial<PolySize> kappaCoeffs_;

        //- 1/kappa ddpolynomial coefficients
        Polynomial<PolySize> oneOverKappaCoeffs_;


    // Private Member Functions

        //- Construct from components
        inline ddpolynomialTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly,
            const Polynomial<PolySize>& oneOverKappaPoly
        );


public:

    // Constructors

        //- Construct as named copy
        inline ddpolynomialTransport(const word&, const ddpolynomialTransport&);

        //- Construct from dictionary
        ddpolynomialTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<ddpolynomialTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<ddpolynomialTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "ddpolynomial<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const ddpolynomialTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend ddpolynomialTransport operator+ <Thermo, PolySize>
        (
            const ddpolynomialTransport&,
            const ddpolynomialTransport&
        );

        friend ddpolynomialTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const ddpolynomialTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const ddpolynomialTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ddpolynomialTransportI.H"

#ifdef NoRepository
    #include "ddpolynomialTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

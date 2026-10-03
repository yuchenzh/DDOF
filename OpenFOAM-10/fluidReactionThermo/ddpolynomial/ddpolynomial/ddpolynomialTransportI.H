/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::ddpolynomialTransport<Thermo, PolySize>::ddpolynomialTransport
(
    const Thermo& t,
    const Polynomial<PolySize>& muCoeffs,
    const Polynomial<PolySize>& kappaCoeffs,
    const Polynomial<PolySize>& oneOverKappaCoeffs
)
:
    Thermo(t),
    muCoeffs_(muCoeffs),
    kappaCoeffs_(kappaCoeffs),
    oneOverKappaCoeffs_(oneOverKappaCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::ddpolynomialTransport<Thermo, PolySize>::ddpolynomialTransport
(
    const word& name,
    const ddpolynomialTransport& pt
)
:
    Thermo(name, pt),
    muCoeffs_(pt.muCoeffs_),
    kappaCoeffs_(pt.kappaCoeffs_),
    oneOverKappaCoeffs_(pt.oneOverKappaCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::ddpolynomialTransport<Thermo, PolySize>>
Foam::ddpolynomialTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<ddpolynomialTransport<Thermo, PolySize>>
    (
        new ddpolynomialTransport<Thermo, PolySize>(*this)
    );
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::ddpolynomialTransport<Thermo, PolySize>>
Foam::ddpolynomialTransport<Thermo, PolySize>::New(const dictionary& dict)
{
    return autoPtr<ddpolynomialTransport<Thermo, PolySize>>
    (
        new ddpolynomialTransport<Thermo, PolySize>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::ddpolynomialTransport<Thermo, PolySize>::mu
(
    const scalar p,
    const scalar T
) const
{
    return muCoeffs_.value(T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::ddpolynomialTransport<Thermo, PolySize>::kappa
(
    const scalar p,
    const scalar T
) const
{

    return 0.5*(kappaCoeffs_.value(T) + 1/oneOverKappaCoeffs_.value(T));
    //return kappaCoeffs_.value(T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::ddpolynomialTransport<Thermo, PolySize>::operator+=
(
    const ddpolynomialTransport<Thermo, PolySize>& pt
)
{
    scalar Y1 = this->Y();
    scalar W1 = this->W();

    Thermo::operator+=(pt);


    if (mag(this->Y()) > small)
    {
        scalar X1 = Y1/W1*this->W();
        scalar X2 = pt.Y()/pt.W()*this->W();
        X1 = X1/(X1 + X2);
        X2 = 1 - X1;
        Y1 /= this->Y();
        scalar Y2 = pt.Y()/this->Y();

 

        muCoeffs_ = Y1*muCoeffs_ + Y2*pt.muCoeffs_;
        kappaCoeffs_ = X1*kappaCoeffs_ + X2*pt.kappaCoeffs_;
        oneOverKappaCoeffs_ = X1*oneOverKappaCoeffs_ + X2*pt.oneOverKappaCoeffs_;
    }
}


template<class Thermo, int PolySize>
inline void Foam::ddpolynomialTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::ddpolynomialTransport<Thermo, PolySize> Foam::operator+
(
    const ddpolynomialTransport<Thermo, PolySize>& pt1,
    const ddpolynomialTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    if (mag(t.Y()) < small)
    {
        return ddpolynomialTransport<Thermo>
        (
            t,
            0,
            pt1.muCoeffs_,
            pt1.kappaCoeffs_
        );
    }
    else
    {
        scalar Y1 = pt1.Y()/t.Y();
        scalar Y2 = pt2.Y()/t.Y();

        // Calculate molar fraction
        scalar X1 = pt1.Y()/pt1.W();
        scalar X2 = pt2.Y()/pt2.W();
        X1 = X1/(X1 + X2);
        X2 = 1 - X1;

        return ddpolynomialTransport<Thermo, PolySize>
        (
            t,
            Y1*pt1.muCoeffs_ + Y2*pt2.muCoeffs_,
            X1*pt1.kappaCoeffs_ + X2*pt2.kappaCoeffs_,
            X1*pt1.oneOverKappaCoeffs_ + X2*pt2.oneOverKappaCoeffs_
        );
    }
}


template<class Thermo, int PolySize>
inline Foam::ddpolynomialTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const ddpolynomialTransport<Thermo, PolySize>& pt
)
{
    return ddpolynomialTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.muCoeffs_,
        pt.kappaCoeffs_,
        pt.oneOverKappaCoeffs_
    );
}


// ************************************************************************* //

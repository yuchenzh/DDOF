/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ddpolynomialTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::ddpolynomialTransport<Thermo, PolySize>::ddpolynomialTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    muCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "muCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    kappaCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "kappaCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    oneOverKappaCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "oneOverKappaCoeffs<" + Foam::name(PolySize) + '>'
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::ddpolynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add
    (
        word("muCoeffs<" + Foam::name(PolySize) + '>'),
        muCoeffs_
    );
    dict.add
    (
        word("kappaCoeffs<" + Foam::name(PolySize) + '>'),
        kappaCoeffs_
    );
    dict.add
    (
        word("oneOverKappaCoeffs<" + Foam::name(PolySize) + '>'),
        oneOverKappaCoeffs_
    );
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ddpolynomialTransport<Thermo, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //

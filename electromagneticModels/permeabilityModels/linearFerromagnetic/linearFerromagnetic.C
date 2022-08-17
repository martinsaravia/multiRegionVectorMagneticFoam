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

#include "linearFerromagnetic.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace permeabilityModels
{
    defineTypeNameAndDebug(linearFerromagnetic, 0);
    addToRunTimeSelectionTable
	(
		permeabilityModel,
		linearFerromagnetic,
		dictionary
	);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::permeabilityModels::linearFerromagnetic::linearFerromagnetic
(
    const word& name,
    const dictionary& permeabilityProperties,
    const volVectorField& B
)
:
    permeabilityModel(name, permeabilityProperties, B),  //Call Base class constructo to initialize its members
    mur0_("mur", dimless, permeabilityProperties_),
    mur_
    (
        IOobject
        (
            name,
            B_.time().timeName(),
            B_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        B_.mesh(),
        mur0_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::permeabilityModels::linearFerromagnetic::read
(
    const dictionary& permeabilityProperties
)
{
    permeabilityModel::read(permeabilityProperties);

    permeabilityProperties_.lookup("mur") >> mur0_;
    mur_ = mur0_;

    return true;
}


// ************************************************************************* //

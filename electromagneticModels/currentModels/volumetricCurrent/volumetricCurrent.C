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

#include "volumetricCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace currentModels
{
    defineTypeNameAndDebug(volumetricCurrent, 0);
    addToRunTimeSelectionTable
	(
		currentModel,
		volumetricCurrent,
		dictionary
	);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::currentModels::volumetricCurrent::volumetricCurrent
(
    const word& name,
    const dictionary& currentProperties,
    const volVectorField& B
)
:
    currentModel(name, currentProperties, B),  //Call Base class constructo to initialize its members
    J0_("J", dimensionSet(0,-2,0,0,0,1,0), currentProperties_),
    J_
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
        J0_
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::currentModels::volumetricCurrent::read
(
    const dictionary& currentProperties
)
{
    currentModel::read(currentProperties);

	  currentProperties_.lookup("J") >> J0_;

    J_ = J0_;

    return true;
}


// ************************************************************************* //

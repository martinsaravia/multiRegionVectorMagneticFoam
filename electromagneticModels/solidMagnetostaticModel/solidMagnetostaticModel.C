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

#include "solidMagnetostaticModel.H"
#include "addToRunTimeSelectionTable.H"

#include "permeabilityModel.H"
#include "magnetizationModel.H"
#include "currentModel.H"

#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidMagnetostaticModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidMagnetostaticModel::solidMagnetostaticModel
(
    const volVectorField& B
)
:
    IOdictionary
    (
        IOobject
        (
            "magneticProperties",
            B.time().constant(),
            B.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    permeabilityModelPtr_
  	(
  		permeabilityModel::New
  		(
  			"mur",
  			subDict("permeability"),
  			B
  		)
  	),

    magnetizationModelPtr_
  	(
  		magnetizationModel::New
  		(
  			"M",
  			subDict("magnetization"),
  			B
  		)
  	),

    currentModelPtr_
    (
      currentModel::New
      (
        "J",
        subDict("current"),
        B
      )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidMagnetostaticModel::~solidMagnetostaticModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::solidMagnetostaticModel::mur() const
{
    return permeabilityModelPtr_->mur();
}

Foam::tmp<Foam::volVectorField>
Foam::solidMagnetostaticModel::M() const
{
    return magnetizationModelPtr_->M();
}

Foam::tmp<Foam::volVectorField>
Foam::solidMagnetostaticModel::J() const
{
    Info<<"    point1 !\n"  << endl;

    return currentModelPtr_->J();
}

void Foam::solidMagnetostaticModel::correct()
{
  permeabilityModelPtr_->correct();
	magnetizationModelPtr_->correct();
  currentModelPtr_->correct();
}


bool Foam::solidMagnetostaticModel::read()
{
   return regIOobject::read() &&
           permeabilityModelPtr_->read(subDict("permeability")) &&
           magnetizationModelPtr_->read(subDict("magnetization")) &&
           currentModelPtr_->read(subDict("current"));
}


// ************************************************************************* //

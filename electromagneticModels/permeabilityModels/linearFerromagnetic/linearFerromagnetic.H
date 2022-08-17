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

Class
    Foam::permeabilityModels::linearFerromagnetic

Description
    An linear ferromagnetic magnetic model.

SourceFiles
    linearFerromagnetic.C

\*---------------------------------------------------------------------------*/

#ifndef linearFerromagnetic_H
#define linearFerromagnetic_H

#include "permeabilityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace permeabilityModels
{

/*---------------------------------------------------------------------------*\
                           Class linearFerromagnetic Declaration
\*---------------------------------------------------------------------------*/

class linearFerromagnetic
:
    public permeabilityModel
{
    // Private data

        dimensionedScalar mur0_;

        volScalarField mur_;


public:

    //- Runtime type information
    TypeName("linearFerromagnetic");


    // Constructors

        //- Construct from components
        linearFerromagnetic
        (
            const word& name,
            const dictionary& permeabilityProperties,
            const volVectorField& B
        );


    //- Destructor
    virtual ~linearFerromagnetic()
    {}


    // Member Functions

        //- Return the relative permeability
        virtual tmp<volScalarField> mur() const
        {
            return mur_;
        }		

        //- Correct the relative permeability (not appropriate, viscosity constant)
        virtual void correct()
        {}

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& permeabilityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace permeabilityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

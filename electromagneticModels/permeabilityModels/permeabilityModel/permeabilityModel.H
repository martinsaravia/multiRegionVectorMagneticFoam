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

Namespace
    Foam::permeabilityModels

Description
    A namespace for various  permeabilityModel implementations.

Class
    Foam::permeabilityModel

Description
    An abstract base class for permeabilityModels.

SourceFiles
    permeabilityModel.C
    permeabilityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef permeabilityModel_H
#define permeabilityModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class permeabilityModel Declaration
\*---------------------------------------------------------------------------*/

class permeabilityModel
{

protected:

    // Protected data

        word name_;
        dictionary permeabilityProperties_;

        const volVectorField& B_;

    // Private Member Functions

        //- Disallow copy construct
        permeabilityModel(const permeabilityModel&);

        //- Disallow default bitwise assignment
        void operator=(const permeabilityModel&);


public:

    //- Runtime type information
    TypeName("permeabilityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            permeabilityModel,
            dictionary,
            (
                const word& name,
                const dictionary& permeabilityProperties,
                const volVectorField& B
            ),
            (name, permeabilityProperties, B)
        );


    // Selectors

        //- Return a reference to the selected permeability model
        static autoPtr<permeabilityModel> New
        (
            const word& name,
            const dictionary& permeabilityProperties,
            const volVectorField& B
        );


    // Constructors

        //- Construct from components
        permeabilityModel
        (
            const word& name,
            const dictionary& permeabilityProperties,
            const volVectorField& B
        );


    //- Destructor
    virtual ~permeabilityModel()
    {}


    // Member Functions

        //- Return the permeability properties dictionary
        const dictionary& properties() const
        {
            return permeabilityProperties_;
        }

        //- Return the relative permeability
        virtual tmp<volScalarField> mur() const = 0;	

        //- Correct the relative permeability
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& permeabilityProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

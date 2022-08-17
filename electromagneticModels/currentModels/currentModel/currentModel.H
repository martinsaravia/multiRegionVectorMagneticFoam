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
    Foam::currentModels

Description
    A namespace for various  currentModel implementations.

Class
    Foam::currentModel

Description
    An abstract base class for currentModels.

SourceFiles
    currentModel.C
    currentModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef currentModel_H
#define currentModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class currentModel Declaration
\*---------------------------------------------------------------------------*/

class currentModel
{

protected:

    // Protected data

        word name_;
        dictionary currentProperties_;

        const volVectorField& B_;

    // Private Member Functions

        //- Disallow copy construct
        currentModel(const currentModel&);

        //- Disallow default bitwise assignment
        void operator=(const currentModel&);


public:

    //- Runtime type information
    TypeName("currentModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            currentModel,
            dictionary,
            (
                const word& name,
                const dictionary& currentProperties,
                const volVectorField& B
            ),
            (name, currentProperties, B)
        );


    // Selectors

        //- Return a reference to the selected current model
        static autoPtr<currentModel> New
        (
            const word& name,
            const dictionary& currentProperties,
            const volVectorField& B
        );


    // Constructors

        //- Construct from components
        currentModel
        (
            const word& name,
            const dictionary& currentProperties,
            const volVectorField& B
        );


    //- Destructor
    virtual ~currentModel()
    {}


    // Member Functions

        //- Return the current properties dictionary
        const dictionary& currentProperties() const
        {
            return currentProperties_;
        }

        //- Return the relative current
        virtual tmp<volVectorField> J() const = 0;

        //- Correct the relative current
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& currentProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

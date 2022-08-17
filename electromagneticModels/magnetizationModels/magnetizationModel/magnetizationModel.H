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
    Foam::magnetizationModels

Description
    A namespace for various  magnetizationModel implementations.

Class
    Foam::magnetizationModel

Description
    An abstract base class for magnetizationModels.

SourceFiles
    magnetizationModel.C
    magnetizationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef magnetizationModel_H
#define magnetizationModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class magnetizationModel Declaration
\*---------------------------------------------------------------------------*/

class magnetizationModel
{

protected:

    // Protected data

        word name_;
        dictionary magnetizationProperties_;

        const volVectorField& B_;

    // Private Member Functions

        //- Disallow copy construct
        magnetizationModel(const magnetizationModel&);

        //- Disallow default bitwise assignment
        void operator=(const magnetizationModel&);


public:

    //- Runtime type information
    TypeName("magnetizationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            magnetizationModel,
            dictionary,
            (
                const word& name,
                const dictionary& magnetizationProperties,
                const volVectorField& B
            ),
            (name, magnetizationProperties, B)
        );


    // Selectors

        //- Return a reference to the selected magnetization model
        static autoPtr<magnetizationModel> New
        (
            const word& name,
            const dictionary& magnetizationProperties,
            const volVectorField& B
        );


    // Constructors

        //- Construct from components
        magnetizationModel
        (
            const word& name,
            const dictionary& magnetizationProperties,
            const volVectorField& B
        );


    //- Destructor
    virtual ~magnetizationModel()
    {}


    // Member Functions

        //- Return the magnetization properties dictionary
        const dictionary& magnetizationProperties() const
        {
            return magnetizationProperties_;
        }

        //- Return the relative magnetization
        virtual tmp<volVectorField> M() const = 0;	

        //- Correct the relative magnetization
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& magnetizationProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

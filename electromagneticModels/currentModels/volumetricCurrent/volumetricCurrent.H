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
    Foam::currentModels::volumetricCurrent

Description
    A volumetric current model.

SourceFiles
    volumetricCurrent.C

\*---------------------------------------------------------------------------*/

#ifndef volumetricCurrent_H
#define volumetricCurrent_H

#include "currentModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace currentModels
{

/*---------------------------------------------------------------------------*\
                           Class volumetricCurrent Declaration
\*---------------------------------------------------------------------------*/

class volumetricCurrent
:
    public currentModel
{
    // Private data

		dimensionedVector J0_;

		volVectorField J_;


public:

    //- Runtime type information
    TypeName("volumetricCurrent");


    // Constructors

        //- Construct from components
        volumetricCurrent
        (
            const word& name,
            const dictionary& currentProperties,
            const volVectorField& B
        );


    //- Destructor
    virtual ~volumetricCurrent()
    {}


    // Member Functions

        //- Return the current
        virtual tmp<volVectorField> J() const
        {
            return J_;
        }

        //- Correct the current (not appropriate, current constant)
        virtual void correct()
        {}

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& currentProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace currentModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

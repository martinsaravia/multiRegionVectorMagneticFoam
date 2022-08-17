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
    Foam::magnetizationModels::hardMagnetization

Description
    A hard magnetization model.

SourceFiles
    hardMagnetization.C

\*---------------------------------------------------------------------------*/

#ifndef hardMagnetization_H
#define hardMagnetization_H

#include "magnetizationModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace magnetizationModels
{

/*---------------------------------------------------------------------------*\
                           Class hardMagnetization Declaration
\*---------------------------------------------------------------------------*/

class hardMagnetization
:
    public magnetizationModel
{
    // Private data

		dimensionedVector M0_;

		volVectorField M_;


public:

    //- Runtime type information
    TypeName("hardMagnetization");


    // Constructors

        //- Construct from components
        hardMagnetization
        (
            const word& name,
            const dictionary& magnetizationProperties,
            const volVectorField& B
        );


    //- Destructor
    virtual ~hardMagnetization()
    {}


    // Member Functions

        //- Return the magnetization
        virtual tmp<volVectorField> M() const
        {
            return M_;
        }

        //- Correct the magnetization (not appropriate, magnetization constant)
        virtual void correct()
        {}

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& magnetizationProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace magnetizationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::electromagneticModel

Description
    Base-class for all electromagnetic models used by the incompressible turbulence
    models.

SourceFiles
    electromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef electromagneticModel_H
#define electromagneticModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class electromagneticModel Declaration
\*---------------------------------------------------------------------------*/

class electromagneticModel
{
    // Private Member Functions

        //- Disallow copy construct
        electromagneticModel(const electromagneticModel&);

        //- Disallow default bitwise assignment
        void operator=(const electromagneticModel&);

public:

    //- Runtime type information
    TypeName("electromagneticModel");


    // Constructors

        //- Construct from components
        electromagneticModel
        ();


    //- Destructor
    virtual ~electromagneticModel();


    // Member Functions

        //- Return the relative permeability
        virtual tmp<volScalarField> mur() const = 0;

        //- Return the magnetization vector
        virtual tmp<volVectorField> M() const = 0;

        //- Return the magnetization vector
        virtual tmp<volVectorField> J() const = 0;

        //- Correct the fields
        virtual void correct() = 0;

        //- Read electromagneticProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

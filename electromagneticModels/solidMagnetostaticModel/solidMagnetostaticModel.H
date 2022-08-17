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
    Foam::solidMagnetostaticModel

Description
    A model for multi-region magnetostatics.

SourceFiles
    solidMagnetostaticModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidMagnetostaticModel_H
#define solidMagnetostaticModel_H

#include "electromagneticModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class permeabilityModel;
class magnetizationModel;  //Forward declaration
class currentModel;

/*---------------------------------------------------------------------------*\
                Class solidMagnetostaticModel Declaration
\*---------------------------------------------------------------------------*/

class solidMagnetostaticModel
:
    public IOdictionary,
    public electromagneticModel
{
    // Private Data

        autoPtr<permeabilityModel> permeabilityModelPtr_;
		    autoPtr<magnetizationModel> magnetizationModelPtr_;
        autoPtr<currentModel> currentModelPtr_;


    // Private Member Functions

        //- Disallow copy construct
        solidMagnetostaticModel(const solidMagnetostaticModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidMagnetostaticModel&);


public:

    //- Runtime type information
    TypeName("solidMagnetostaticModel");


    // Constructors

        //- Construct from components
        solidMagnetostaticModel
        (
            const volVectorField& B
        );


    //- Destructor
    virtual ~solidMagnetostaticModel();


    // Member Functions

        //- Return the relative permeability
        virtual tmp<volScalarField> mur() const;

		    //- Return the magnetization
        virtual tmp<volVectorField> M() const;

        //- Return the current
        virtual tmp<volVectorField> J() const;

        //- Correct the fields
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

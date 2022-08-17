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
    Foam::CVIFvPatchVectorField

Description
    This boundary condition supplies a magnetic interface constraint.

Usage
    \table
        Property     | Description             | Required    | Default value
        value        | Patch face values       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            magneticInterface;
        value           uniform 0;  // Example for scalar field usage
    }
    \endverbatim

SourceFiles
   CVIFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef CVIFvPatchVectorField_H
#define CVIFvPatchVectorField_H

/* #include "fvPatchField.H" */
#include "magneticInterfaceFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CVIFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class CVIFvPatchVectorField
:
    public magneticInterfaceFvPatchVectorField
{
  // Private data

      //- Name of the neighbour patch
      word nbrPatchName_;

      //- Name of the neighbour mesh
      word nbrMeshName_;

      //- Gradient internal coefficients
      Field<vector> gradICoeffs_;

      //- Gradient boundary coefficients
      Field<vector> gradBCoeffs_;

public:

    //- Runtime type information
    TypeName("CVI");


    // Constructors

        //- Construct from patch and internal field
        CVIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        /* //- Construct from patch, internal field and value
        CVIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const vector& value
        ); */

        //- Construct from patch, internal field and dictionary
        CVIFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given CVIFvPatchVectorField
        //  onto a new patch
        CVIFvPatchVectorField
        (
            const CVIFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        CVIFvPatchVectorField(const CVIFvPatchVectorField&) = delete;

        //- Construct as copy setting internal field reference
        CVIFvPatchVectorField
        (
            const CVIFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new CVIFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            virtual tmp<Field<vector>> gradICoeffs() const
            {
                return gradICoeffs_;
            }

            virtual tmp<Field<vector>> gradBCoeffs() const
            {
                return gradBCoeffs_;
            }


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            /* using magneticInterfaceFvPatchVectorField::updateCoeffs; */

            // Update Coefficients
            virtual void updateCoeffs();



        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<vector>&) {}

        virtual void operator=(const fvPatchField<vector>&) {}
        virtual void operator+=(const fvPatchField<vector>&) {}
        virtual void operator-=(const fvPatchField<vector>&) {}
        virtual void operator*=(const fvPatchScalarField&) {}
        virtual void operator/=(const fvPatchScalarField&) {}

        virtual void operator+=(const Field<vector>&) {}
        virtual void operator-=(const Field<vector>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const vector&) {}
        virtual void operator+=(const vector&) {}
        virtual void operator-=(const vector&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

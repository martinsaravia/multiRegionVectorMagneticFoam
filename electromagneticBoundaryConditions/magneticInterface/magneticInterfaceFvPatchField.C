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

#include "magneticInterfaceFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
Foam::magneticInterfaceFvPatchField<Type>::magneticInterfaceFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    nbrPatchName_("none"),
    nbrMeshName_("none")
{}

template<class Type>
Foam::magneticInterfaceFvPatchField<Type>::magneticInterfaceFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    fvPatchField<Type>(p, iF, value),
    nbrPatchName_("none"),
    nbrMeshName_("none")
{}

template<class Type>
Foam::magneticInterfaceFvPatchField<Type>::magneticInterfaceFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<Type>(p, iF, dict, valueRequired),
    nbrPatchName_("none"),
    nbrMeshName_("none")
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
      FatalErrorIn
      (
          "magneticInterfaceFvPatchField::"
          "magneticInterfaceFvPatchField\n"
          "(\n"
          "    const fvPatch& p,\n"
          "    const DimensionedField<Type, volMesh>& iF,\n"
          "    const dictionary& dict\n"
          ")\n"
      )   << "\n    patch type '" << p.type()
          << "' not type '" << mappedPatchBase::typeName << "'"
          << "\n    for patch " << p.name()
          << exit(FatalError);
    }

  fvPatchField<Type>::operator=(Field<Type>("value", dict, p.size()));
}

template<class Type>
Foam::magneticInterfaceFvPatchField<Type>::magneticInterfaceFvPatchField
(
    const magneticInterfaceFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    nbrPatchName_(ptf.nbrPatchName_),
    nbrMeshName_(ptf.nbrMeshName_)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::magneticInterfaceFvPatchField<Type>::magneticInterfaceFvPatchField
(
    const magneticInterfaceFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    nbrPatchName_(ptf.nbrPatchName_),
    nbrMeshName_(ptf.nbrMeshName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Template specialiation for the vector variant of valueInternalCoeffs
// inline must be added in order to avoid the compiler to complain againts
//
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::snGrad() const
{
    Info << "WARNING: CALLING not implemented snGrad in magneticInterfaceFvPatchField... " << endl;
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::valICoeffs() const
{
    Info << "WARNING: CALLING not implemented valICoeffs in magneticInterfaceFvPatchField... " << endl;
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::valBCoeffs() const
{
    Info << "WARNING: CALLING not implemented valBCoeffs() in magneticInterfaceFvPatchField"<< endl;
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::gradICoeffs() const
{
    Info << "WARNING: CALLING not implemented gradICoeffs() in magneticInterfaceFvPatchField"<< endl;
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::gradBCoeffs() const
{
    Info << "WARNING: CALLING not implemented gradBCoeffs() in magneticInterfaceFvPatchField"<< endl;
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


template <class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    Info<<"hete"<<endl;
    return valICoeffs();
}


template <class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    Info<<"here"<<endl;
    return valBCoeffs();
}


template <class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::gradientInternalCoeffs() const
{
    return gradICoeffs();
}


template <class Type>
Foam::tmp<Foam::Field<Type>>
Foam::magneticInterfaceFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return gradBCoeffs();
}


template<class Type>
void Foam::magneticInterfaceFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    //this->writeEntry("value", os);  // line for OpenFOAM 6
    writeEntry(os, "value", *this); // line for OpenFOAM 7
}


// ************************************************************************* //

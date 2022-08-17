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

#include "MVIFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "electromagneticConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::MVIFvPatchVectorField::MVIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    magneticInterfaceFvPatchVectorField(p, iF),
    nbrPatchName_("none"),
    nbrMeshName_("none"),
    gradICoeffs_(0),
    gradBCoeffs_(0),
    valICoeffs_(0),
    valBCoeffs_(0),
    snGradCoeffs_(this->patch().nf())
{}


/* Foam::MVIFvPatchVectorField::MVIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const vector& value
)
:
    magneticInterfaceFvPatchVectorField(p, iF, value),
    nbrPatchName_("none"),
    nbrMeshName_("none"),
    gradICoeffs_(0),
    gradBCoeffs_(0)
{} */


Foam::MVIFvPatchVectorField::MVIFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    magneticInterfaceFvPatchVectorField(p, iF, dict),
    nbrPatchName_("none"),
    nbrMeshName_("none"),
    gradICoeffs_(0),
    gradBCoeffs_(0),
    valICoeffs_(0),
    valBCoeffs_(0),
    snGradCoeffs_(this->patch().nf())  // Just Fill with something to start (not the correct value)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
      FatalErrorIn
      (
          "magneticInterfaceFvPatchField::"
          "magneticInterfaceFvPatchField\n"
          "(\n"
          "    const fvPatch& p,\n"
          "    const DimensionedField<vector, volMesh>& iF,\n"
          "    const dictionary& dict\n"
          ")\n"
      )   << "\n    patch type '" << p.type()
          << "' not type '" << mappedPatchBase::typeName << "'"
          << "\n    for patch " << p.name()
          << exit(FatalError);
    }
  magneticInterfaceFvPatchVectorField::operator=(Field<vector>("value", dict, p.size()));
}


Foam::MVIFvPatchVectorField::MVIFvPatchVectorField
(
    const MVIFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    magneticInterfaceFvPatchVectorField(ptf, p, iF, mapper),
    nbrPatchName_(ptf.nbrPatchName_),
    nbrMeshName_(ptf.nbrMeshName_),
    gradICoeffs_(ptf.gradICoeffs_),
    gradBCoeffs_(ptf.gradBCoeffs_),
    valICoeffs_(ptf.valICoeffs_),
    valBCoeffs_(ptf.valBCoeffs_),
    snGradCoeffs_(ptf.snGradCoeffs_)
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


Foam::MVIFvPatchVectorField::MVIFvPatchVectorField
(
    const MVIFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    magneticInterfaceFvPatchVectorField(ptf, iF),
    gradICoeffs_(ptf.gradICoeffs_),
    gradBCoeffs_(ptf.gradBCoeffs_),
    valICoeffs_(ptf.valICoeffs_),
    valBCoeffs_(ptf.valBCoeffs_),
    snGradCoeffs_(ptf.snGradCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MVIFvPatchVectorField::updateCoeffs()
{
if (this->updated())
{
    return;
}

// Since we're inside initEvaluate/evaluate there might be processor
// comms underway. Change the tag we use.
int oldTag = UPstream::msgType();
UPstream::msgType() = oldTag+1;

// Get the coupling information from the mappedPatchBase
// Basically retrieve nbrMesh and nbrPatch
const mappedPatchBase& mpp = refCast<const mappedPatchBase>
(
    MVIFvPatchVectorField::patch().patch()
);

// Get this and neigubour meshes
//const fvMesh& thisMesh = this->patch().boundaryMesh().mesh();
const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());

// Get the neighbour patch
const fvPatch& nbrPatch = nbrMesh.boundary()[mpp.samplePolyPatch().index()];

// Force recalculation of mapping and schedule
const distributionMap& distMap = mpp.map();

// Get the index of the neighbour patch
//const label samplePatchi = mpp.samplePolyPatch().index();
//const label nbrPatchID =
//  nbrMesh.boundaryMesh().findPatchID(mpp.samplePatch());

// Get the neighbour Patch fields as MVI patch fields
const MVIFvPatchVectorField& nbrPatchFieldA =
    refCast<const MVIFvPatchVectorField>
    (
      nbrPatch.lookupPatchField<volVectorField, vector>("A")
    );
const fvPatchField<vector>& nbrPatchFieldM =
      nbrPatch.lookupPatchField<volVectorField, vector>("M");
const fvPatchField<vector>& nbrPatchFieldB =
      nbrPatch.lookupPatchField<volVectorField, vector>("B");
const fvPatchField<scalar>& nbrPatchFieldChi =
      nbrPatch.lookupPatchField<volScalarField, scalar>("chi");

// Swap to obtain full local values of neighbour internal field
tmp<Field<vector>> nbrPatchIntFieldA(new Field<vector>(nbrPatchFieldA.size(), Zero));
tmp<Field<vector>> nbrPatchIntFieldM(new Field<vector>(nbrPatchFieldA.size(), Zero));
tmp<Field<vector>> nbrPatchIntFieldB(new Field<vector>(nbrPatchFieldA.size(), Zero));
tmp<Field<scalar>> nbrPatchIntFieldChi(new Field<scalar>(nbrPatchFieldA.size(), 0.0));

// Fill local neighbour data with internal fields
nbrPatchIntFieldA.ref() = nbrPatchFieldA.patchInternalField();
nbrPatchIntFieldM.ref() = nbrPatchFieldM.patchInternalField();
nbrPatchIntFieldB.ref() = nbrPatchFieldB.patchInternalField();
nbrPatchIntFieldChi.ref() = nbrPatchFieldChi.patchInternalField();

// Distribute
mpp.distribute(nbrPatchIntFieldA.ref());
mpp.distribute(nbrPatchIntFieldM.ref());
mpp.distribute(nbrPatchIntFieldB.ref());
mpp.distribute(nbrPatchIntFieldChi.ref());

// Get the A and M internal fields associated to the this patch
Field<vector> basePatchIntFieldA = this->patchInternalField();
Field<vector> basePatchIntFieldM =
  this->patch().template lookupPatchField<volVectorField, vector>("M").patchInternalField();
Field<vector> basePatchIntFieldB =
  this->patch().template lookupPatchField<volVectorField, vector>("B").patchInternalField();
Field<scalar> basePatchIntFieldChi =
  this->patch().template lookupPatchField<volScalarField, scalar>("chi").patchInternalField();

// Return base patch face normals (unit vectors)
Field<vector> n = this->patch().nf();

// Surface current vector
Field<vector> deltaM0 = nbrPatchIntFieldM - basePatchIntFieldM;
Field<vector> deltaChiB = nbrPatchIntFieldChi * nbrPatchIntFieldB
                        - basePatchIntFieldChi * basePatchIntFieldB;

const dimensionedScalar& mu0 = constant::electromagnetic::mu0;
Field<vector> K = mu0.value() * (deltaM0 + deltaChiB) ^ n ;


// Cell center to face center distances
Field<scalar> deltaEe(mag(nbrPatch.delta()));
distMap.distribute(deltaEe);
Field<scalar> deltaCe(mag(this->patch().delta()));
Field<scalar> deltaEC(deltaEe + deltaCe);

// Coefficients (for the laplacian)
gradICoeffs_ = -pTraits<vector>::one * (scalar(1.0)/deltaEC);
gradBCoeffs_ = (nbrPatchIntFieldA - deltaEe * K)/deltaEC;

// Average the field to half (correct, it is not ok, use weights)
valICoeffs_ = pTraits<vector>::one *0.5;
valBCoeffs_ =  nbrPatchFieldA * 0.5;


// Assign values to the fvPathField
Field<vector> nGradA = (scalar(-1.0)/deltaEC)*basePatchIntFieldA + gradBCoeffs_;
snGradCoeffs_ = nGradA;  // Normal gradient

// Assign the values to the field (values to give good B)
this->operator==(basePatchIntFieldA + nGradA * deltaCe);

// Values for continuity of A (bad B )
//this->operator==(0.5*basePatchIntFieldA + valBCoeffs_ );

magneticInterfaceFvPatchVectorField::updateCoeffs();

UPstream::msgType() = oldTag;
}


void Foam::MVIFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    //this->writeEntry("value", os);  // line for OpenFOAM 6
    writeEntry(os, "value", *this); // line for OpenFOAM 7
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        MVIFvPatchVectorField
    );
}

// ************************************************************************* //

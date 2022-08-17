Info<< "Creating fields..." << endl;

// List of pointers for the fields of each region
PtrList<solidMagnetostaticModel> solidModels(allRegions.size());

PtrList<volVectorField> As(allRegions.size());

PtrList<volScalarField> murs(allRegions.size());
PtrList<volScalarField> mus(allRegions.size());
PtrList<volScalarField> chis(allRegions.size());
PtrList<volVectorField> Js(allRegions.size());

PtrList<volVectorField> Bs(allRegions.size());
PtrList<volVectorField> Brs(allRegions.size());
PtrList<volVectorField> Hs(allRegions.size());
PtrList<volVectorField> Hrs(allRegions.size());
PtrList<volVectorField> Ms(allRegions.size());
PtrList<volVectorField> MIs(allRegions.size());
PtrList<volVectorField> MIrs(allRegions.size());

PtrList<volVectorField> MTs(allRegions.size());
PtrList<volVectorField> MTrs(allRegions.size());

PtrList<volVectorField> Fls(allRegions.size());
PtrList<volVectorField> Fks(allRegions.size());
PtrList<volVectorField> Fkrs(allRegions.size());
PtrList<volVectorField> Fk2s(allRegions.size());
PtrList<volVectorField> Fk2rs(allRegions.size());
// PtrList<surfaceVectorField> Fbs(allRegions.size());
// PtrList<surfaceVectorField> Fbrs(allRegions.size());
// PtrList<surfaceVectorField> Fb2s(allRegions.size());
// PtrList<surfaceVectorField> Fb2rs(allRegions.size());

PtrList<volTensorField> Ts(allRegions.size());
PtrList<volTensorField> Trs(allRegions.size());

const dimensionedScalar& mu0 = constant::electromagnetic::mu0;

// Identity tensor
const dimensionedTensor T1
        (
            "T1",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            tensor(1,0,0, 0,1,0, 0,0,1)
        );

// Populate field pointer lists
forAll(allRegions, i)
{
    Info<< "	Reading magnetic properties for region: "
        << allRegions[i].name() << endl;


    Info<< "		Adding A to region" << endl;
    As.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "A",
                runTime.timeName(),
                allRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            allRegions[i]
        )
    );


    Info<< "		Adding B to region" << endl;
    Bs.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "B",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::curl(As[i])
        )
    );

    Info<< "		Adding solidMagnetostaticModel to region" << endl;
    solidModels.set
    (
        i,
        new solidMagnetostaticModel(Bs[i])
    );

    Info<< "    	Adding mur to region" << endl;
    murs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mur",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
			solidModels[i].mur()
        )
    );

    Info<< "    	Adding mu to region" << endl;
    mus.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mu",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
			mu0 * murs[i]
        )
    );

    Info<< "		Initializing chi to zero for region" << endl;
	chis.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"chi",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			(murs[i] - 1) / (murs[i] * mu0)
		)
	);

    Info<< "		Creating J for region" << endl;
	Js.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"J",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            solidModels[i].J()
		)
	);


    Info<< "		Adding B to region" << endl;
    Brs.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Br",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::curl(As[i])
        )
    );

    Info<< "		Creating H for region" << endl;
	Hs.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"H",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            (1 / mus[i]) * Bs[i]
		)
	);

    Info<< "		Creating H for region" << endl;
    Hrs.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Hr",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (1 / mus[i]) * Brs[i]
        )
    );


	Info<< "		Adding M to region" << endl;
	Ms.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "M",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
			solidModels[i].M()
        )
    );

    Info<< "		Creating Mi for region" << endl;
	MIs.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"MI",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            chis[i] * Bs[i]
		)
	);

    Info<< "		Creating Mri for region" << endl;
    MIrs.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "MIr",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            chis[i] * Brs[i]
        )
    );

    Info<< "		Creating Mi for region" << endl;
	MTs.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"MT",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            Ms[i] + chis[i] * Bs[i]
		)
	);

    Info<< "		Creating Mri for region" << endl;
    MTrs.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "MTr",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ms[i] + chis[i] * Brs[i]
        )
    );

    Info<< "		Creating the Lorentz Force for region" << endl;
	Fls.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"Fl",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            mu0 * Js[i] ^ Hs[i]
		)
	);
    Info<< "		Creating the Kelvin Force for region" << endl;
	Fks.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"Fk",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            mu0 * MTs[i] & fvc::grad(Hs[i])
		)
	);

    Info<< "		Creating the Kelvin Force for region" << endl;
	Fkrs.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"Fkr",
				runTime.timeName(),
				allRegions[i],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
            mu0 * (MTrs[i]) & fvc::grad(Hrs[i])
		)
	);

    Info<< "		Creating the Kelvin Force for region" << endl;
    Fk2s.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Fk2",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (MTs[i] & fvc::grad(Bs[i])) + (MTs[i] ^ fvc::curl(Bs[i]))
        )
    );

    Info<< "		Creating the Kelvin Force for region" << endl;
    Fk2rs.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "Fk2r",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (MTrs[i] & fvc::grad(Brs[i])) + (MTrs[i] ^ fvc::curl(Brs[i]))
        )
    );

    // Info<< "		Creating the Force for the boundary" << endl;
	// Fbs.set
	// (
	// 	i,
	// 	new surfaceVectorField
	// 	(
	// 		IOobject
	// 		(
	// 			"Fb",
	// 			runTime.timeName(),
	// 			allRegions[i],
	// 			IOobject::NO_READ,
	// 			IOobject::AUTO_WRITE
	// 		),
    //         allRegions[i],
    //         dimensionedVector(dimForce/dimArea, Zero)
	// 	)
	// );
    //
    // Info<< "		Creating the Force for the boundary" << endl;
	// Fbrs.set
	// (
	// 	i,
	// 	new surfaceVectorField
	// 	(
	// 		IOobject
	// 		(
	// 			"Fbr",
	// 			runTime.timeName(),
	// 			allRegions[i],
	// 			IOobject::NO_READ,
	// 			IOobject::AUTO_WRITE
	// 		),
    //         allRegions[i],
    //         dimensionedVector(dimForce/dimArea, Zero)
	// 	)
	// );
    //
    //
    // Info<< "		Creating the Force for the boundary" << endl;
	// Fb2s.set
	// (
	// 	i,
	// 	new surfaceVectorField
	// 	(
	// 		IOobject
	// 		(
	// 			"Fb2",
	// 			runTime.timeName(),
	// 			allRegions[i],
	// 			IOobject::NO_READ,
	// 			IOobject::AUTO_WRITE
	// 		),
    //         allRegions[i],
    //         dimensionedVector(dimForce/dimArea, Zero)
	// 	)
	// );
    //
    // Info<< "		Creating the Force for the boundary" << endl;
	// Fb2rs.set
	// (
	// 	i,
	// 	new surfaceVectorField
	// 	(
	// 		IOobject
	// 		(
	// 			"Fb2r",
	// 			runTime.timeName(),
	// 			allRegions[i],
	// 			IOobject::NO_READ,
	// 			IOobject::AUTO_WRITE
	// 		),
    //         allRegions[i],
    //         dimensionedVector(dimForce/dimArea, Zero)
	// 	)
	// );

    Info<< "		Adding T to region" << endl;
    Ts.set
    (
        i,
        new volTensorField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (1.0 / mu0) * ( Bs[i] * Bs[i] - 0.5 * (Bs[i] & Bs[i]) * T1)
        )
    );

    Info<< "		Adding T to region" << endl;
    Trs.set
    (
        i,
        new volTensorField
        (
            IOobject
            (
                "Tr",
                runTime.timeName(),
                allRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (1.0 / mu0) * ( Brs[i] * Brs[i] - 0.5 * (Brs[i] & Brs[i]) * T1)
        )
    );

	Info<< "Fields created !\n" << endl;
}

// Create a function object for writing the residuals of A (incomplete)
// dictionary residualDict();
// residualDict.add("fields", "A");
// const functionObjects::residual& residual(residualDict);

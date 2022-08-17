// Mur scaling variant, not very stable


scalar ppf(0.0); // Parameter Proportional Factor
scalar iStop(runTime.endTime().value()-iCutOff); // Iteration of Stop
scalar iDelta(iStop-iCutIn);

// Scale chi field with the iteration number in order to help convergence
if ( runTime.value() > iCutIn && runTime.value() <= iStop)
{
    ppf = 1.0/(iDelta)*runTime.value() - iCutIn/(iDelta);
}
else if(runTime.value() > iStop)
{
    ppf = 1.0;
}
Info << "Mur Scaling Factor = " << ppf << endl;


PtrList<volScalarField> murs(allRegions.size());
PtrList<volScalarField> chis(allRegions.size());

Info<< "    	Adding mur to region" << endl;
forAll(allRegions, i)
{
Info<< "    	Mur for region "
    << allRegions[i].name()
    << " is: "
    << ppf * (solidModels[i].mur()-1) + 1 << endl;

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
                IOobject::	NO_READ,
                IOobject::AUTO_WRITE
            ),
            ppf * (solidModels[i].mur()-1) + 1
        )
    );
}


forAll(allRegions, i)
{
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
            (murs[i]-1)/(murs[i]*mu0)
        )
    );
}

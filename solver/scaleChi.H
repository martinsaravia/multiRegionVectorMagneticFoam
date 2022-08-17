// This code creates a new chi field on every Iteration since we need to scale with ppf
// the boundary condition since it is a function of chi. If we do not update chi
// here, then the boundary condition is always reading the full chi value (since
// it is created in createFields.H without scaling)

scalar ppf = 0.0; // Parameter Proportional Factor
scalar iStop(runTime.endTime().value()-chiIterCutOff); // Iteration of Stop
scalar iDelta(iStop-chiIterCutIn);

// Scale chi field with the iteration number in order to help convergence
if (runTime.value() > chiIterCutIn && runTime.value() <= iStop)
{
    ppf = 1.0/(iDelta)*runTime.value() - chiIterCutIn/(iDelta);
}
else if(runTime.value() > iStop)
{
    ppf = 1.0;
    scaleChi = false;
}

forAll(allRegions, i)
{
    volScalarField& chii = chis[i];
    chii = ppf * (murs[i]-1)/(murs[i]*mu0);
}

Info << "Scaling chi with factor: " << ppf << endl;

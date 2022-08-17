Info<< "--> Creating meshes..." << endl;
regionProperties rp(runTime);  // Create list of region on a per physics basis ("fluid", "solid")
const wordList allRegionsNames(rp["allRegions"]); // Works.

PtrList<fvMesh> allRegions(allRegionsNames.size());

forAll(allRegionsNames, i)
{
	Info<< "	Creating mesh for region: " << allRegionsNames[i]
		<< " for time = " << runTime.timeName()  << endl;

	allRegions.set
	(
		i,
		new fvMesh
		(
			IOobject
			(
				allRegionsNames[i],
				runTime.timeName(),
				runTime,
				IOobject::MUST_READ
			)
		)
	);
}

Info<< "--> Meshes created !\n" << endl;

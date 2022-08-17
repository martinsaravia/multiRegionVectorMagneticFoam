Info<< "Creating magnetic controls..." << endl;

// We need to create a dict from the fvSolution file because
// it is not possible to access to the fvSolution IOdictionary
// because the classes do not return a non-const reference.
IOdictionary solDict
(
	IOobject
	(
	"fvSolution",
	runTime.system(),
	runTime,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
	)
);

dictionary& vmfControl = solDict.subDict("VMF");
dictionary& solverControl = solDict.subDict("solvers").subDict("A");
//dictionary& simpleDict = solDict.subDict("SIMPLE");

Info<< "Magnetic controls created !\n" << endl;

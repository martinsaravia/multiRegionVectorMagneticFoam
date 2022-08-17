
// Notes: Trying several options in the RHS, the only stable are the options wich use
// a reconstructed B, i.e.
// 	fvc::div(chii*2.0*skew(fvc::reconstruct(fvc::snGrad(Ai) * allRegions[i].magSf())))
//  fvc::curl(chii * Bri)
// The following non-reconstruted option diverge, at leas for the DYnamag case
// 	fvc::curl(chii * Bi)
// 	fvc::div(chii*2.0*skew(fvc::grad(Ai)))
//2 * fvc::div(chii*skew(fvc::reconstruct(fvc::snGrad(Ai) * allRegions[i].magSf())))
// Relaxation is implemented but not used in the RHS currently
Bi.relax();
Bri.relax();

tmp<fvVectorMatrix> tAEqn
	(
	   fvm::laplacian(Ai)
	==
	   mu0*(
				- Ji
				- fvc::curl(chii * Bri)
				- fvc::curl(Mi)
		  )
	);

fvVectorMatrix& AEqn = tAEqn.ref();

AEqn.solve(solverControl);

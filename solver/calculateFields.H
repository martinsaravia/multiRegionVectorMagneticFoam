const surfaceScalarField& magSfi = allRegions[i].magSf();

// Fields in terms of non-reduced fields
Bi = fvc::curl(Ai);
//Bi = fvc::surfaceIntegrate(Sfi ^ fvc::interpolate(Ai));
MIi = chii * Bi;
MTi = Mi + MIi;
Hi = (1/mu0) * Bi - MTi;

// Fields in terms of reduced fields
Bri = 2.0*(*skew(fvc::reconstruct(fvc::snGrad(Ai) * magSfi)));
MIri = chii * Bri;
MTri = Mi + MIri;
Hri = (1/mu0) * Bri - MTri;

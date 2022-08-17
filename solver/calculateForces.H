
const surfaceVectorField& Sfi = allRegions[i].Sf();
surfaceVectorField nf(Sfi / magSfi);
const fvBoundaryMesh& boundary = allRegions[i].boundary();

// Kelvin force
Fki = mu0 * MTi & fvc::grad(Hi);  // Volumetric Force from H
Fk2i = (MTi & fvc::grad(Bi)) + (MTi ^ fvc::curl(Bi));  // Volumetric force from B
dimensionedVector FkT = fvc::domainIntegrate(Fki);  // Total Kelvin Force in the region
dimensionedVector Fk2T = fvc::domainIntegrate(Fk2i);  // Total Kelvin Force in the region
Info<< "Fk: " << FkT.value()  << endl;
Info<< "Fk2: " << Fk2T.value()  << endl;

// Stress tensor
Ti = mu0 * (Hi * Hi - 0.5 * (Hi & Hi) * T1); // Maxwell Stress Tensor
const volTensorField::Boundary& Tib= Ti.boundaryField(); // MST at the boundary
dimensionedVector divTi = fvc::domainIntegrate(fvc::div(Ti)); // Integral of the normal times the stress tensor divergence
Info<< "divT: "<< divTi.value()  << endl;

// Boundary pressure
volVectorField::Boundary& Fkbref = Fki.boundaryFieldRef();
volVectorField::Boundary& Fk2bref = Fk2i.boundaryFieldRef();
surfaceScalarField Mnfi(fvc::interpolate(MTi) & nf);  // Normal magnetization
surfaceScalarField Mtfi(mag(fvc::interpolate(MTi) - Mnfi * nf));  // Tangential magnetization
surfaceVectorField Fb(mu0 * 0.5 * pow(Mnfi , 2) * nf); // Boundary magnetic pressure
surfaceVectorField Fb2(-mu0 * 0.5 * pow(Mtfi , 2) * nf); // Boundary magnetic pressure
forAll(boundary, patchi)
{
    // Assign the calculated boundary values to the Kelvin forde boundary. I does not perturb
    // the internal field since we never use interpolation for the volume force.
    // This print does not give the correct force since it is parallelized (integrate
    // in paraview after reconstruction using only the surface (not the internal volume)
    Fkbref[patchi] == Fb.boundaryField()[patchi];
    Fk2bref[patchi] == Fb2.boundaryField()[patchi];
    Info<<"Fb: "<< boundary[patchi].name() << " is "<< sum(Fb.boundaryField()[patchi] * magSfi[patchi])<<endl;
    Info<<"Fb2: "<< boundary[patchi].name() << " is "<< sum(Fb2.boundaryField()[patchi] * magSfi[patchi])<<endl;
}


// Forces in terms of reduced fields

// Kelvin force
Fkri = mu0 * (MTri & fvc::grad(Hri));
Fk2ri = (MTri & fvc::grad(Bri)) + (MTri ^ fvc::curl(Bri));
dimensionedVector FkrT = fvc::domainIntegrate(Fkri);
dimensionedVector Fk2rT = fvc::domainIntegrate(Fk2ri);
Info<< "Fkr: " << FkrT.value()  << endl;
Info<< "Fk2r: " << Fk2rT.value()  << endl;

// Stress tensor
Tri = (1/mui) * (Bri * Bri - 0.5 *  (Bri & Bri) * T1);  // Maxwell Stress Tensor
const volTensorField::Boundary& Trib= Tri.boundaryField(); // MST at the boundary
dimensionedVector divTri = fvc::domainIntegrate(fvc::div(Tri));
Info<< "divTr: "  << divTri.value()  << endl;

// Boundary pressure
volVectorField::Boundary& Fkrbref = Fkri.boundaryFieldRef();
volVectorField::Boundary& Fk2rbref = Fk2ri.boundaryFieldRef();
surfaceScalarField Mnfri(fvc::interpolate(MTri) & nf);
surfaceScalarField Mtfri(mag(fvc::interpolate(MTri) - Mnfri * nf));
surfaceVectorField Fbr(mu0 * 0.5 * pow(Mnfri , 2) * nf); // Boundary magnetic pressure
surfaceVectorField Fb2r(-mu0 * 0.5 * pow(Mtfri , 2) * nf); // Boundary magnetic pressure
forAll(boundary, patchi)
{
    Fkrbref[patchi] == Fbr.boundaryField()[patchi];
    Fk2rbref[patchi] == Fb2r.boundaryField()[patchi];
    Info<<"Fbr: "<< boundary[patchi].name() << " is "<< sum(Fbr.boundaryField()[patchi] * magSfi[patchi])<<endl;
    Info<<"Fb2r: "<< boundary[patchi].name() << " is "<< sum(Fb2r.boundaryField()[patchi] * magSfi[patchi])<<endl;
}

// Lorentz force
Fli = mu0 * Ji ^ Hi;













// Brown boundary force (Method 1)
// const surfaceVectorField::Boundary& Sfbi = allRegions[i].Sf().boundaryField();
// const surfaceScalarField::Boundary& magSfbi = allRegions[i].magSf().boundaryField();
// Calculate the normal surface field





// Normal magnetization jump (assumming air is on one side)



// Loop through the patches of the force field to assig the correct value






// // Brown Boundary force (method 2)
// const surfaceVectorField::Boundary& Sfbi2 = allRegions[i].Sf().boundaryField();
// //const fvBoundaryMesh& boundary = allRegions[i].boundary();
// surfaceVectorField Mfi(fvc::interpolate(Mi+MIri));
// Fbi = fvc::interpolate(Mi+MIri);
// const polyBoundaryMesh& bMeshi = allRegions[i].boundaryMesh();
// forAll(boundary, patchi)
// {
//     scalarField Mni = Mfi.boundaryField()[patchi] & boundary[patchi].nf();
//     vectorField FB =  mu0.value() * 0.5 * pow(Mni, 2) * Sfbi2[patchi];
//     Info<<"Fb2: "<< boundary[patchi].name() << " is "<< sum(FB)<<endl;
// }
























//
//
// forAll(boundary, patchi)
// {
//     vectorField FT1 = Sfbi[patchi] & T1bi[patchi];
//     vectorField FT2 = Sfbi[patchi] & T2bi[patchi];
//     Info<<"T1Patch: "<< boundary[patchi].name() << " is "<< sum(FT1)<<endl;
//     Info<<"T2Patch: "<< boundary[patchi].name() << " is "<< sum(FT2)<<endl;
//     scalarField Mni = Mfi.boundaryField()[patchi] & boundary[patchi].nf();
//     vectorField FB =  mu0.value() * 0.5 * pow(Mni, 2) * Sfbi[patchi];
//     Info<<"Fb: "<< boundary[patchi].name() << " is "<< sum(FB)<<endl;
//     Info<<"MaxMfi: "<< boundary[patchi].name() << " is "<< max(Mfi)<<endl;
//     Info<<"MaxMni: "<< boundary[patchi].name() << " is "<< max(Mni)<<endl;
//     Info<<"MinMni: "<< boundary[patchi].name() << " is "<< min(Mni)<<endl;
// }


// Print the force information to the log


// dimensionedVector force7 = fvc::domainIntegrate(fvc::div(T3i));
// Info<< "The divergence of the MST for region "<< allRegions[i].name() <<" is: " << force7  << nl<< endl;


// Boundary force



//Fbi = fvc::interpolate(Mi+Mrii);
//
// Fbi = 0.5 * pow(Mfi.boundaryField() & nfbi, 2) * Sfbi;



//const polyBoundaryMesh& bMeshi = allRegions[i].boundaryMesh();



// dimensionedVector force3 = fvc::domainIntegrate(Fk3i);
// Info<< "Fk3: " << force3.value()  << endl;
// dimensionedVector force4 = fvc::domainIntegrate(Fk4i);
// Info<< "Fk4: " << force4.value()  << endl;

//Fk2i = mui * (Mrii + Mi) & GHT.T();
// Fk3i = mu0 * (Mii + Mi) & fvc::reconstruct(fvc::snGrad(Hi) * allRegions[i].magSf());
// Fk4i = mu0 * (Mrii + Mi) & GHTr.T();
//Fk3i = mu0 * Hi * fvc::div(Hi);
//Fk4i = fvc::curl(chii*Bi) ^ Bi;

// Stress tensors
// Ti = (1.0 / mu0) * (Bi * Bi - 0.5 * (Bi & Bi) * T1);
// T2i = (1.0 / mu0) * (Bri * Bri - 0.5 * (Bri & Bri) * T1) ;
// T2i = Bi * Hi - 0.5 * mu0 * (Hi & Hi) * T1 ;
// T2i = mu0 * (Mii + Mi + Hi) * Hi - 0.5 * mu0 * (Hi & Hi) * T1 ;
// T3i = (1.0 / (mu0*muri)) * Bi * Bi  - 0.5 * mu0 * (1.0 / (pow(mu0*muri, 2))) * (Bi & Bi) * T1 ;
// T3i = (1.0 / (mu0*muri)) * Bi * Bi  - 0.5 * mu0 * (1.0 / (pow(mu0*muri, 2))) * (Bi & Bi) * T1 ;
// Ti = Bi * Hi - 0.5 * mu0 * (Hi & Hi) * T1 ;
// T2i = Bri * Hri - 0.5 * mu0 * (Hri & Hri) * T1 ;

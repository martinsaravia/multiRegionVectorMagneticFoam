/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) Original authors
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    multiRegionVectorMagneticFoam

Description of version v1.0
    Magnetostatics solver based on the vector potential formulation using multi regions.
	This version uses two approaches, in one approach you can scale the chi field,
	in the other you can impose it at the fist iteration. Both schemes appear to work ok.
	The solver uses under-relaxation of the magnetic field in order to stabilizes
	the solution; othwerise oscilations appear and the solution may converge to a
	surprinsingly stable solution with opposite sign.

Authors
    Martin Saravia

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "OSspecific.H"
#include "electromagneticConstants.H"
#include "regionProperties.H" // Allows creation of object to hold region properties
//#include "fvOptions.H" //
#include "solidMagnetostaticModel.H"
#include "simpleControl.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
		#include "setRootCaseLists.H" // Add some list functionality to the old setRootCase
		#include "createTime.H"
		#include "createMeshes.H"
		#include "createFields.H"
		#include "createMagneticControls.H"
		#include "readMagneticControls.H"

		// Use the simple solver to simplify convergence control (only 1 control for all regions)
		fvMesh& mesh = allRegions[0];
		simpleControl simple(mesh);

		int nCorrNonOrth(simple.nCorrNonOrth());

		/* Store the initial iter of the magnetic field
		(probably could be replaced by something better) */
		forAll(allRegions, i)
		{
			Bs[i].storePrevIter();
			Brs[i].storePrevIter();
		}

		// Outer loop for ramping chi
		while (simple.loop(runTime))
	    {

			Info<< "Iteration: " << runTime.value() << nl << endl;

			// Scale chi
			if (scaleChi)
			{
				#include "scaleChi.H"
			}

			// Set the tolerance for the last two iterations (access through the mesh)
			if (correctFinal)
			{
				if (runTime.value() >= runTime.endTime().value() - finalIterations)
				{
					Info<< "Switching to final correction solution parameters..." << nl << endl;
					solverControl.set("tolerance", finalTolerance);
					// Attempt to modify the fvSolution dict and write it to file
					// but I can get it to write in the correct location (region/system)
					/* simpleDict.set("nNonOrthogonalCorrectors", 2);
					solDict.regIOobject::write(true); // Call the base class method */
					nCorrNonOrth = finalCorrectors;
					correctFinal = false;
				}
			}
			for (int k = 1; k<=nCorrNonOrth; k++)
			{
				// Solve regions
				Info<< "Non-orthogonal corrector: " << k << nl << endl;
				forAll(allRegions, i)
				{
		        	Info<< nl <<"--> Solving region: " << allRegions[i].name() << endl;
					#include "setRegionFields.H"
					#include "solveRegion.H"
					#include "calculateFields.H"
					if (calcForces)
					{
						Info<< nl <<"--> Calculating forces ..."<< endl;
						#include "calculateForces.H"
					}
				}
			}

			runTime.write();

			Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
				<< nl << endl;

	}
	return 0;
}

// Choose chi Scaling
Switch scaleChi = vmfControl.lookupOrDefault<Switch>("scaleChi", "on");
int chiIterCutIn = vmfControl.lookupOrDefault<int>("chiIterCutIn", 2);
int chiIterCutOff = vmfControl.lookupOrDefault<int>("chiIterCutOff", 5);
Switch calcForces = vmfControl.lookupOrDefault<Switch>("calcForces", "on");

// Last correction steps
Switch correctFinal= vmfControl.lookupOrDefault<Switch>("correctFinal", "on");
int finalIterations = vmfControl.lookupOrDefault<int>("finalIterations", 5);
scalar finalTolerance = vmfControl.lookupOrDefault<scalar>("finalTolerance", 1E-10);
int finalCorrectors= vmfControl.lookupOrDefault<int>("finalCorrectors", 50);

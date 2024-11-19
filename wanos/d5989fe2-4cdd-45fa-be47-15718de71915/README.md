# Conversion_CG-to-AA
**Conversion_CG-to-AA** WaNo converts the CG model from the 3Dprint WaNo to an all-atom representation.
The final result is a MD simualtion-ready, all-inclusive LAMMPS datafile for the polymer network. This will be later used to run the relaxation, as well as docking simulations.

## Input
The input required would be the reference CG conformation as a LAMMPSdata file. A sample is included as fRESIN_porous.

## Settings
No user-defined settings are needed.

## Output
The output of the WaNo is an all-atom representation of the polymer network built from the CG model, in LAMMPSdata format. The forcefield parameters are included to run the relaxation.

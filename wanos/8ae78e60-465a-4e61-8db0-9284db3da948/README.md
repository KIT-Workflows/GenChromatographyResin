# LMP-GENERIC
Lammps-generic WaNo for Simstack is a WaNo designed to perform molecular dynamics (MD) simulations with LAMMPS package.

LMP-GENERIC takes LAMMPS input files to run MD simulations optionally initialized by the input data file. Such LAMMPS data files can be generated by LAMMPS from alternative simulations or prepared by the lammps-interface (LMP-INTER) WaNo.

**lammps-generic** operates in two exclusive modes:
 - *basic scripting* is a step-by-step process where a LAMMPS input script file is generated using the simstack GUI tool (Work-in-progress)
 - *advanced scripting* takes readily available LAMMPS input files

The final result are the log file, the intermediate data file, and the trajectory file of the MD simualtion. These data files allow for quick analysis based on "compute", "variable" and "thermo" outputs from within the LAMMPS simulation, as well as visualization and post-processing based on the trajectory.

## Input
If needed, the **lammps-generic** WaNo takes a readily available LAMMPS data file as the initial state of the MD simulation.

In advanced scripting mode, the **lammps-generic** WaNo requires a prescripted LAMMPS input file.

## Output
The output of the **LMP-GENERIC** Wano is a series LAMMPS outputs including:
 - *fLMPdataMED* as an output data file representing the intermediate / final simulation status
 - *fLMPlog* as an output log file including all the thermo outputs from the simulation run
 - *fLMPtrj* in the LAMMPStrajectory format as combined snapshots from the MD simulation, to be read by the visualization tools such as OVITO or VMD. The "rerun" command in LAMMPS can be utilized to review all data generated by the simulation.

 > In advanced scripting mode the log, data, and trajectory files have to be specified by exact names in italic in the previous list when using "log", "write_data", and "dump" commands in LAMMPS, otherwise the intermediate results may not be properly saved. For example, the output trajectory file name in the "log" command in LAMMPS should be "log fLMPlog".

## Example
In the WaNo a sample input file named sampleinput starts a MD simualtion from a clean sheet (without specifying the input data file) and performs a short simulation of grain boundary movements for a monolayer Alumium atoms using the EAM forcefield. All log, data, and trajectory files are written by respective LAMMPS commands following the guidelines in the Output section.

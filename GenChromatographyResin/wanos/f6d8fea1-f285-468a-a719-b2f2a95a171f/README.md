# LMP-INTER
Lammps-interface WaNo for Simstack is a WaNo designed to serve the first part in the IL@MOF workflow.

LMP-INTER takes .cif files from the MOF database as a starting point. Such .cif file can be acquired across multiple sources as CoRE MOF or Materials project. As .cif descriptions usually do not correspond to MD ready atomic coordinates and partial charges, **openbabel** and **egulp** are used in serial to generate a unitcell for the MOF in P1 symmetry while partial charges are assigned using the QEq method.

**Lammps-interface** is then employed to generate from the .cif file with QEq charges a supercell of desired dimensions. Also force-field parameters are automatically determined in correspondance with either UFF or UFF4MOF force fields.

The final result is a MD simualtion-ready, all-inclusive LAMMPS datafile for the MOF. This will be later used to generate the IL@MOF system, or as a base for GCMC simulated adsorption of gas molecules.

## Input
The sole input required would be the MOF structure as a .cif file. A HKUST1.cif is included in the WaNo for demonstration purposes.

## Settings
There are checkbox switches for either populating the MOF unitcell by P1 representation or designate partial charges based on QEq methods. 
 - If the .cif already corresponds to a MOF unitcell understandable by **lammps-interface** then this is mandatory.
 - If the .cif file already includes charge information, then reassignement via QEq is not necessary. Alternatively, if the MOF is to be modelled as a charge-neutral framework, the QEq checkbox should also be de-selected.

Since for most MOF unitcells, the dimensions of the all-atom MD simulation typically lies within 10nm x 10nm x 10nm, the number of unitcells can be controlled by setting up the supercell parameters. For HKUST1 in P1 representation, a 3 x 3 x 3 supercell corresponds to around 16K atoms with a Cartesian simulation box of 8nm x 8nm x 8nm, where simulations around 1ns timescale can be completed within around 10hrs on a generic GPU.

## Output
The output of the LMP-INTER Wano is a LAMMPS datafile with complete information to run a MD simulation as *data.fMOF*.

Additionally, the force-field parameters are separately stored in *FFparam* for further processing in the IL@MOF workflow.

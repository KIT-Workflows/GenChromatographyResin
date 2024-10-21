# polynet_input
**polynet_input** is a WaNo which provides all data needed to run the 3D printing simulation of a polymer network.
The polymer network is printed using 3 ingredients: HEMA as the basic building block of the backbone in the form of HEMA-HEMA dimers, EGDMA as a crosslinking agent, and TRP on the pedant sites, branching off the HEMA.

## Input
The input files are HEMA, EGDMA and TRP data files parametrized with OPLS-AA from LigParGen.
In particular, the EGDMAhalf is generated to provide all-atom conformations in the all-atom conversion.

## Settings
No user-defined settings are needed.

## Output
The WaNo passes the required datafiles to other WaNos in the workflow.
# CG-Polymerisation
**CG-Polymerisation** WaNo for Simstack is a WaNo designed to coarse-grained (CG) molecular dynamics simulations for direct laser writing of a polymer network.
In the CG representation, (typically acrylate-based) monomers can be activated to enable initiation, propagation and quenching reactions as described in https://doi.org/10.1038/s41467-022-29749-9.

## Input
Number of monomers, crosslinkers and ligands.
Dimensions of the simulation box.

## Settings
A *boolRandomSeeds* click box can be enabled to apply time-evolved seeds for the random number generator. When switched on, the CG polymer networks should be unique across multiple runs.

## Output
The final result is a lammps data file named *fRESIN_porous*, which serves as a template for the all-atom conversion in the chromatography resin generating workflow.

# Simstack Workflow: GenChromatographyResin for generating atomistic models for chromatography resins

# Table of contents
- [Aims](#Aims)
- [Installation](#Installation)
- [How it works](#How-it-works)
- [Contacts](#Contacts)

# Aims
[Go to the top](#Table-of-contents)

This workflow generates the atomistic representation of a polymer netowrk through coarse-grained 3D printing protocol, followed by the CG-to-All-Atom conversion.

This workflow is used to generate polymer network structures in [Ballweg et al.](https://doi.org/10.1016/j.chroma.2024.465089).

> [!IMPORTANT]
> The coarse-grained simulation to generate the topology of the 3D printed polymer network requires modified LAMMPS compilation to include the capability for generating bond with a certain bond type with an offset timestep. Currently, such LAMMPS compilation is hosted in the shared software folder on the INT-Nano cluster in int.kit.edu. A later update is planned to provide a docker for the entire software environment.

# Installation
[Go to the top](#Table-of-contents)

### Requirements

In order to run this Simstack workflow the following software are needed:
- Simstack 2021-12-15 (or newer version)
- LAMMPS 5Jun19 - modified for 3D printing

> [!NOTE]
> Simstack client can be installed and setup by following the [Simstack documentation](https://simstack.readthedocs.io/en/latest/)

> [!NOTE]
> The original results in [Ballweg et al.](https://doi.org/10.1016/j.chroma.2024.465089) were obtained using the software versions: Simstack (2021-12-15) and LAMMPS (5Jun19-modified).

# How it works 
[Go to the top](#Table-of-contents)

### Generation of coarse-grained polymer network structure and conversion to atomistic representation: inputs/outputs

The GenChromatographyResin workflow consists of 5 Simstack WaNos in a linear sequence of:
 1. Input_Polymer-Building-Blocks
 2. CG-Polymerisation
 3. Conversion_CG-to-AA
 4. MD-Relaxation
 5. Output_AA-Structure

The intermediate outputs from the previous WaNo is directly taken as the input of the next, rendering an automated protocol to generate the all-atom structure of a polymer network.

As inputs,
-  the building block types can be selected from the dropdown menu and number of corresponding building blocks are defined by the user in WaNo **Input_Polymer-Builidng-Blocks**
-  the system dimensions are defined by the user in WaNo **CG-Polymerisation** in the unit of nanometers

As output,
- The all-atom polymer network structure is downloadable in MOL2 data formats, which can subsequently be used as input to Maestro Schrödinger software in the follow-through [KNIME workflow](https://github.com/Ahmedkhalil-Mama/KNIME-Workflow-Binding-Poses-and-Energies-Calculation-of-Linear-Peptides).

A complete example, starting with 690 DHPMA monomers, 260 EGDMA crosslinkers, and 500 Q ligands with 4 -CH2- spacer units, in a simulation box size of 15x15x24 cubic nanometers, is included with all intermediate results in the [Example folder](Example/). The final MOL2 format output of the GenChromatographyResin is located in the [data of WaNo **Output_AA-Structure**](Example/2-205-01-05-17h35m04s-GenChromatographyResin_1.2c/exec_directories/2025-01-13-06h54m09s-Output_AA-structure/fPolyREF_conv.mol2).

# Contacts
If you have any question regarding this Simstack workflow, you can contact us:
-  liu.modan@kit.edu
-  wolfgang.wenzel@kit.edu
-  matthias.franzreb@kit.edu

# Supported Software

DP-GEN provides interfaces to various molecular dynamics (MD) and first-principles (FP) calculation software packages. This page summarizes all the supported software and their key features in the DP-GEN workflow. When used with these software packages, they must build integration with DeePMD-kit.

The calculation results from first-principles software are stored in **dpdata** formats (see [dpdata formats](https://docs.deepmodeling.com/projects/dpdata/en/master/formats.html) for more information).

## Machine Learning Potentials (MLP)

Machine learning potential engines are specified using the `mlp_engine` parameter and are used in the training phase to construct deep potential models.

### DeePMD-kit

{dargs:argument}`mlp_engine <run_jdata/mlp_engine>: dp`

[DeePMD-kit](https://github.com/deepmodeling/deepmd-kit) is a package for constructing deep potential models for molecular dynamics simulations. DeePMD-kit is the only supported machine learning potential software in DP-GEN and serves as the foundation for all deep potential training and inference. It provides the neural network models that are used by MD engines for force field calculations during the exploration phase.

## Molecular Dynamics (MD) Engines

MD engines are used in the exploration phase (`01.model_devi`) to generate candidate structures for potential improvement. These engines are specified using the `model_devi_engine` parameter.

### LAMMPS

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>: lammps`

[LAMMPS](https://www.lammps.org/) (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics code with a focus on materials modeling. LAMMPS serves as the primary MD engine for structure exploration in DP-GEN and integrates seamlessly with DeePMD-kit through the USER-DEEPMD package. The integration allows LAMMPS to use deep potential models for highly accurate force field calculations during molecular dynamics simulations. LAMMPS simulations in DP-GEN can be automatically restarted if they fail, providing robust handling of long-running molecular dynamics calculations.

LAMMPS integrated with DeePMD-kit can be installed following the [easy install guide](https://docs.deepmodeling.com/projects/deepmd/en/stable/install/easy-install.html) or by [building from source](https://docs.deepmodeling.com/projects/deepmd/en/stable/install/install-lammps.html). When configuring LAMMPS, ensure that the DeePMD-kit plugin is properly installed and that the {dargs:argument}`command <run_mdata/model_devi/command>` in the machine file points to the correct LAMMPS executable with DeePMD support.

### Amber

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>: amber`

[Amber](https://ambermd.org/) is a suite of biomolecular simulation programs primarily used for biological systems. Amber in DP-GEN is designed for biomolecular simulations and integrates with DeePMD-kit models for enhanced force field accuracy in biological systems. This integration is particularly useful for complex biomolecular systems requiring high-precision force field calculations.

AmberTools needs to enable the DeePMD-kit integration - refer to [Amber documentation](https://ambermd.org/) for installation details. The {dargs:argument}`command <run_mdata/model_devi/command>` in the machine file should be set to the path to sander for proper execution.

### Calypso

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>: calypso`

[Calypso](http://www.calypso.cn/) is a crystal structure prediction software that uses particle swarm optimization algorithms. Calypso in DP-GEN automates crystal structure generation and exploration for materials discovery workflows. It interfaces with deep potential models to efficiently sample the configuration space of crystal structures.

Configuration requires proper setup of the Calypso executable and ensuring compatibility with the DeePMD-kit interface for structure evaluation.

### Gromacs

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>: gromacs`

[Gromacs](http://www.gromacs.org/) is a versatile package for molecular dynamics simulations, particularly popular for biochemical molecules. Gromacs integration with DP-GEN enables high-performance biomolecular system exploration using deep potential models. The integration provides enhanced sampling techniques optimized for biological systems.

Gromacs integrated with DeePMD-kit should be installed following the [installation guide](https://docs.deepmodeling.com/projects/deepmd/en/stable/install/install-gromacs.html). Ensure that the {dargs:argument}`command <run_mdata/model_devi/command>` in the machine file is properly configured for the Gromacs executable.

## First-Principles (FP) Calculation Software

FP software packages are used in the labeling phase (`02.fp`) to calculate accurate energies and forces for selected structures. These are specified using the `fp_style` parameter. The input/output format is **dpdata** format - refer to [dpdata formats](https://docs.deepmodeling.com/projects/dpdata/en/master/formats.html) for all supported formats.

### VASP

{dargs:argument}`fp_style <run_jdata/fp_style>: vasp`

[VASP](https://www.vasp.at/) (Vienna Ab initio Simulation Package) is a package for performing ab initio quantum mechanical calculations using density functional theory (DFT). VASP is widely used in materials science for high-accuracy energy and force calculations in DP-GEN workflows. It integrates seamlessly with DP-GEN through automated input file generation and output parsing for training data preparation.

When using VASP, ensure that the INCAR file specifies KSPACING and KGAMMA parameters, and properly configure the {dargs:argument}`fp_pp_path <run_jdata[fp_style=vasp]/fp_pp_path>` and {dargs:argument}`fp_pp_files <run_jdata[fp_style=vasp]/fp_pp_files>` for pseudopotential files. Optional Custodian integration can be enabled for enhanced error handling.

### Gaussian

{dargs:argument}`fp_style <run_jdata/fp_style>: gaussian`

[Gaussian](https://gaussian.com/) is a computational chemistry software package for electronic structure modeling. Gaussian in DP-GEN is primarily used for molecular system calculations and provides flexible input file generation for gas-phase and solution-phase modeling. The integration supports automatic setup of Gaussian input files with specified keywords and multiplicity settings.

When using Gaussian, the {dargs:argument}`command <run_mdata/fp/command>` in the machine file should be set to `g16 < input || :` to ensure proper execution and error handling. Configure the keywords parameter appropriately for your computational level of theory.

### CP2K

{dargs:argument}`fp_style <run_jdata/fp_style>: cp2k`

[CP2K](https://www.cp2k.org/) is a quantum chemistry and solid state physics software package for atomistic simulations. CP2K in DP-GEN enables large-scale DFT calculations using mixed Gaussian and plane-wave basis sets with linear scaling methods. It is particularly effective for combined QM/MM simulations and materials modeling.

Configure CP2K input templates carefully and ensure proper basis set and pseudopotential specifications for your system. The integration supports both molecular and materials systems.

### Siesta

{dargs:argument}`fp_style <run_jdata/fp_style>: siesta`

[SIESTA](https://departments.icmab.es/leem/siesta/) is an electronic structure code for large systems using density functional theory. SIESTA in DP-GEN is efficient for extended systems and large-scale materials calculations due to its localized atomic orbital basis sets and linear scaling algorithms. It provides reduced computational cost for calculations involving structures with thousands of atoms.

Ensure proper configuration of basis sets and exchange-correlation functionals in the SIESTA input files. The code is particularly suitable for calculations requiring good scaling with system size.

### ABACUS

{dargs:argument}`fp_style <run_jdata/fp_style>: abacus`

[ABACUS](https://abacus.ustc.edu.cn/) (Atomic-orbital Based Ab-initio Computation at UStc) is a DFT software package with support for both plane-wave and localized atomic orbital basis sets. ABACUS in DP-GEN provides versatile DFT calculations with excellent efficiency for various system types from molecules to materials. It offers optimization for high-performance computing and flexibility for different calculation requirements.

Configure the INPUT file properly with appropriate basis sets and k-point sampling. ABACUS supports both plane-wave and localized atomic orbital basis sets, providing flexibility for different calculation requirements.

### PWSCF (Quantum ESPRESSO)

{dargs:argument}`fp_style <run_jdata/fp_style>: pwscf`

[Quantum ESPRESSO](https://www.quantum-espresso.org/) is an integrated suite of open-source computer codes for electronic-structure calculations and materials modeling at the nanoscale. PWSCF (PWmat Self-Consistent Field) integration in DP-GEN provides reliable DFT calculations with extensive pseudopotential libraries and well-established computational methods. It serves as an excellent open-source alternative to proprietary DFT codes.

Ensure proper pseudopotential selection and k-point convergence for your calculations. The extensive documentation and large user community provide good support for troubleshooting.

### PWmat

{dargs:argument}`fp_style <run_jdata/fp_style>: pwmat`

[PWmat](http://www.pwmat.com/) is a plane-wave density functional theory package with advanced features for materials simulation. PWmat in DP-GEN enables efficient DFT calculations with GPU acceleration capabilities, providing high-performance implementation for materials property predictions. It is particularly effective for calculations requiring significant computational resources.

Configure GPU settings appropriately if available, and ensure proper setup of the PWmat executable and license. The GPU acceleration can significantly reduce calculation times for large systems.

### Amber DPRc

{dargs:argument}`fp_style <run_jdata/fp_style>: amber/diff`

Amber DPRc is a specialized interface for using Amber within the first-principles calculation workflow. This interface in DP-GEN enables biomolecular system labeling where sander is used to label the structure rather than perform energy calculations. It provides specialized capabilities for systems requiring Amber-based structural analysis and labeling.

This fp_style only supports being used with {dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>` set to `amber`. The {dargs:argument}`command <run_mdata/fp/command>` in the machine file should be set to the path to sander, and the [dpamber](https://github.com/njzjz/dpamber) package must be installed in the environment.

### Custom

{dargs:argument}`fp_style <run_jdata/fp_style>: custom`

The custom interface allows users to integrate their own first-principles calculation software with DP-GEN. This interface provides an extensible framework for integrating specialized or proprietary codes and enables research-specific calculation methods. It supports prototype testing of new software packages not directly supported by DP-GEN.

Users must provide input and output file format definitions and specify the custom script in the machine file. Forward and backward file handling must be properly defined for successful integration with the DP-GEN workflow.

For detailed configuration examples and parameters, refer to the [parameter documentation](../run/example-of-param.md) and [examples directory](https://github.com/deepmodeling/dpgen/tree/master/examples).
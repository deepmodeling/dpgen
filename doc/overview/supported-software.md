# Supported Software

DP-GEN provides interfaces to various molecular dynamics (MD) and first-principles (FP) calculation software packages. This page documents all supported software and their integration with the DP-GEN workflow. When used with DP-GEN, these software packages must build integration with DeePMD-kit.

Calculation results from first-principles software are stored in **dpdata** formats. See the [dpdata formats documentation](https://docs.deepmodeling.com/projects/dpdata/en/master/formats.html) for detailed format specifications.

For detailed configuration examples and parameters, refer to the [parameter documentation](../run/example-of-param.md) and [examples directory](https://github.com/deepmodeling/dpgen/tree/master/examples).

## Machine Learning Potentials (MLP)

Machine learning potential engines are specified using the `mlp_engine` parameter for the training phase.

### DeePMD-kit

{dargs:argument}`mlp_engine <run_jdata/mlp_engine>`: `dp`

[DeePMD-kit](https://github.com/deepmodeling/deepmd-kit) is a deep potential modeling package for molecular dynamics simulations. DeePMD-kit is the only supported machine learning potential software and serves as the foundation for all deep potential training and inference. The trained neural network models are used by molecular simulation engines for accurate force field calculations during structure exploration.

## Molecular Simulation Engines

Molecular simulation engines are used in the exploration phase (`01.model_devi`) to generate candidate structures. These engines are specified using the `model_devi_engine` parameter.

### LAMMPS

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>`: `lammps`

[LAMMPS](https://www.lammps.org/) (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics code focused on materials modeling. LAMMPS serves as the primary molecular simulation engine for structure exploration and integrates with DeePMD-kit through the USER-DEEPMD package. This integration enables accurate force field calculations using deep potential models during molecular dynamics simulations. LAMMPS simulations can be automatically restarted if they fail, providing robust handling of long-running calculations.

LAMMPS with DeePMD-kit integration can be installed via [easy install](https://docs.deepmodeling.com/projects/deepmd/en/stable/install/easy-install.html) or by [building from source](https://docs.deepmodeling.com/projects/deepmd/en/stable/install/install-lammps.html). Ensure the {dargs:argument}`command <run_mdata/model_devi/command>` in the machine file points to the correct LAMMPS executable with DeePMD support.

### Amber

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>`: `amber`

[Amber](https://ambermd.org/) is a suite of biomolecular simulation programs primarily used for biological systems. Amber enables biomolecular simulations with DeePMD-kit models for enhanced force field accuracy in biological systems. This integration is particularly useful for complex biomolecular systems requiring high-precision force calculations.

AmberTools must enable DeePMD-kit integration - refer to [Amber documentation](https://ambermd.org/) for installation details. Set the {dargs:argument}`command <run_mdata/model_devi/command>` in the machine file to the path to sander.

### Calypso

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>`: `calypso`

[Calypso](http://www.calypso.cn/) is a crystal structure prediction software using particle swarm optimization algorithms. Calypso automates crystal structure generation and exploration for materials discovery workflows. It interfaces with deep potential models to efficiently sample crystal structure configuration space.

Configure the Calypso executable properly and ensure compatibility with the DeePMD-kit interface for structure evaluation.

### Gromacs

{dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>`: `gromacs`

[Gromacs](http://www.gromacs.org/) is a versatile molecular dynamics package particularly popular for biochemical molecules. Gromacs enables high-performance biomolecular system exploration using deep potential models. The integration provides enhanced sampling techniques optimized for biological systems.

Install Gromacs with DeePMD-kit integration following the [installation guide](https://docs.deepmodeling.com/projects/deepmd/en/stable/install/install-gromacs.html). Configure the {dargs:argument}`command <run_mdata/model_devi/command>` in the machine file properly for the Gromacs executable.

## First-Principles (FP) Calculation Software

FP software packages are used in the labeling phase (`02.fp`) to calculate accurate energies and forces for selected structures. These are specified using the `fp_style` parameter. Input/output formats follow **dpdata** standards - see [dpdata formats](https://docs.deepmodeling.com/projects/dpdata/en/master/formats.html) for complete format specifications.

### VASP

{dargs:argument}`fp_style <run_jdata/fp_style>`: `vasp`

[VASP](https://www.vasp.at/) (Vienna Ab initio Simulation Package) performs ab initio quantum mechanical calculations using density functional theory (DFT). VASP is widely used in materials science for high-accuracy energy and force calculations. It integrates seamlessly through automated input file generation and output parsing for training data preparation.

Configure INCAR files with KSPACING and KGAMMA parameters, and properly set {dargs:argument}`fp_pp_path <run_jdata[fp_style=vasp]/fp_pp_path>` and {dargs:argument}`fp_pp_files <run_jdata[fp_style=vasp]/fp_pp_files>` for pseudopotential files. Optional Custodian integration can be enabled for enhanced error handling.

### Gaussian

{dargs:argument}`fp_style <run_jdata/fp_style>`: `gaussian`

[Gaussian](https://gaussian.com/) is a computational chemistry software package for electronic structure modeling. Gaussian is primarily used for molecular system calculations and provides flexible input file generation for gas-phase and solution-phase modeling. The integration supports automatic setup of input files with specified keywords and multiplicity settings.

Set the {dargs:argument}`command <run_mdata/fp/command>` in the machine file to `g16 < input || :` to ensure proper execution and error handling. Configure the keywords parameter appropriately for your computational level of theory.

### CP2K

{dargs:argument}`fp_style <run_jdata/fp_style>`: `cp2k`

[CP2K](https://www.cp2k.org/) is a quantum chemistry and solid state physics software package for atomistic simulations. CP2K enables large-scale DFT calculations using mixed Gaussian and plane-wave basis sets with linear scaling methods. It is particularly effective for combined QM/MM simulations and materials modeling.

Configure input templates carefully and ensure proper basis set and pseudopotential specifications. The integration supports both molecular and materials systems.

### SIESTA

{dargs:argument}`fp_style <run_jdata/fp_style>`: `siesta`

[SIESTA](https://siesta-project.org/siesta/) is an electronic structure code for large systems using density functional theory. SIESTA is efficient for extended systems and large-scale materials calculations due to its localized atomic orbital basis sets and linear scaling algorithms. It provides reduced computational cost for calculations involving structures with thousands of atoms.

Configure basis sets and exchange-correlation functionals properly in input files. SIESTA is particularly suitable for calculations requiring good scaling with system size.

### ABACUS

{dargs:argument}`fp_style <run_jdata/fp_style>`: `abacus`

[ABACUS](https://abacus.ustc.edu.cn/) (Atomic-orbital Based Ab-initio Computation at UStc) is a DFT software package with support for both plane-wave and localized atomic orbital basis sets. ABACUS provides versatile DFT calculations with excellent efficiency for various system types from molecules to materials. It offers optimization for high-performance computing and flexibility for different calculation requirements.

Configure the INPUT file properly with appropriate basis sets and k-point sampling. ABACUS supports both plane-wave and localized atomic orbital basis sets, providing flexibility for different calculation requirements.

### PWSCF (Quantum ESPRESSO)

{dargs:argument}`fp_style <run_jdata/fp_style>`: `pwscf`

[Quantum ESPRESSO](https://www.quantum-espresso.org/) is an integrated suite of open-source computer codes for electronic-structure calculations and materials modeling at the nanoscale. PWSCF (PWmat Self-Consistent Field) integration provides reliable DFT calculations with extensive pseudopotential libraries and well-established computational methods. It serves as an excellent open-source alternative to proprietary DFT codes.

Ensure proper pseudopotential selection and k-point convergence for calculations. The extensive documentation and large user community provide good support for troubleshooting.

### CPMD (Quantum ESPRESSO)

{dargs:argument}`fp_style <run_jdata/fp_style>`: `cpx`

CPMD (Car-Parrinello Molecular Dynamics) is another integration from [Quantum ESPRESSO](https://www.quantum-espresso.org/) suite, specifically the `cp.x` executable. This interface is well-suited for calculating energies and forces of various types of structures, even large molecules, using Car-Parrinello method. It utilizes a plane-wave basis set and pseudopotentials, optimized for efficient molecular dynamics simulations.

To use it with DP-GEN, one must prepare a template for the `cp.x` program configured for an electron relaxation procedure to obtain accurate ground-state properties, which is usually achieved after several tens to several hundred minimization steps. Be cautious to set `iprint` equal to `nstep` to ensure the output contains the final relaxed structure and its properties. Other key configuration requirements for the template are: both `ATOMIC_POSITIONS` and `CELL_PARAMETERS` cards must contain the `{angstrom}` keyword; the `%POSITIONS%` string must be in place of actual atomic coordinates, while the `%CELL%` string must be in place of actual cell vectors. It is also recommended to use absolute path for pseudopotentials. Ensure proper pseudopotential selection and electron dynamics configuration for reliable results.

### PWmat

{dargs:argument}`fp_style <run_jdata/fp_style>`: `pwmat`

[PWmat](http://www.pwmat.com/) is a plane-wave density functional theory package with advanced features for materials simulation. PWmat enables efficient DFT calculations with GPU acceleration capabilities, providing high-performance implementation for materials property predictions. It is particularly effective for calculations requiring significant computational resources.

Configure GPU settings appropriately if available, and ensure proper setup of the PWmat executable and license. GPU acceleration can significantly reduce calculation times for large systems.

### Amber DPRc

{dargs:argument}`fp_style <run_jdata/fp_style>`: `amber/diff`

Amber DPRc is a specialized interface for using Amber within the first-principles calculation workflow. This interface enables biomolecular system labeling where sander is used to label the structure rather than perform energy calculations. It provides specialized capabilities for systems requiring Amber-based structural analysis and labeling.

This fp_style only supports being used with {dargs:argument}`model_devi_engine <run_jdata/model_devi_engine>` set to `amber`. Set the {dargs:argument}`command <run_mdata/fp/command>` in the machine file to the path to sander, and ensure the [dpamber](https://github.com/njzjz/dpamber) package is installed in the environment.

### Custom

{dargs:argument}`fp_style <run_jdata/fp_style>`: `custom`

The custom interface allows users to integrate their own first-principles calculation software. This interface provides an extensible framework for integrating specialized or proprietary codes and enables research-specific calculation methods. It supports prototype testing of new software packages not directly supported.

Users must provide input and output file format definitions following [dpdata formats](https://docs.deepmodeling.com/projects/dpdata/en/master/formats.html) and specify the custom script in the machine file. Forward and backward file handling must be properly defined for successful integration with the workflow.

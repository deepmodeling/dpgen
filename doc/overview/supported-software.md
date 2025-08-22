# Supported Software

DP-GEN provides interfaces to various molecular dynamics (MD) and first-principles (FP) calculation software packages. This page summarizes all the supported software and their key features in the DP-GEN workflow.

## Molecular Dynamics (MD) Engines

MD engines are used in the exploration phase (`01.model_devi`) to generate candidate structures for potential improvement. These engines are specified using the `model_devi_engine` parameter.

### LAMMPS

[LAMMPS](https://www.lammps.org/) (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics code with a focus on materials modeling.

**Key Features:**
- High-performance parallel MD simulations
- Extensive force field support
- Deep potential integration via DeePMD-kit
- Flexible input scripting language
- Wide adoption in materials science community

**Usage in DP-GEN:**
- Primary MD engine for structure exploration
- Supports various ensemble types (NVT, NPT, NVE)
- Customizable simulation parameters
- Efficient scaling on HPC systems

### Amber

[Amber](https://ambermd.org/) is a suite of biomolecular simulation programs primarily used for biological systems.

**Key Features:**
- Specialized for biomolecular systems
- DPRc (Deep Potential Reactive) integration
- Advanced sampling methods
- Extensive force field libraries for biological molecules

**Usage in DP-GEN:**
- Biomolecular system exploration
- Used with DPRc models for reactive force fields
- Requires specific setup with dpamber package

### Calypso

[Calypso](http://www.calypso.cn/) is a crystal structure prediction software that uses particle swarm optimization algorithms.

**Key Features:**
- Crystal structure prediction
- Variable composition optimization
- Interface with various DFT codes
- Evolutionary algorithms for structure searching

**Usage in DP-GEN:**
- Automated structure generation
- Crystal structure exploration
- Materials discovery workflows

### Gromacs

[Gromacs](http://www.gromacs.org/) is a versatile package for molecular dynamics simulations, particularly popular for biochemical molecules.

**Key Features:**
- High-performance MD simulations
- Optimized for biological systems
- GPU acceleration support
- Advanced analysis tools
- Free and open-source

**Usage in DP-GEN:**
- Biomolecular system exploration
- Enhanced sampling techniques
- Integration with deep potential models

## First-Principles (FP) Calculation Software

FP software packages are used in the labeling phase (`02.fp`) to calculate accurate energies and forces for selected structures. These are specified using the `fp_style` parameter.

### VASP

[VASP](https://www.vasp.at/) (Vienna Ab initio Simulation Package) is a package for performing ab initio quantum mechanical calculations using density functional theory (DFT).

**Key Features:**
- Plane-wave basis set DFT calculations
- Hybrid functionals support
- Spin-orbit coupling
- Advanced electronic structure methods
- Widely used in materials science

**Usage in DP-GEN:**
- High-accuracy energy and force calculations
- Supports KSPACING and KGAMMA settings
- Pseudopotential file management
- Optional Custodian integration for error handling

### Gaussian

[Gaussian](https://gaussian.com/) is a computational chemistry software package for electronic structure modeling.

**Key Features:**
- Gaussian basis set calculations
- Wide range of quantum chemical methods
- Molecular property calculations
- Excited state calculations
- Thermochemical analysis

**Usage in DP-GEN:**
- Molecular system calculations
- Gas-phase and solution-phase modeling
- Flexible input file generation

### CP2K

[CP2K](https://www.cp2k.org/) is a quantum chemistry and solid state physics software package for atomistic simulations.

**Key Features:**
- Mixed Gaussian and plane-wave basis sets
- Linear scaling DFT methods
- Molecular dynamics capabilities
- Free and open-source
- Efficient parallel implementation

**Usage in DP-GEN:**
- Large-scale DFT calculations
- Combined QM/MM simulations
- Materials and molecular systems

### SIESTA

[SIESTA](https://departments.icmab.es/leem/siesta/) is an electronic structure code for large systems using density functional theory.

**Key Features:**
- Localized atomic orbital basis sets
- Linear scaling algorithms
- Large system calculations (1000+ atoms)
- Materials modeling focus
- Open-source availability

**Usage in DP-GEN:**
- Large-scale materials calculations
- Efficient for extended systems
- Reduced computational cost for big structures

### ABACUS

[ABACUS](https://abacus.ustc.edu.cn/) (Atomic-orbital Based Ab-initio Computation at UStc) is a DFT software package with support for both plane-wave and localized atomic orbital basis sets.

**Key Features:**
- Plane-wave and LCAO basis sets
- Advanced DFT implementations
- Machine learning integration
- High-performance computing optimization
- Open-source development

**Usage in DP-GEN:**
- Versatile DFT calculations
- Native machine learning potential support
- Efficient for various system types

### PWSCF (Quantum ESPRESSO)

[Quantum ESPRESSO](https://www.quantum-espresso.org/) is an integrated suite of open-source computer codes for electronic-structure calculations and materials modeling at the nanoscale.

**Key Features:**
- Plane-wave pseudopotential DFT
- Advanced electronic structure methods
- Materials modeling capabilities
- Open-source and freely available
- Large user community

**Usage in DP-GEN:**
- Reliable DFT calculations
- Extensive pseudopotential libraries
- Well-established computational methods

### PWmat

[PWmat](http://www.pwmat.com/) is a plane-wave density functional theory package with advanced features for materials simulation.

**Key Features:**
- Plane-wave DFT calculations
- GPU acceleration
- Advanced algorithms
- Materials science applications
- High-performance implementation

**Usage in DP-GEN:**
- Efficient DFT calculations
- GPU-accelerated computations
- Materials property predictions

### Amber DPRc

Amber DPRc is a specialized interface for using Deep Potential Reactive (DPRc) models within the Amber MD package for reactive force field simulations.

**Key Features:**
- Reactive force field capabilities
- Chemical bond breaking/forming
- Specialized for reactive systems
- Integration with Amber ecosystem

**Usage in DP-GEN:**
- Reactive system labeling
- Chemical reaction modeling
- Biomolecular reactive processes

**Note:** This fp_style only supports being used with `model_devi_engine` set to `amber`.

### Custom

The custom interface allows users to integrate their own first-principles calculation software with DP-GEN.

**Key Features:**
- User-defined input/output formats
- Flexible script integration
- Custom file handling
- Extensible framework

**Usage in DP-GEN:**
- Integration of specialized or proprietary codes
- Research-specific calculation methods
- Prototype testing of new software

**Requirements:**
- User must provide input and output file format definitions
- Custom script must be specified in machine file
- Forward and backward file handling definition

## Choosing the Right Software

### For MD Engines:

- **LAMMPS**: Best general-purpose choice for materials systems
- **Amber**: Recommended for biomolecular systems, especially with DPRc
- **Calypso**: Ideal for crystal structure prediction workflows
- **Gromacs**: Good for biomolecular systems with specific requirements

### For FP Calculations:

- **VASP**: Excellent for materials science applications, most widely tested
- **Quantum ESPRESSO (PWSCF)**: Good open-source alternative to VASP
- **CP2K**: Best for large systems requiring efficient scaling
- **ABACUS**: Good for mixed basis set calculations and ML integration
- **Gaussian**: Preferred for molecular systems and gas-phase calculations
- **SIESTA**: Efficient for very large extended systems
- **PWmat**: Good for GPU-accelerated calculations
- **Amber DPRc**: Essential for reactive biomolecular systems
- **Custom**: For specialized or proprietary calculation methods

## Configuration Notes

1. **Software Installation**: Ensure all selected software packages are properly installed and accessible on your computing environment.

2. **Machine Files**: Configure appropriate machine files with correct paths and resource allocations for each software.

3. **Input Files**: Prepare template input files and parameter sets appropriate for your system and calculation requirements.

4. **Compatibility**: Some combinations work better together (e.g., Amber + Amber DPRc for reactive biomolecular systems).

5. **Performance**: Consider computational cost and scaling when choosing software for large-scale calculations.

For detailed configuration examples and parameters, refer to the [parameter documentation](../run/example-of-param.md) and [examples directory](https://github.com/deepmodeling/dpgen/tree/master/examples).
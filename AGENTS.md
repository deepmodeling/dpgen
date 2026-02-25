# DP-GEN: Deep Potential Generator

DP-GEN is a Python scientific computing package for generating deep learning-based interatomic potential energy models. It integrates with HPC systems, molecular dynamics software (LAMMPS, Gromacs), and ab-initio calculation software (VASP, PWSCF, etc.).

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Bootstrap and Install

- **Environment Setup**:
  - Ensure Python 3.9+ is available: `python --version`
  - Create virtual environment: `python -m venv dpgen_env && source dpgen_env/bin/activate`
  - **Preferred**: Install with uv: `uv pip install -e .` -- takes 1-5 minutes due to scientific dependencies. Set timeout to 10+ minutes.
  - Alternative: `pip install -e .` (fallback if uv not available)
  - Released version: `uv pip install dpgen` or `pip install dpgen`
  - Test installation: `dpgen -h`

### Core Dependencies Installation Times

- **Dependencies**: Include large scientific packages (numpy, pymatgen, ASE, etc.)
- **Installation Time**: Full installation takes 1-5 minutes with uv, 3-10 minutes with pip. Use timeout of 10+ minutes.
- If network timeouts occur, retry with: `uv pip install --timeout 600 -e .` or `pip install --timeout 600 -e .`
- For faster testing: `uv pip install --no-deps -e .` or `pip install --no-deps -e .` (installs without dependencies, will fail at runtime)

### Testing

- **Unit Tests**: `python -m unittest discover tests -v` -- takes 2-5 minutes. Set timeout to 10+ minutes.
- **Coverage**: `coverage run --source=./dpgen -m unittest -v && coverage report` -- takes 3-8 minutes. Set timeout to 15+ minutes.
- **Quick Test**: `python -m unittest tests.test_load_file -v` -- takes \<1 second (will fail without dependencies but tests framework)
- **CLI Test**: `python -m unittest tests.test_cli -v` -- tests all dpgen subcommands

### Build Process

- **No Traditional Build**: This is a pure Python package using setuptools
- **Documentation Build**: `cd doc && make html` -- takes 2-5 minutes. Set timeout to 10+ minutes.
- **Package Build**: `python -m build` -- takes 1-3 minutes

## Development Requirements

### Commit Messages and PR Titles

- **Use Semantic Commit Messages**: Follow conventional commit format for all commits and PR titles
  - `feat:` for new features
  - `fix:` for bug fixes
  - `docs:` for documentation changes
  - `test:` for test additions/modifications
  - `refactor:` for code refactoring
  - `chore:` for maintenance tasks
  - Examples: `feat: add comprehensive GitHub Copilot instructions`, `fix: resolve timeout in dependency installation`

### Package Management

- **Prefer uv**: Use `uv` for Python dependency management when available
  - Installation: `uv pip install -e .`
  - Adding dependencies: `uv add package-name`
  - Fallback to `pip` only when `uv` is not available in the environment

## Main Workflows and Commands

### dpgen run - Core Workflow

- **Purpose**: Main deep potential generation with iterative training/exploration/labeling
- **Usage**: `dpgen run param.json machine.json`
- **Key Files**:
  - `param.json`: System parameters, training config, exploration settings
  - `machine.json`: HPC/computation resource configuration
- **Process**: Creates iterations (iter.000001, iter.000002, etc.) with 00.train, 01.model_devi, 02.fp subdirectories

### dpgen init\_\* - Data Preparation

- **init_bulk**: `dpgen init_bulk param.json [machine.json]` -- bulk systems
- **init_surf**: `dpgen init_surf param.json [machine.json]` -- surface systems
- **init_reaction**: `dpgen init_reaction param.json [machine.json]` -- reactive systems
- **Purpose**: Generate initial training data for deep potential models

### dpgen autotest - Validation

- **Usage**: `dpgen autotest make|run|post param.json [machine.json]`
- **Three Phases**:
  - `make`: Set up calculation tasks
  - `run`: Execute calculations
  - `post`: Analyze results
- **Supports**: VASP, ABACUS, DeepMD, MEAM, EAM force fields

### dpgen simplify - Dataset Optimization

- **Usage**: `dpgen simplify param.json machine.json`
- **Purpose**: Reduce dataset size while maintaining quality

### Other Commands

- **collect**: `dpgen collect JOB_DIR OUTPUT` -- gather generated data
- **db**: `dpgen db param.json` -- database operations
- **gui**: `dpgen gui` -- web interface (requires dpgui package)

## Configuration Files

### Parameter Files (param.json)

- **Location**: `examples/` directory contains templates for different scenarios
- **Key Examples**:
  - `examples/run/deprecated/dp2.x-lammps-vasp/param.json` -- VASP with LAMMPS
  - `examples/run/dp2.x-lammps-gaussian/param.json` -- Gaussian calculations
  - `examples/run/dp-calypso-vasp/param.json` -- CALYPSO with VASP
  - `examples/run/ch4/param.json` -- CH4 system example
- **Format**: JSON with extensive nested parameters for system, training, exploration settings

### Machine Files (machine.json)

- **Location**: `examples/machine/` directory
- **Purpose**: Define computational resources (HPC systems, cloud platforms)
- **Examples**:
  - `examples/machine/slurm/` -- Slurm cluster configurations
  - `examples/machine/pbs/` -- PBS system setup
  - `examples/machine/lebesgue/` -- Cloud platform integration

## Validation

### Always Run Before Committing

- **Full Dependencies Required**: All validation requires `uv pip install -e .` or `pip install -e .` (1-5 min install)
- `python -m unittest discover tests -v` -- full test suite (requires dependencies)
- `dpgen -h && dpgen run -h && dpgen autotest -h` -- CLI validation (requires installation)

### Working Without Full Installation

- **Code Exploration**: All source code can be read and edited without installation
- **Local Import Test**: `python -c "import sys; sys.path.insert(0, '.'); import dpgen"` (0.03s)
- **JSON Examples**: Configuration files in `examples/` can be examined without running
- **Test Framework**: `time python -m unittest tests.test_load_file -v` shows test structure (0.15s, will fail without deps)
- **File Structure**: Use `find examples -name "*.json"` to explore configuration templates

## Repository Structure

### Key Directories

- `dpgen/` -- Main source code
  - `generator/` -- Core DP-GEN functionality (`dpgen run`)
  - `auto_test/` -- Testing framework (`dpgen autotest`)
  - `data/` -- Data preparation (`dpgen init_*`)
  - `simplify/` -- Dataset optimization
  - `main.py` -- CLI entry point
- `tests/` -- Unit tests organized by module
- `examples/` -- Configuration templates and examples
- `doc/` -- Sphinx documentation

### Important Files

- `pyproject.toml` -- Modern Python packaging configuration
- `.github/workflows/test.yml` -- CI pipeline (Python 3.9, 3.12)
- `dpgen/main.py` -- CLI command definitions and entry points
- `examples/init/surf.json` -- Surface initialization example
- `examples/run/ch4/param.json` -- Complete run parameter example

## Common Development Tasks

### Adding New Features

- **Main CLI**: Edit `dpgen/main.py` to add new subcommands
- **Core Logic**: Extend modules in `dpgen/generator/`, `dpgen/auto_test/`, etc.
- **Tests**: Add corresponding tests in `tests/` directory
- **Examples**: Provide configuration examples in `examples/`

### Code Quality

- **Linting**: Uses ruff configuration in pyproject.toml
- **Import Style**: isort with black profile
- **Documentation**: Numpy-style docstrings, Sphinx for docs

### Working with Scientific Dependencies

- **Heavy Dependencies**: Package integrates with VASP, LAMMPS, Gaussian, etc.
- **File I/O**: Uses dpdata for structure file handling
- **HPC Integration**: dpdispatcher for job submission/management
- **Configuration**: Complex JSON schemas for scientific computing workflows

## Timing Expectations

- **Installation**: 1-5 minutes with uv, 3-10 minutes with pip (timeout 10+ minutes)
- **Full Test Suite**: 2-5 minutes (timeout 10+ minutes)
- **Documentation Build**: 2-5 minutes (timeout 10+ minutes)
- **Quick Import Test**: \<1 second
- **CI Pipeline**: Runs on Python 3.9 and 3.12, full cycle ~10-15 minutes

## Error Patterns

### Common Installation Issues

- **Network Timeouts**: Scientific packages are large, increase pip timeout
- **Missing System Dependencies**: Some packages require system libraries
- **Version Conflicts**: Pin to specific Python versions (3.9-3.12)

### Runtime Issues

- **Missing Dependencies**: Many features require specific scientific software
- **Configuration Errors**: JSON files have complex nested schemas
- **HPC Connectivity**: Machine files require proper authentication setup

Always test with minimal examples from `examples/` directory before using custom configurations.

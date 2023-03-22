## Autotest Overview: Autotest for Deep Generator
Suppose that we have a potential (can be DFT, DP, MEAM ...), `autotest` helps us automatically calculate M properties on N configurations. The folder where the `autotest` runs is called the working directory of `autotest`. Different potentials should be tested in different working directories.

A property is tested in three steps: `make`, `run` and `post`. `make` prepares all computational tasks that are needed to calculate the property. For example to calculate EOS, `make` prepares a series of tasks, each of which has a scaled configuration with certain volume, and all necessary input files necessary for starting a VASP, ABACUS, or LAMMPS calculations. `run` sends all the computational tasks to remote computational resources defined in a machine configuration file like `machine.json`, and automatically collects the results when remote calculations finish. `post` calculates the desired property from the collected results.

### Relaxation

The relaxation of a structure should be carried out before calculating all other properties:
```bash
dpgen autotest make relax.json
dpgen autotest run relax.json machine.json
dpgen autotest post relax.json
```
If, for some reasons, the main program terminated at stage `run`, one can easily restart with the same command.
`relax.json` is the parameter file. An example for `deepmd` relaxation is given as:
```json
{
        "structures":   ["confs/mp-*"],
        "interaction": {
                "type":         "deepmd",
                "model":        "frozen_model.pb",
                "type_map":     {"Al": 0, "Mg": 1}
        },
        "relaxation": {}
}
```

where the key `structures` provides the structures to relax. `interaction` is provided with `deepmd`, and other options are `vasp`, `abacus`, `meam`...

### Task type
There are now six task types implemented in the package: `vasp`, `abacus`, `deepmd`, `meam`, `eam_fs`, and `eam_alloy`. An `inter.json` file in json format containing the interaction parameters will be written in the directory of each task after `make`. We give input examples of the `interaction` part for each type below:

**VASP**:

The default of `potcar_prefix` is "".
```json
	"interaction": {
		"type":		"vasp",
		"incar":	"vasp_input/INCAR",
		"potcar_prefix":"vasp_input",
		"potcars":	{"Al": "POTCAR.al", "Mg": "POTCAR.mg"}
	}
```
**ABACUS**:

The default of `potcar_prefix` is "". The path of potcars/orb_files/deepks_desc is `potcar_prefix` + `potcars`/`orb_files`/`deepks_desc`/`deepks_model`.
```json
	"interaction": {
		"type":		"abacus",
		"incar":	"abacus_input/INPUT",
		"potcar_prefix":"abacus_input",
		"potcars":	{"Al": "pseudo_potential.al", "Mg": "pseudo_potential.mg"},
		"orb_files": {"Al": "numerical_orb.al", "Mg": "numerical_orb.mg"},
		"atom_masses": {"Al": 26.9815, "Mg":24.305},
		"deepks_desc": "jle.orb",
		"deepks_model": "model.ptg"
	}
```
**deepmd**:

**Only 1** model can be used in autotest in one working directory.

```json
	"interaction": {
		"type":		 "deepmd",
		"model":	 "frozen_model.pb",
		"type_map":      {"Al": 0, "Mg": 1}
	}
```
**meam**:

Please make sure the [USER-MEAMC package](https://lammps.sandia.gov/doc/Packages_details.html#pkg-user-meamc) has already been installed in LAMMPS.
```json
	"interaction": {
		"type":		 "meam",
		"model":	 ["meam.lib","AlMg.meam"],
		"type_map":      {"Al": 1, "Mg": 2}
	}
```
**eam_fs & eam_alloy**:

Please make sure the [MANYBODY package](https://lammps.sandia.gov/doc/Packages_details.html#pkg-manybody) has already been installed in LAMMPS
```json
	"interaction": {
		"type":		 "eam_fs (eam_alloy)",
		"model":	 "AlMg.eam.fs (AlMg.eam.alloy)",
		"type_map":      {"Al": 1, "Mg": 2}
	}
```

### Property type

Now the supported property types are `eos`, `elastic`, `vacancy`, `interstitial`, `surface`, and `gamma`. Before property tests, `relaxation` should be done first or the relaxation results should be present in the corresponding directory `confs/mp-*/relaxation/relax_task`. A file named `task.json` in json format containing the property parameter will be written in the directory of each task after `make` step. Multiple property tests can be performed simultaneously.

## Make run and post

There are three operations in auto test package, namely `make`, `run`, and `post`. Here we take `eos` property as an example for property type.

### Make
The `INCAR`, `POSCAR`, `POTCAR` input files for VASP or `in.lammps`, `conf.lmp`, and the interatomic potential files for LAMMPS will be generated in the directory `confs/mp-*/relaxation/relax_task` for relaxation or `confs/mp-*/eos_00/task.[0-9]*[0-9]` for EOS. The `machine.json` file is not needed for `make`. Example:
```bash
dpgen autotest make relaxation.json
```

### Run
The jobs would be dispatched according to the parameter in `machine.json` file and the calculation results would be sent back. Example:
```bash
dpgen autotest run relaxation.json machine.json
```

### Post
The post process of calculation results would be performed. `result.json` in json format will be generated in `confs/mp-*/relaxation/relax_task` for relaxation and `result.json` in json format and `result.out` in txt format in `confs/mp-*/eos_00` for EOS. The `machine.json` file is also not needed for `post`. Example:
```bash
dpgen autotest post relaxation.json
```

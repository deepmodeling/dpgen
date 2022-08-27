## Task-type

There are now five task types implemented in the package: `vasp`, `deepmd`, `meam`, `eam_fs`, and `eam_alloy`. An `inter.json` file in json format containing the interaction parameters will be written in the directory of each task. The input examples of the `"interaction"` part of each type can be found below:

### VASP: 
    
The default of `potcar_prefix` is "".
```json
	"interaction": {
		"type":		"vasp",
		"incar":	"vasp_input/INCAR",
		"potcar_prefix":"vasp_input",
		"potcars":	{"Al": "POTCAR.al", "Mg": "POTCAR.mg"}
	}
```
### deepmd:

**Only 1** model can be used in autotest in one working directory and the default `"deepmd_version"` is **1.2.0**.

```json
	"interaction": {
		"type":		 "deepmd",
		"model":	 "frozen_model.pb", 
		"type_map":      {"Al": 0, "Mg": 1},
                "deepmd_version":"1.2.0"
	}
```
### meam:
Please make sure the [USER-MEAMC package](https://lammps.sandia.gov/doc/Packages_details.html#pkg-user-meamc) has already been installed in LAMMPS.
```json
	"interaction": {
		"type":		 "meam",
		"model":	 ["meam.lib","AlMg.meam"],
		"type_map":      {"Al": 1, "Mg": 2}
	}
```
### eam_fs & eam_alloy:
Please make sure the [MANYBODY package](https://lammps.sandia.gov/doc/Packages_details.html#pkg-manybody) has already been installed in LAMMPS
```json
	"interaction": {
		"type":		 "eam_fs (eam_alloy)", 
		"model":	 "AlMg.eam.fs (AlMg.eam.alloy)", 
		"type_map":      {"Al": 1, "Mg": 2}
	}
```
 

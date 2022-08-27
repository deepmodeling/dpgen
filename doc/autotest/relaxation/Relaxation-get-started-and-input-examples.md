## Relaxation get started and input examples

The relaxation of a structure should be carried out before calculating all other properties. 

First, we need input parameter file and we name it `relax.json` here. All the relaxation calculations should be taken either by `VASP`, `ABACUS`, or `LAMMPS`. Here are two input examples for `VASP` and `LAMMPS` respectively. 

An example of the input file for relaxation by VASP:

```json
{
    "structures":            ["confs/std-*"],
    "interaction": {
            "type":           "vasp",
            "incar":          "vasp_input/INCAR",
            "potcar_prefix":  "vasp_input",
            "potcars":       {"Al": "POTCAR.al"}
	},
    "relaxation": {
            "cal_type":       "relaxation",
            "cal_setting":   {"relax_pos":       true,
                              "relax_shape":     true,
                              "relax_vol":       true,
                              "ediff":           1e-6,
                              "ediffg":         -0.01,
                              "encut":           650,
                              "kspacing":        0.1,
                              "kgamma":          false}
	}
}
```

Key words | data structure | example | description
---|---|---|---
**structures** | List of String | ["confs/std-*"] | path of different structures
**interaction** | Dict | See above | description of the task type and atomic interaction
**type** | String | "vasp" | task type
**incar** | String | "vasp_input/INCAR" | path for INCAR file in vasp
potcar_prefix | String | "vasp_input" | prefix of path for POTCAR file in vasp, default = ""
**potcars** | Dict | {"Al": "POTCAR.al"} | key is element type and value is potcar name
**relaxation** | Dict | See above | calculation type and setting for relaxation
cal_type  | String | "relaxation" or "static" | calculation type
cal_setting | Dict | See above | calculation setting
relax_pos | Boolean | true | relax atomic position or not, default = true for relaxation
relax_shape | Boolean | true | relax box shape or not, default = true for relaxation
relax_vol | Boolean | true | relax box volume or not, default = true for relaxation
ediff | Float | 1e-6 | set `EDIFF` parameter in INCAR files
ediffg | Float | -0.01 | set `EDIFFG` parameter in INCAR files
encut | Int | 650 | set `encut` parameter in INCAR files
kspacing | Float | 0.1 | set `KSPACING` parameter in INCAR files
kgamma | Boolean | false | set `KGAMMA` parameter in INCAR files

An example of the input file for relaxation by LAMMPS:

```json
{
    "structures":         ["confs/std-*"],
    "interaction": {
            "type":        "deepmd",
            "model":       "frozen_model.pb",
            "in_lammps":   "lammps_input/in.lammps",
            "type_map":   {"Al": 0}
	},
    "relaxation": {
            "cal_setting":{"etol": 0,
                           "ftol": 1e-10,
                           "maxiter": 5000,
                           "maximal": 500000}
	}
}
```
**Other key words different from vasp:**

Key words | data structure | example | description
---|---|---|---
**model** | String or List of String | "frozen_model.pb" | model file for atomic interaction
in_lammps | String | "lammps_input/in.lammps" | input file for lammps commands
**type_map** | Dict | {"Al": 0} | key is element type and value is type number. DP starts from 0, others starts from 1
etol | Float | 0 | stopping tolerance for energy
ftol | Float | 1e-10 | stopping tolerance for force
maxiter | Int | 5000 | max iterations of minimizer
maxeval | Int | 500000 | max number of force/energy evaluations

For LAMMPS relaxation and all the property calculations, **package will help to generate `in.lammps` file for user automatically** according to the property type. We can also make the final changes in the `minimize` setting (`minimize etol ftol maxiter maxeval`) in `in.lammps`. In addition, users can apply the input file for lammps commands in the `interaction` part. For further information of the LAMMPS relaxation, we refer users to [minimize command](https://lammps.sandia.gov/doc/minimize.html).



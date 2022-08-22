In some cases, we want to refine the calculation results of a property based on previous results by using different convergence criteria like `EDIFF` and `EDIFFG` or higher `ENCUT`. If the parameter of `init_from_suffix` and `output_suffix` are both provided in the input file, `refine` would start based on the results in `init_from_suffix` directory and output the results to `output_suffix` directory. Otherwise, the calculation results would be output to the default suffix `00`. An example of the input file is given below:
```json
{
    "structures":       ["confs/std-*"],
    "interaction": {
        "type":          "deepmd",
        "model":         "frozen_model.pb",
        "deepmd_version":"1.2.0",
        "type_map":     {"Al": 0}
    },
    "properties": [
        {
        "type":             "vacancy",
        "init_from_suffix": "00",
        "output_suffix":    "01",
        "cal_setting":     {"input_prop":  "lammps_input/lammps_high"}
        }
        ]
}
```

In this example, `refine` would output the results to `vacancy_01` based on the previous results in `vacancy_00` by using a different input commands file for lammps.
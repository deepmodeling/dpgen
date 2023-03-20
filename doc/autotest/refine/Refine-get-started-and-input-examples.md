## Refine get started and input examples

Sometimes we want to refine the calculation of a property from previous results. For example, when higher convergence criteria `EDIFF` and `EDIFFG` are necessary in VASP, the new VASP calculation is desired to start from the previous output configuration, rather than starting from scratch.


An example of the input file `refine.json` is given below:

```json
{
    "structures":       ["confs/std-*"],
    "interaction": {
        "type":          "deepmd",
        "model":         "frozen_model.pb",
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

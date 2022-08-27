## Reproduce-get-started-and-input-examples

Some times we want to reproduce the initial results with the same configurations for cross validation. This version of auto-test package can accomplish this successfully in all property types except for `Elastic`. An input example for using `deepmd` to reproduce the `VASP` Interstitial results is given as below:
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
        "type":             "interstitial",
        "reproduce":        true,
        "init_from_suffix": "00",
        "init_data_path":   "../vasp/confs",
        "reprod_last_frame":       false
        }
        ]
}
```

`reproduce` denotes whether to do `reproduce` or not and the default value is False. 

`init_data_path` is the path of VASP or LAMMPS initial data to be reproduced. `init_from_suffix` is the suffix of the initial data and the default value is "00". In this case, the VASP Interstitial results are stored in `../vasp/confs/std-*/interstitial_00` and the reproduced Interstitial results would be in `deepmd/confs/std-*/interstitial_reprod`. 

`reprod_last_frame` denotes if only the last frame is used in reproduce. The default value is True for eos and surface, but is False for vacancy and interstitial.

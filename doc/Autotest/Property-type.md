Now the supported property types are `eos`, `elastic`, `vacancy`, `interstitial`, and `surface`. Before property tests, `relaxation` should be done first or the relaxation results should be present in the corresponding directory `confs/mp-*/relaxation/relax_task`. A file named `task.json` in json format containing the property parameter will be written in the directory of each task. Multiple property tests can be performed simultaneously and are written in the `"properties"` part of the input file. An example of `EOS` and `Elastic` tests can be given as follows (please refer to [Property](https://github.com/deepmodeling/dpgen/wiki/Property:-get-started-and-input-examples) for further information of the property parameters):
```json
"properties": [
		{
                        "type":         "eos",
			"vol_start":    0.8,
			"vol_end":	1.2,
			"vol_step":	0.01
		},
		{
                        "type":         "elastic",
			"norm_deform":	2e-2,
			"shear_deform": 5e-2
		}
        ]
```
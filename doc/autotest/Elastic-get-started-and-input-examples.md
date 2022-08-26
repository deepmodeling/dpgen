## Elastic-get-started-and-input-examples

Here we calculate the mechanical properties which include elastic constants (C11 to C66), bulk modulus Bv, shear modulus Gv, Youngs modulus Ev, and Poission ratio Uv of a certain crystal structure.

#### An example of the input file for Elastic by deepmd:

```json
{
	"structures":	["confs/mp-*","confs/std-*","confs/test-*"],
	"interaction": {
		"type":		"deepmd",
                "model":        "frozen_model.pb",
		"type_map":	{"Al": 0, "Mg": 1}
	},
	"properties": [
            {
                "type":         "elastic",
                "norm_deform":	2e-2,
	       "shear_deform":  5e-2
	    }
        ]
}
```

Here the default values of `norm_deform` and `shear_deform` are **2e-3** and **5e-3**, respectively. A list of `norm_strains` and `shear_strains` would be generated as below:

```bash
[-norm_def, -0.5 * norm_def, 0.5 * norm_def, norm_def]
[-shear_def, -0.5 * shear_def, 0.5 * shear_def, shear_def]
```

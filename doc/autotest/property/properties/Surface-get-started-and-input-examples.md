## Surface get started and input examples

`Surface` calculates the surface energy. We need to give the information of `min_slab_size`, `min_vacuum_size`, `max_miller` (default value is 2), and `pert_xz` which means perturbations in xz and will help work around vasp bug.

#### An example of the input file for Surface by deepmd:

```json
{
	"structures":	"confs/mp-*",
	"interaction": {
		"type":		"deepmd",
                "model":        "frozen_model.pb",
		"type_map":	{"Al": 0, "Mg": 1}
	},
	"properties": [
            {
                "type":           "surface",
                "min_slab_size":  10,
                "min_vacuum_size":11,
                "max_miller":     2,
                "cal_type":       "static"
	    }
        ]
}
```

`Surface` calculates the surface energy. We need to give the information of `min_slab_size`, `min_vacuum_size`, `max_miller` (default value is 2), and `pert_xz` which means perturbations in xz and will help work around vasp bug. If `static-opt` parameter is given and is `True`, the static calculations of surface energies would be performed.

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
                "pert_xz":        0.01,
                "max_miller":     1,
                "static-opt":     false 
	    }
        ]
}
```

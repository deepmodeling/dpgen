`Vacancy` calculates the energy difference when removing an atom from the crystal structure. We only need to give the information of `supercell` to help calculate the vacancy energy and the default value of `supercell` is [1, 1, 1].

#### An example of the input file for Vacancy by deepmd:

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
                "type":         "vacancy",
                "supercell":	[1, 1, 1]
	    }
        ]
}
```
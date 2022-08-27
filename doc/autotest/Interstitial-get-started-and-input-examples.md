## Interstitial-get-started-and-input-examples

`Interstitial` calculates the energy difference when adding an atom into the crystal structure. We need to give the information of `supercell` (default value is [1, 1, 1]) and `insert_ele` list for the element types of the atoms added in.

#### An example of the input file for Interstitial by deepmd:

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
                "type":         "interstitial",
                "supercell":	[1, 1, 1],
                "insert_ele":   ["Al","Mg"],
                "conf_filters": {"min_dist": 1.5} 
	    }
        ]
}
```

We add a `conf_filters` parameter in `properties` part and this parameter can help to eliminate undesirable structure which can render rather difficult convergence in calculations. In the example above, **"min_dist": 1.5** means if the smallest atomic distance in the structure is less than 1.5 angstrom, the configuration would be eliminated and not used in calculations.

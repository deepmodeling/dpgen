## EOS get started and input examples

Equation of State (EOS) here calculates the energies of the most stable structures as a function of volume. Users can refer to Figure 4 of the [dpgen CPC paper](https://www.sciencedirect.com/science/article/pii/S001046552030045X?via%3Dihub) for more information of EOS.

#### An example of the input file for EOS by VASP:

```json
{
	"structures":	["confs/mp-*","confs/std-*","confs/test-*"],
	"interaction": {
		"type":		"vasp",
		"incar":	"vasp_input/INCAR",
                "potcar_prefix":"vasp_input",
		"potcars":	{"Al": "POTCAR.al", "Mg": "POTCAR.mg"}
	},
	"properties": [
		{
                        "type":         "eos",
			"vol_start":	10,
			"vol_end":	30,
			"vol_step":	0.5,
                        "change_box":   true
		}
        ]
}
```

`vol_start` is the starting volume per atom in Å^3/atom, `vol_step` is the increasing step of volume and the biggest volume is smaller than `vol_end`. In the above example, 40 tasks would be generated as `task.000000` to `task.000039` with the volume `10.00, 10.50, 11.00, ..., 29.50` Å^3/atom, respectively.

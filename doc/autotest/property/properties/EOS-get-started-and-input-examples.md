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
         "vol_start":    0.9,
         "vol_end":      1.1,
         "vol_step":     0.01
        }
        ]
}
```

`vol_start` is the starting relative volume per atom in Å^3/atom, `vol_step` is the increasing step of the ratio of relative volume, and the biggest volume is smaller than `vol_end` times the volume of equilibrium structure.

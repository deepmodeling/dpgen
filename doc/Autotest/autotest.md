Suppose that we have a potential (can be DFT, DP, MEAM ...), autotest helps us automatically calculate M porperties on N configurations. The folder where the autotest runs is called the autotest's working directory. Different potentials should be tested in different working directories.

A property is tested in three stages: make, run and post. make prepare all computational tasks that are needed to calculate the property. For example to calculate EOS, autotest prepare a series of tasks, each of which has a scaled configuration with certain volume, and all necessary input files necessary for starting a VAPS or LAMMPS relaxation. run sends all the computational tasks to remote computational resources defined in a machine configuration file like machine.json, and automatically collect the results when remote calculations finish. post calculates the desired property from the collected results.

Relaxation
The relaxation of a structure should be carried out before calculating all other properties:

dpgen autotest make equi.json 
dpgen autotest run relax.json machine.json
dpgen autotest post equi.json 
If, for some reason, the main program terminated at stage run, one can easily restart with the same command. relax.json is the parameter file. An example for deepmd relaxation is given as:

{
	"structures":	"confs/mp-*",
	"interaction": {
		"type":		"deepmd",
		"model":	"frozen_model.pb",
                "type_map":     {"Al": 0, "Mg": 1}
	},
	"relaxation": {
	}
}
where the key structures provides the structures to relax. interaction is provided with deepmd, and other options are vasp, eam, meam...

Yuzhi:

We should notice that the interaction here should always be considered as a unified abstract class, which means that we should avoid repeating identifing which interaction we're using in the main code.
The structures here should always considered as a list, and the wildcard should be supported by using glob. Before all calculations , there is a stage where we generate the configurations.
The outputs of the relaxation are stored in the mp-*/00.relaxation directory.

ls mp-*
mp-1/relaxation  mp-2/relaxation  mp-3/relaxation
Other properties
Other properties can be computed in parallel:

dpgen autotest make properties.json 
dpgen autotest run properties.json machine.json
dpgen autotest post properties.json 
where an example of properties.json is given by

{
	"structures":	"confs/mp-*",
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
			"vol_step":	0.5
		},
		{
                        "type":         "elastic",
			"norm_deform":	2e-2,
			"shear_deform": 5e-2
		}
        ]
}
The dpgen packed all eos and elastic task and sends them to corresponding computational resources defined in machine.json. The outputs of a property, taking eos for example, are stored in

ls mp-*/ | grep eos
mp-1/eos_00  mp-2/eos_00  mp-3/eos_00
where 00 are suffix of the task.

Refine the calculation of a property
Some times we want to refine the calculation of a property from previous results. For example, when higher convergence criteria EDIFF and EDIFFG are necessary, and the new VASP calculation is desired to start from the previous output configration, rather than starting from scratch.

dpgen autotest make refine.json 
dpgen autotest run refine.json machine.json
with refine.json

{
	"properties": {
		"eos" : {
			"init_from_suffix":	"00",
                        "output_suffix":        "01",
			"vol_start":	10,
			"vol_end":	30,
			"vol_step":	0.5
		}
	}	
}
Configuration filter
Some times the configurations automatically generated are problematic. For example, the distance between the interstitial atom and the lattic is too small, then these configurations should be filtered out. One can set filters of configurations by

{
	"properties": {
		"intersitital" : {
			"supercell":	[3,3,3],
			"insert_atom":	["Al"],
			"conf_filters": [
				{  "min_dist": 2 }
			] 
		}
	}	
}
Footer


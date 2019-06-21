<span style="font-size:larger;">dpgen Manual</span>
========

# Table of contents
- [About dpgen](#About-dpgen)
- [genenrator](#generator)
- [auto_test](#auto_test)
  - [How to write `param.json`](##How-to-write-`param.json`)
  - [How to write `machine.json`](##How-to-write-`machine.json`)
- [Troubleshooting](#Troubleshooting)
# About dpgen
The deep potential generator
# generator
# auto_test
At this step, we assume that you have prepared some graph files like `graph.*.pb` and the particular pseudopotential `POTCAR`.
Then you only need one command to achieve automatic testing of physical properties.
```
python run.py param.json machine.json
```
## How to write `param.json`
We take Cu as an example to show the parameter settings of `param.json`.

The first part is the fundamental setting for particular alloy system. 
```
    "_comment": "models",
    "potcar_map" : {
	"Cu" : "/somewhere/POTCAR",
        "Zr" : "/elsewhere/POTCAR"
    },
    "conf_dir":"the folder of configuration",
    "task_type":"deepmd",
    "task":"eos",
```
You need to add the specified paths of necessary `POTCAR` files in "potcar_map". The different `POTCAR` paths are separated by commas.
Then you also need to add the folder path of particular configuration, which contains `POSCAR` file. For your convenience, we recommend that you use `gen_confs.py` to generate configurations by the following command.
```
python gen_confs.py Cu
```
It will store the various configurations of the given element or alloy in **confs** folder.

So far, there are  3 test types for "task_type" (i.e. "vasp", "deepmd" and "meam") and 7 test items for "task" as follows.

>00.**equi**:(default task) the equilibrium state, return energy and volume per atom

>01.**eos**: the equation of state

>02.**elastic**: the elasticity like Young's module

>03.**vacancy**: the vacancy defect 

>04.**interstitial**: the interstitial defect

>05.**surf**: the surface energy

>06.**phonon**:(beta version) the phonon specturm

If you want to test all items, you can choose "all".

The second part is the computational settings for vasp and lammps. The most important setting is to add the folder path of deepmd model and supply the corresponding element type map.
```
"vasp_params":	{
	"ecut":		520,
	"ediff":	1e-6,
	"kspacing":	0.32,
	"kgamma":	false,
	"npar":		1,
	"kpar":		1,
	"_comment":	" that's all "
    },
    "deepmd_model_dir":	"the folder of deepmd model",
    "deepmd_type_map":	[
	"Cu"
    ],
    "meam_potfile_dir":	"meam",
    "meam_type_map":	[
	"Al", "Si", "Mg", "Cu", "Fe"
    ],
    "meam_potfile":	[
	"library.meam",
	"AlSiMgCuFe.meam"
    ],
    "meam_param_type":	[
	"AlS", "SiS", "MgS", "CuS", "FeS"
    ],
```
The last part is the optional settings for various tasks mentioned above. You can change the parameters according to actual needs.
```
    "_comment":"00.equi",
    "store_stable":true,

    "_comment": "01.eos",
    "vol_start":	6,
    "vol_end":		16,
    "vol_step":		0.5,
    "store_fix":false,

    "_comment": "02.elastic",
    "norm_deform":	2e-2,
    "shear_deform":	5e-2,
    
    "_comment":"03.vacancy",
    "supercell":[2,2,2],

    "_comment":"04.interstitial",
    "insert_ele":["Cu"],
    "reprod-opt":false,

    "_comment": "05.surface",
    "min_slab_size":	10,
    "min_vacuum_size":	11,
    "_comment": "pert xz to work around vasp bug...",
    "pert_xz":		0.01,
    "max_miller": 2,
    "static-opt":false,
    "store_relax":false,    

    "_comment":"06.phonon",
    "supercell_matrix":[2,2,2],
    "band":"0 1 0  0.5 1 0.5  0.375 0.75 0.375  0  0  0  0.5 0.5 0.5",

    "_comment":	"that's all"
```

## How to write `machine.json`
(Same as generator)
You need to set up your working environment through this file, which is divided into three parts, **deepmd**, **lammps** and **vasp**.
Taking **deepmd** as an example, you need to add the folder path of **deepmd** and supply the settings of train_machine and train_resource. The settings of **lammps** and **vasp** are similar to it.
```
    "deepmd_path":	"the folder of deepmd",
    "train_machine":	{
	"machine_type":	"slurm",
	"hostname" :	"localhost",
	"port" :	22,
	"username":	"username",
    	"password":     "password",
	"work_path" :	"the path of workplace",
	"_comment" :	"that's all"
    },
    "train_resources":	{
	"numb_node":	1,
	"numb_gpu":	1,
	"task_per_node":7,
	"source_list":	[ "the path of deepmd source" ],
	"module_list":	[ ],
	"time_limit":	"23:0:0",
	"mem_limit":	32,
	"_comment":	"that's all"
    },
```
# Troubleshooting

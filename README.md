# Workflow for DP-GEN


## Table of Contents

- [Workflow for DP-GEN](#workflow-for-dp-gen)
  * [Table of Contents](#table-of-contents)
  * [Basic structure of DPGEN](#basic-structure-of-dpgen)
  * [User flows](#user-flows)
  * [Preparing Data](#preparing-data)
    + [param.json](#paramjson)
  * [Main Process of Generator](#main-process-of-generator)
    + [Basics](#basics)
    + [machine.json](#machinejson)
    + [param.json](#paramjson-1)
    + [Test case : Methane](#test-case---methane)
  * [Doing Auto_test](#doing-auto-test)
    + [param.json](#paramjson-2)
  * [Appendix and FAQ](#appendix-and-faq)


## Installation
One can download the source code of dpgen by 
```bash
git clone https://github.com/deepmodeling/dpgen.git 
```
then use `setup.py` to install the module
```bash
cd dpgen
python setup.py install --user
```
With this command, the dpgen executable is install to `$HOME/.local/bin/dpgen`. You may want to export the `PATH` by
```bash
export PATH=$HOME/.local/bin/dpgen:$PATH
```
To test if the installation is successful, you may execute
```bash
dpgen -h
```
and if everything works, it gives
```
DeepModeling
------------

Version: 0.1.0
Path:    /home/wanghan/.local/lib/python3.6/site-packages/dpgen-0.1.0-py3.6.egg/dpgen
Date:    Jun 27, 2019

usage: dpgen [-h] {init,run,test} ...

dpgen is a convenient script that uses DeepGenerator to prepare initial data,
drive DeepMDkit and analyze results. This script works based on several sub-
commands with their own options. To see the options for the sub-commands, type
"dpgen sub-command -h".

positional arguments:
  {init,run,test}
    init           dpgen initial data preparation tools.
    run            Runing DeepMD with generator model.
    test           auto test for deep potential.

optional arguments:
  -h, --help       show this help message and exit

Author: DeepGenTeam Version: 0.1.0 Last updated: 2019.06.26
```



## Basic structure of DPGEN

One can easily download the code of DPGEN and go into the main folder:

Then you will see three sub-folders, namely
```
data		generator	auto_test
```

+ Folder `data`: As the first step, you may prepare the initital data here.

+ Folder `generator`: The main exploration process will take place here.
+ Folder `auto_test` : You can use this module to undertake relevant tests, based on the  model trained in previous process and Pymatgen.


It's necessary to mention that, the various and specific configurations for the implementation of DPGEN is realised by modifying some jsons file. After that, all work can be automatically done by DPGEN! We wil give an explicit illustration for the parameters needed in the following part.

## Preparing Data
The main code of this part is `gen.py`.
You may get the instructions of the codes, by:
```
python gen.py -h
```
We prepare some successful examples in the sub-folder `jsons`, next we take `jsons/al.fcc.111.json` as an example.

1. `python gen.py jsons/al.fcc.111.json 1`
you will create a foder al.fcc.01x01x01/00.place_ele/sys-0004, go to it and that will be the file for VASP relaxation;
you may do vasp calculation in that folder, the obtained CONTCAR in that folder is the only necessary one for the next step.
2. `python gen.py jsons/al.fcc.111.json 2` the code will collect relaxed configurations, perturb the system with 0.03 (3%) for box and 0.01 angstrom for atomic position (may be modified in the json file);
3. `python gen.py jsons/al.fcc.111.json 3` for perturbed systems, set up vasp jobs and run a 10-step short MD (you may run it in the corresponding al..fcc.01x01x01/02.md/sys-0004 folder);
and
4. `python gen.py jsons/al.fcc.111.json 4` collect vasp md confs, make the data prepared for the DeePMD training.

For Al, we prepared data for fcc.111, hcp.111, sc.111, and diamond.111 to kick off the initial training.
For later DeePMD exploration, we prepared fcc.222, hcp.222, sc.222, and diamond.222, these only serve as initial structures for DeePMD but did not undergo a VASP calculation. Starting from this step, only those picked up from the exploration step are sent to VASP (essentially around 0.003% of the confs explored by DeePMD are sent to VASP).
### param.json 

To be finished.
## Main Process of Generator

### Basics
The main code of this part is `run.py`

The command is simply `python run.py PARAM MACHINE`, where `PARAM` and `MACHINE` are both json files.

The whole process will contain a series of iterations, sorted by some rules such as heating the system to certain temperature.

In each iteration, there are three stages of work, namely, `00.train  01.model_devi  02.fp`.

+ 00.train: DP-GEN will train several (default 4) models based on initial and generated data. The only difference between these models is the random seed for neural network initialization.

+ 01.model_devi : represent for model-deviation. DP-GEN will use models obtained from 00.train to run Molecular Dynamics(default LAMMPS). Larger deviation for structure properties (default is force of atoms) means less accuracy of the models. Using this criterion, a few fructures will be selected and put into next stage `02.fp` for more accurate calculation based on First Principles.

+ 02.fp : Selected structures will be calculated by first principles methods(default VASP). DP-GEN will obtain some new data and put them together with initial data and data generated in previous iterations. After that a new training will be set up and DP-GEN will enter next iteration!

DP-GEN identifies the current stage by a record file, `record.dpgen`, which will be created and upgraded by codes.Each line contains two number: the first is index of iteration, and the second ,ranging from 0 to 9 ,records which stage in each iteration is currently running.

0,1,2 correspond to make_train, run_train, post_train. DP-GEN will write scripts in `make_train`, run the task by specific machine in `run_train` and collect result in `post_train`. The record for model_devi and fp stage follows similar rules.

The setting of json files is the most important and time-cosuming issue.If everything is set up correctly, the whole iteration is supposed to automatically done, without any extra work!



### machine.json

When switching into a new machine, you may modifying the `machine.json`, according to the actual circumstance. Once you have finished, the `machine.json` can be re-used for any DP-GEN tasks without any extra efforts.

An example for `machine.json` is:
```
{
    "deepmd_path":	"/sharedext4/local/deepmd-kit-0.12.4/",
    "train_machine":	{
	"machine_type":	"slurm",
	"hostname" :	"localhost",
	"port" :	22,
	"username":	"root",
	"work_path" :	"/sharedext4/generator/example/deep.gen/generator/ch4/",
	"_comment" :	"that's all"
    },
    "train_resources":	{
	"numb_node":	1,
	"numb_gpu":	1,
	"task_per_node":8,
	"partition" : "GPU-H",
	"exclude_list" : [],
	"source_list":	[ "/sharedext4/local/deepmd-kit-0.12.4/bin/activate" ],
	"module_list":	[ ],
	"time_limit":	"23:0:0",
	"mem_limit":	32,
	"_comment":	"that's all"
    },

    "lmp_command":	"/sharedext4/softwares/lammps/bin/lmp_serial",
    "model_devi_group_size":	1,
    "_comment":		"model_devi on localhost",
    "model_devi_machine":	{
	"machine_type":	"slurm",
	"hostname" :	"localhost",
	"port" :	22,
	"username":	"root",
	"work_path" :	"/sharedext4/generator/example/deep.gen/generator/ch4/",
	"_comment" :	"that's all"
    },
    "_comment": " if use GPU, numb_nodes(nn) should always be 1 ",
    "_comment": " if numb_nodes(nn) = 1 multi-threading rather than mpi is assumed",
    "model_devi_resources":	{
	"numb_node":	1,
	"numb_gpu":	0,
	"task_per_node":1,
	"source_list":	["/sharedext4/local/deepmd-kit-0.12.4/bin/lammps.activate" ],
	"module_list":	[ ],
	"time_limit":	"19:0:0",
	"mem_limit":	32,
	"partition" : "GPU-All",
	"_comment":	"that's all"
    },

    "_comment":         "fp on localhost ",
    "fp_command":       "/sharedext4/vasp/vasp.5.4.4/bin/vasp_std",
    "fp_group_size":    1,
    "fp_machine":       {
        "machine_type": "slurm",
        "hostname" :    "localhost",
        "port" :        22,
        "username":     "root",
        "work_path" :   "/sharedext4/generator/example/deep.gen/generator/ch4/",
        "_comment" :    "that's all"
    },
    "fp_resources":     {
        "numb_node":    1,
        "task_per_node":1,
        "numb_gpu":     0,
	"exclude_list" : [],
        "source_list":  ["/sharedext4/softwares/source/vasp_cpu.activate" ],
        "module_list":  [],
	"with_mpi" : 1,
	"time_limit":   "0:10:0",
	"partition" : "GPU-All",
        "_comment":     "that's all"
    },


    "_comment":		" that's all "
}
```
For convenience, we use TASK to represent train or model_devi or fp, for similar parameters in the json file.

+ `deepmd_path` (string), is the installed directory of DeepMD-Kit, which should contain `bin lib include`.

+ `TASK_machine` (dict), modifies the settings of the machine for TASK. 

+ `TASK_resources` (dict), modify the resources needed for calculation. In the dict, `numb_node`(integer) is node count required for the job.`numb_gpu`(integer) specializes whether the job apply for the GPU resources. `task_per_node` (integer) is number of tasks to be launched. `source_list` (list) contains the environment needed for certain job. For example, if "env" is in the list, 'source env' will be written in the script for slurm. `module_list` (list) If "moduleA" is in the list, "module load moduleA" will be written.`time_limit`(string) is the maximum time permitted for the job.`mem_limit`(string) is the maximum memory permitted to apply for.


+ `TASK_command` (string), is the command for call TASK.

+ `TASK_group_size` (integer) DP-GEN will put these jobs together in one submitting script.



### param.json 

In param.json, you can specialize the task as you expect.


```
{
    "type_map":		["H","C"],
    "mass_map":		[1, 12],

    "init_data_prefix":	"/sharedext4/generator/example/deep.gen/data/",

    "init_data_sys":	[
	"/sharedext4/generator/example/deep.gen/data/CH4.POSCAR.01x01x01/02.md/sys-0004-0001/deepmd"
			],
    "init_batch_size":	[
	8
			],
    "sys_configs":	[
    ["/sharedext4/generator/example/deep.gen/data/CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00000[0-4]/POSCAR"],
    ["/sharedext4/generator/example/deep.gen/data/CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00000[5-9]/POSCAR"],
    ["/sharedext4/generator/example/deep.gen/data/CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00001*/POSCAR"],
    ["/sharedext4/generator/example/deep.gen/data/CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00002*/POSCAR"]
	],
    "_comment":		"0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25",
    "sys_batch_size":	[
	8, 8, 8, 8
    ],


    "_comment":		" 00.train ",
    "numb_models":	4,
    "train_param":	"input.json",
    "default_training_param" : {
	"_comment": " model parameters",
	"use_smooth":		true,
	"sel_a":		[16,4],
	"rcut_smth":		0.50,
	"rcut":			5,
	"filter_neuron":	[10, 20, 40],
	"filter_resnet_dt":	false,
	"n_axis_neuron":	12,
	"n_neuron":		[120,120,120],
	"resnet_dt":		true,
	"coord_norm":		true,
	"type_fitting_net":	false,

	"_comment": " traing controls",
	"systems":		[],
	"set_prefix":		"set",
	"stop_batch":		20000,
	"batch_size":		1,
	"start_lr":		0.001,
	"decay_steps":		100,
	"decay_rate":		0.95,
	"seed":			0,

	"start_pref_e":		0.02,
	"limit_pref_e":		2,
	"start_pref_f":		1000,
	"limit_pref_f":		1,
	"start_pref_v":		0.0,
	"limit_pref_v":		0.0,

	"_comment": " display and restart",
	"_comment": " frequencies counted in batch",
	"disp_file":		"lcurve.out",
	"disp_freq":		1000,
	"numb_test":		4,
	"save_freq":		1000,
	"save_ckpt":		"model.ckpt",
	"load_ckpt":		"model.ckpt",
	"disp_training":	true,
	"time_training":	true,
	"profiling":		false,
	"profiling_file":	"timeline.json",

	"_comment":		"that's all"
    },

    "_comment":		" 01.model_devi ",
    "_comment": "model_devi_skip: the first x of the recorded frames",
    "model_devi_dt":		0.002,
    "model_devi_skip":		0,
    "model_devi_f_trust_lo":	0.05,
    "model_devi_f_trust_hi":	0.15,
    "model_devi_e_trust_lo":	1e10,
    "model_devi_e_trust_hi":	1e10,
    "model_devi_clean_traj":	true,
    "model_devi_jobs":	[
	{"sys_idx": [0],
	"temps": [  100], "press": [1.0], "trj_freq": 10, "nsteps": 300,  "ensemble": "npt", "_idx": "00"},
	{"sys_idx": [1],
	"temps": [  100], "press": [1.0], "trj_freq": 10, "nsteps": 1000,  "ensemble": "npt", "_idx": "01"},
	{"sys_idx": [2],
	"temps": [  100], "press": [1.0], "trj_freq": 10, "nsteps": 3000,  "ensemble": "npt", "_idx": "02"},
	{"sys_idx": [3],
	"temps": [  100], "press": [1.0], "trj_freq": 10, "nsteps": 3000,  "ensemble": "npt", "_idx": "03"}
	],


    "_comment":		" 02.fp ",
    "fp_style":		"vasp",
    "shuffle_poscar":	false,
    "fp_task_max":	30,
    "fp_task_min":	5,
    "fp_pp_path":	"/sharedext4/generator/example/deep.gen/data/ch4/",
    "fp_pp_files":	["POTCAR"],
    "fp_params":	{
	"_comment": "given in unit depending on the fp method",
	"ecut":		400,
	"ediff":	1e-6,
	"kspacing":	2,
	"_comment": "gauss, mp:N(methfessel-paxton:order by default order=1), fd(Fermi-Dirac)",
	"smearing":	"gauss",
	"sigma":	0.05,
	"_comment": "only for vasp, can be NONE, SCAN, TPSS, RTPSS, M06L or MBJ",
	"metagga":	"NONE",
	"npar":		4,
	"kpar":		1,
	"_comment":	" that's all "
    },
    "_comment":		" that's all "
}
```

+ `type_map` (list of string ), is the atom type you want to explore.

+ `mass_map` (list of float), contains the corresponding standard atom weights.

+ `init_data_prefix` (string) is the location of folder containg initial data.

+ `init_data_sys` (list of string) is the exact directory of initial data.

+ `init_batch_size` (list of integer) Each number in the list is the batch_size for training of corresponding system in init_data_sys. Therefore, the size of these two lists should be the same. One recommend rule for setting is that init_batch_size mutiply number of atoms in one structure should be larger than 32.

+ `sys_configs` (list of list) contains the directory of structure to be explored in iterations.Default file format is POSCAR.

+ `sys_batch_size` (list of integer) Each number in the list is the batch_size for training of corresponding system in sys_configs. Setting ules are similar to init_batch_size

+ `num_models` : (integer) is the number of models to be trained in 00.train

+ `default_training_param` (dict): is training parameters for DeepMD in 00.train. You can find instructions from here:  https://github.com/deepmodeling/deepmd-kit

+ `model_devi_f_trust_lo` and `model_devi_f_trust_hi` (float) are the lower bound and upper bound of force for the selection in 01.model_devi. One recommended setting for the lower bound is twice of the traning error.

+ `model_devi_clean_traj`(boolin), decides whether to clean traj folders in Moldecular-Dynamics since they are too large.

+ `model_devi_jobs` (list of dict), contains the settings for model deviation. Each dict in the list correspond one iteration.In the dict, `sys_idx`(list), choose which system to be tested, the index here corresponds exactly to the 'sys_configs'. `temp` and `press` (list) are the temperature and pressure in MD. `trj_freq` is the frequecy of trajectory saved in MD. `nsteps` is the running steps of MD. `ensembles`(string), determines which ensemble used in MD, options include "npt" and "nvt".

+ `fp_style`(string): chooses the software for First Principles.

+ `fp_task_max` and `fp_task_min` are the maximum and minimum numbers of structures to calculate in 02.fp.

+ `fp_pp_path` and `fp_pp_files` determine the location psuedo-potential file to be used for 02.fp.

+ `fp_params` are parameters for 02.fp.


### Test case : Methane
We have provided an easy case of modelling methane (CH4) in $100K$ for you in folder generator. In the following part, we will take it as an example to give an instuction of DP-GEN.

We have prepared data for CH4 in $50K$ using the tools in folder `data` and put them in the folder `data/CH4.POSCAR.01.01.01/` which contains 300 frames of initial data. 

Next we'll try explore CH4 in $100K$ with DP-GEN by calculating faily a few new frames of data.

The two jsons mentioned above have been placed in folder `generator/ch4/`.You need to modify `machine.json` and such following keys in `param.json` : `init_data_prefix`, `init_data_sys`, `sys_configs` and `fp_pp_path` to run the example successfully on your own machine .

(Noted by Yuzhi: We'd better use relative path here for convenience.)

Once the jsons have been set correctly, you may simply run DP-GEN by `python run.py ch4/param.json ch4/machine.json`! 


## Doing Auto_test
At this step, we assume that you have prepared some graph files like `graph.*.pb` and the particular pseudopotential `POTCAR`.

The main code of this step is 
```
dpgen test PARAM MACHINE
```
where `PARAM` and `MACHINE` are both json files. `MACHINE` is the same as above.

The whole program contains a series of tasks shown as follows. In each task, there are three stages of work, generate, run and compute.
+ `00.equi`:(default task) the equilibrium state

+ `01.eos`: the equation of state

+ `02.elastic`: the elasticity like Young's module

+ `03.vacancy`: the vacancy formation energy

+ `04.interstitial`: the interstitial formation energy

+ `05.surf`: the surface formation energy

### param.json 

We take Al as an example to show the parameter settings of `param.json`.
The first part is the fundamental setting for particular alloy system. 
```
    "_comment": "models",
    "potcar_map" : {
	"Al" : "/somewhere/POTCAR"
    },
    "conf_dir":"confs/Al/std-fcc",
    "key_id":"API key of Material project",
    "task_type":"deepmd",
    "task":"eos",
```
You need to add the specified paths of necessary `POTCAR` files in "potcar_map". The different `POTCAR` paths are separated by commas.
Then you also need to add the folder path of particular configuration, which contains `POSCAR` file. 
```
"confs/[element or alloy]/[std-* or mp-**]"
std-*: standard structures, * can be fcc, bcc, hcp and so on.
mp-**: ** means Material id from Material Project.
```
Usually, if you add the relative path of POSCAR as the above format,
`dpgen test` will check the existence of such file and automatically downloads the standard and existed configurations of the given element or alloy from Materials Project and stores them in **confs** folder, which needs the API key of Materials project.

+ `task_type` contains 3 optional types for testing, i.e. **vasp**, **deepmd** and **meam**.
+ `task` contains 7 options, **equi**, **eos**, **elastic**, **vacancy**, **interstitial**, **surf** and **all**. The option **all** can do all the tasks. 

It is worth noting that the subsequent tasks need to rely on the calculation results of the equilibrium state, so it is necessary to give priority to the calculation of the equilibrium state while testing. And due to the stable consideration, we recommand you to test the equilibrium state of **vasp** before other tests.

The second part is the computational settings for vasp and lammps. The most important setting is to add the folder path `model_dir` of **deepmd** model and supply the corresponding element type map. Besides, `dpgen test` also is able to call common lammps packages, such as **meam**. 
```json
"vasp_params":	{
	"ecut":		650,
	"ediff":	1e-6,
	"kspacing":	0.1,
	"kgamma":	false,
	"npar":		1,
	"kpar":		1,
	"_comment":	" that's all "
    },
    "lammps_params":    {
        "model_dir":"somewhere/example/Al_model",
        "type_map":["Al"],
        "model_name":false,
        "model_param_type":false
    },
```
The last part is the optional settings for various tasks mentioned above. You can change the parameters according to actual needs.
```
    "_comment":"00.equi",
    "store_stable":true,
```
+ `store_stable`:(boolean) whether to store the stable energy and volume

```
    "_comment": "01.eos",
    "vol_start":	12,
    "vol_end":		22,
    "vol_step":		0.5,
```
+ `vol_start`, `vol_end` and `vol_step` determine the volumetric range and accuracy of the **eos**.

```
    "_comment": "02.elastic",
    "norm_deform":	2e-2,
    "shear_deform":	5e-2,
```
+ `norm_deform` and `shear_deform` are the scales of material deformation. 
This task uses the stress-strain relationship to calculate the elastic constant.

```
    "_comment":"03.vacancy",
    "supercell":[3,3,3],
```
+ `supercell`:(list of integer) the supercell size used to generate vacancy defect and interstitial defect
```
    "_comment":"04.interstitial",
    "insert_ele":["Al"],
    "reprod-opt":false,
```
+ `insert_ele`:(list of string) the elements used to generate point interstitial defect
+ `repord-opt`:(boolean) whether to reproduce trajectories of interstitial defect

```
    "_comment": "05.surface",
    "min_slab_size":	10,
    "min_vacuum_size":	11,
    "_comment": "pert xz to work around vasp bug...",
    "pert_xz":		0.01,
    "max_miller": 2,
    "static-opt":false,
    "relax_box":false,    
```
+ `min_slab_size` and `min_vacuum_size` are the minimum size of slab thickness  and  the vacuume width.
+ `pert_xz` is the perturbation through xz direction used to compute surface energy.
+ `max_miller` (integer) is the maximum miller index
+ `static-opt`:(boolean) whether to use atomic relaxation to compute surface energy. if false, the structure will be relaxed.
+ `relax_box`:(boolean) set true if the box is relaxed, otherwise only relax atom positions.



## Appendix and FAQ

# DP-GEN Manual

## Table of Contents

 * [DP-GEN Manual](#dp-gen-manual)
      * [Table of Contents](#table-of-contents)
      * [About DP-GEN](#about-dp-gen)
         * [Highlighted features](#highlighted-features)
         * [Code structure and interface](#code-structure-and-interface)
         * [License and credits](#license-and-credits)
      * [Download and Install](#download-and-install)
      * [Init: Preparing Initial Data](#init-preparing-initial-data)
         * [Init_bulk](#init_bulk)
         * [Init_surf](#init_surf)
      * [Run: Main Process of Generator](#run-main-process-of-generator)
      * [Test: Evaluating performances of model](#test-evaluating-performances-of-model)
      * [Set up machine](#set-up-machine)
      * [Troubleshooting](#troubleshooting)
      * [License](#license)


## About DP-GEN
DP-GEN (Deep Generator)  is a software written in Python, delicately designed to generate a deep learning based model of interatomic potential energy and force field. DP-GEN is depedent on DeepMD-kit (https://github.com/deepmodeling/deepmd-kit/blob/master/README.md). With highly scalable interface with common softwares for molecular simulation, DP-GEN is capable to  automatically prepare scripts and maintain job queues on HPCs (High Performance Cluster) and analyze results
### Highlighted features
+ **Accurate and efficient**: DP-GEN is capable to sample tens of million structures and select only a few (< 0.01%) for first principles calculation. DP-GEN will finally obtain a uniformly accurate model.
+ **User-friendly and automatic**: Users may install and run DP-GEN easily. Once succusefully running, DP-GEN may dispatch and handle all jobs on HPCs, and thus there's no need for any personal effort. 
+ **Highly scalable**: With modularized code structures, users and developers can easily extend DP-GEN for their most relevant needs. DP-GEN currently supports for HPC systems (Slurm, PBS, LSF and cloud machines ), Deep Potential interface with DeePMD-kit, MD interface with LAMMPS  and *ab-initio* calculation interface with VASP, PWSCF, and Gaussian. We're sincerely welcome and embraced to users' contributions, with more possibilities and cases to use DP-GEN.

### Code structure and interface
+ dpgen:
    * data: source codes for preparing initial data of bulk and surf systems.

    * generator: source codes for main process of deep generator.

    * auto_test : source code for undertaking materials property analysis.
    * remote : source code for automatically submiting scripts,maintaining job queues and collecting results. 
    * database
+ examples : providing example JSON files.

+ tests : unittest tools for developers.

One can easily run DP-GEN with :
```
dpgen TASK PARAM MACHINE
```

where TASK is the key word, PARAM and MACHINE are both JSON files.

Options for TASK:
* `init_bulk` : Generating initial data for bulk systems. 
* `init_surf` : Generating initial data for bulk systems.
* `run` : Main process of Generator.
* `test`: Evaluating the performance of models.

### License and credits

## Download and Install 
One can download the source code of dpgen by 
```bash
git clone https://github.com/deepmodeling/dpgen.git 
```
then use `setup.py` to install the module
```bash
cd dpgen
pip install --user .
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

usage: dpgen [-h] {init_surf,init_bulk,run,test,db} ...

dpgen is a convenient script that uses DeepGenerator to prepare initial data,
drive DeepMDkit and analyze results. This script works based on several sub-
commands with their own options. To see the options for the sub-commands, type
"dpgen sub-command -h".

positional arguments:
  {init_surf,init_bulk,run,test,db}
    init_surf           dpgen initial data preparation tools for surface
                        systems.
    init_bulk           dpgen initial data preparation tools for bulk systems.
    run                 Runing DeepMD with generator model.
    test                auto test for deep potential.
    db                  data collection from dpgen.

optional arguments:
  -h, --help            show this help message and exit
```
+ Option `init`: As the first step, you may prepare the initital data here.

+ Option `run`: The main exploration process will take place here.

+ Option `test` : You can use this module to undertake relevant tests, based on the  model trained in previous process and Pymatgen.

+ Option `db` :  a convinent tool to collect data from dpgen iterations and save it into .json file, which can be directly saved into 
                 `MongoDB` or `ElasticSearch` database





## Init: Preparing Initial Data

### Init_bulk

You may prepare initial data for bulk systems with VASP by:

```bash
dpgen init_bulk PARAM MACHINE 
```

Basically `init_bulk` can be devided into four parts , denoted as `stages` in `PARAM`:
1. Relax in folder `00.place_ele`
2. Pertub and scale in folder `01.scale_pert`
3. Run a shor AIMD in folder `02.md`
4. Collect data in folder `02.md`.
 
All stages must be **in order**. One doesn't need to run all stages. For example, you may run stage 1 and 2, generating supercells as starting point of exploration in `dpgen run`.

Following is an example for `PARAM`, which generates data from a typical structure hcp.
```json
{
    "stages" : [1,2,3,4],
    "cell_type":    "hcp",
    "latt":     4.479,
    "super_cell":   [2, 2, 2],
    "elements":     ["Mg"],
    "potcars":      ["....../POTCAR"],
    "relax_incar": "....../INCAR_metal_rlx",
    "md_incar" : "....../INCAR_metal_md",
    "scale":        [1.00],
    "skip_relax":   false,
    "pert_numb":    2,
    "md_nstep" : 5,
    "pert_box":     0.03,
    "pert_atom":    0.01,
    "coll_ndata":   5000,
    "_comment":     "that's all"
}
```

If you want to specify a structure as starting point for `init_bulk`, you may set in `PARAM` as follows.

```json
"from_poscar":	true,
"from_poscar_path":	"....../C_mp-47_conventional.POSCAR",
```
The following table gives explicit descriptions on keys in `PARAM`. 

The bold notation of key (such aas **type_map**) means that it's a necessary key. 

 Key  | Type          | Example                                                      | Discription                                                      |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
| **stages** | List of Integer | [1,2,3,4] | Stages for `init_bulk` 
| **Elements** | List of String | ["Mg"] | Atom types
|  cell_type | String  | "hcp" | Specifying which typical structure to be generated. **Options** include fcc, hcp, bcc, sc, diamond.
| latt | Float | 4.479 | Lattice constant for single cell.
| from_poscar | Boolean | True | Deciding whether to use a given poscar as the beginning of relaxation. If it's true, keys (`cell_type`, `latt`) will be aborted. Otherwise, these two keys are **necessary**.
| from_poscar_path | String | "....../C_mp-47_conventional.POSCAR" | Path of POSCAR
| relax_incar | String | "....../INCAR" | Path of INCAR for relaxation in VASP. **Necessary** if `stages` include 1.
| md_incar | String |  "....../INCAR" | Path of INCAR for MD in VASP. **Necessary** if `stages` include 3.| 
| **scale** | List of float | [0.980, 1.000, 1.020] | Scales for transforming cells.
| **skip_relax** | Boolean | False | If it's true, you may directly run stage 2 (pertub and scale) using an unrelaxed POSCAR.
| **pert_numb** | Integer | 30 | Number of pertubations for each POSCAR. 
| **md_nstep** | Integer | 10 | Steps of AIMD in stage 3. If it's not equal to settings via `NSW` in `md_incar`, DP-GEN will follow `NSW`.
| **pert_box** | Float | 0.03 | 
| **pert_atom** | Float | 0.01 | 
| **coll_ndata** | Integer | 5000 | 
### Init_surf




## Run: Main Process of Generator


You may call the main process by:
`dpgen run PARAM MACHINE`.


The whole process of generator will contain a series of iterations, succussively undertaken in order such as heating the system to certain temperature.

In each iteration, there are three stages of work, namely, `00.train  01.model_devi  02.fp`.

+ 00.train: DP-GEN will train several (default 4) models based on initial and generated data. The only difference between these models is the random seed for neural network initialization.

+ 01.model_devi : represent for model-deviation. DP-GEN will use models obtained from 00.train to run Molecular Dynamics(default LAMMPS). Larger deviation for structure properties (default is force of atoms) means less accuracy of the models. Using this criterion, a few fructures will be selected and put into next stage `02.fp` for more accurate calculation based on First Principles.

+ 02.fp : Selected structures will be calculated by first principles methods(default VASP). DP-GEN will obtain some new data and put them together with initial data and data generated in previous iterations. After that a new training will be set up and DP-GEN will enter next iteration!

DP-GEN identifies the current stage by a record file, `record.dpgen`, which will be created and upgraded by codes.Each line contains two number: the first is index of iteration, and the second ,ranging from 0 to 9 ,records which stage in each iteration is currently running.

0,1,2 correspond to make_train, run_train, post_train. DP-GEN will write scripts in `make_train`, run the task by specific machine in `run_train` and collect result in `post_train`. The records for model_devi and fp stage follow similar rules.


In `PARAM`, you can specialize the task as you expect.


```json
{
  "type_map": [
    "H",
    "C"
  ],
  "mass_map": [
    1,
    12
  ],
  "init_data_prefix": "....../init/",
  "init_data_sys": [
    "CH4.POSCAR.01x01x01/02.md/sys-0004-0001/deepmd"
  ],
  "init_batch_size": [
    8
  ],
  "sys_configs_prefix": "....../init/",
  "sys_configs": [
    [
      "CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00000*/POSCAR"
    ],
    [
      "CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00001*/POSCAR"
    ]
  ],
  "sys_batch_size": [
    8,
    8,
    8,
    8
  ],
  "_comment": " that's all ",
  "numb_models": 4,
  "train_param": "input.json",
  "default_training_param": {
    "_comment": "that's all",
    "use_smooth": true,
    "sel_a": [
      16,
      4
    ],
    "rcut_smth": 0.5,
    "rcut": 5,
    "filter_neuron": [
      10,
      20,
      40
    ],
    "filter_resnet_dt": false,
    "n_axis_neuron": 12,
    "n_neuron": [
      100,
      100,
      100
    ],
    "resnet_dt": true,
    "coord_norm": true,
    "type_fitting_net": false,
    "systems": [],
    "set_prefix": "set",
    "stop_batch": 40000,
    "batch_size": 1,
    "start_lr": 0.001,
    "decay_steps": 200,
    "decay_rate": 0.95,
    "seed": 0,
    "start_pref_e": 0.02,
    "limit_pref_e": 2,
    "start_pref_f": 1000,
    "limit_pref_f": 1,
    "start_pref_v": 0.0,
    "limit_pref_v": 0.0,
    "disp_file": "lcurve.out",
    "disp_freq": 1000,
    "numb_test": 4,
    "save_freq": 1000,
    "save_ckpt": "model.ckpt",
    "load_ckpt": "model.ckpt",
    "disp_training": true,
    "time_training": true,
    "profiling": false,
    "profiling_file": "timeline.json"
  },
  "model_devi_dt": 0.002,
  "model_devi_skip": 0,
  "model_devi_f_trust_lo": 0.05,
  "model_devi_f_trust_hi": 0.15,
  "model_devi_clean_traj": true,
  "model_devi_jobs": [
    {
      "sys_idx": [
        0
      ],
      "temps": [
        100
      ],
      "press": [
        1.0
      ],
      "trj_freq": 10,
      "nsteps": 300,
      "ensemble": "nvt",
      "_idx": "00"
    },
    {
      "sys_idx": [
        1
      ],
      "temps": [
        100
      ],
      "press": [
        1.0
      ],
      "trj_freq": 10,
      "nsteps": 3000,
      "ensemble": "nvt",
      "_idx": "01"
    }
  ],
  "fp_style": "vasp",
  "shuffle_poscar": false,
  "fp_task_max": 20,
  "fp_task_min": 1,
  "fp_pp_path": "....../methane/",
  "fp_pp_files": [
    "POTCAR"
  ],
  "fp_incar": "....../INCAR_methane"
}
```

The following table gives explicit descriptions on keys in `PARAM`. 

The bold notation of key (such aas **type_map**) means that it's a necessary key. 

 Key  | Type          | Example                                                      | Discription                                                      |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
| *#Basics*
| **type_map** | List of string | ["H", "C"] | Atom types
| **mass_map** | List of float |  [1, 12] | Standard atom weights.
| *#Data*
 | init_data_prefix | String | "/sharedext4/.../data/" | Prefix of initial data directories 
 | ***init_data_sys*** | List of string|["CH4.POSCAR.01x01x01/.../deepmd"] |Directories of initial data. You may use either absolute or relative path here.
 | ***sys_format*** | String | "vasp/poscar" | Format of initial data. It will be `vasp/poscar` if not set.
 | **init_batch_size**   | String of integer     | [8]                                                            | Each number is the batch_size of corresponding system  for training in `init_data_sys`. One recommended rule for setting the `sys_batch_size` and `init_batch_size` is that `batch_size` mutiply number of atoms ot the stucture should be larger than 32. If set to `auto`, batch size will be 32 divided by number of atoms. |
  | sys_configs_prefix | String | "/sharedext4/.../data/" | Prefix of `sys_configs`
 | **sys_configs**   | List of list of string         | [<br />["/sharedext4/.../POSCAR"], <br />["....../POSCAR"]<br />] | Containing directories of structures to be explored in iterations.Wildcard characters are supported here. |
| ***sys_batch_size***      | List of integer   | [8, 8]                                                 | Each number  is the batch_size for training of corresponding system in `sys_configs`. If set to `auto`, batch size will be 32 divided by number of atoms. |
| *#Training*
| **numb_models**      | Integer      | 4 (recommend)                                                           | Number of models to be trained in `00.train`. |
| **default_training_param** | Dict | {<br />... <br />"use_smooth": true, <br/>"sel_a": [16, 4], <br/>"rcut_smth": 0.5, <br/>"rcut": 5, <br/>"filter_neuron": [10, 20, 40], <br/>...<br />} | Training parameters for `deepmd-kit` in `00.train`. <br /> You can find instructions from here: (https://github.com/deepmodeling/deepmd-kit)..<br /> We commonly let `stop_batch` = 200 * `decay_steps`. |
| *#Exploration*
| **model_devi_dt** | Float | 0.002 (recommend) | Timestep for MD | 
| **model_devi_skip** | Integer | 0 | Number of structures skipped for fp in each MD  
| **model_devi_f_trust_lo** | Float | 0.05 | Lower bound of forces for the selection.
 | **model_devi_f_trust_hi** | Float | 0.15 | Upper bound of forces for the selection
| **model_devi_e_trust_lo**  | Float | 1e10                                                         | Lower bound of energies for the selection. Recommend to set them a high number, since forces provide more precise information. Special cases such as energy minimization may need this. |
| **model_devi_e_trust_hi**  | Float | 1e10                                                         | Upper bound of energies for the selection. |
| **model_devi_clean_traj**  | Boolean | true                                                         | Deciding whether to clean traj folders in MD since they are too large. |
| **model_devi_jobs**        | [<br/>{<br/>"sys_idx": [0], <br/>"temps": <br/>[100],<br/>"press":<br/>[1],<br/>"trj_freq":<br/>10,<br/>"nsteps":<br/> 1000,<br/> "ensembles": <br/> "nvt" <br />},<br />...<br />] | List of dict | Settings for exploration in `01.model_devi`. Each dict in the list corresponds to one iteration. The index of `model_devi_jobs` exactly accord with index of iterations |
| **model_devi_jobs["sys_idx"]**    | List of integer           | [0]                                                          | Systems to be selected as the initial structure of MD and be explored. The index corresponds exactly to the `sys_configs`. |
| **model_devi_jobs["temps"]**  | List of integer | [50, 300] | Temperature (**K**) in MD 
| **model_devi_jobs["press"]**  | List of integer | [1,10] | Pressure (**Bar**) in MD 
| **model_devi_jobs["trj_freq"]**   | Integer  | 10                                                            | Frequecy of trajectory saved in MD.                   |
| **model_devi_jobs["nsteps"]**     | Integer      | 3000                                                         | Running steps of MD.                                  |
| **model_devi_jobs["ensembles"]** | String             | "nvt"                                    | Determining which ensemble used in MD, **options** include “npt” and “nvt”. |
| *#Labeling*
| ***use_clusters*** | Boolean | false | If set to `true`, clusters will be taken instead of the whole system. This option does not work with DeePMD-kit 0.x.
| ***cluster_cutoff***| Float | 3.5 | The cutoff radius of clusters if `use_clusters` is set to `true`.
| **fp_style** | string                | "vasp"                                                       | Software for First Principles. **Options** include “vasp”, “pwscf” and “gaussian” up to now. |
| **fp_task_max** | Integer            | 20                                                           | Maximum of  structures to be calculated in `02.fp` of each iteration. |
| **fp_task_min**     | Integer        | 5                                                            | Minimum of structures to calculate in `02.fp` of each iteration. |
| **fp_pp_path**   | String           | "/sharedext4/.../ch4/"                                       | Directory of psuedo-potential file to be used for 02.fp exists. |
| **fp_pp_files**    | List of string         | ["POTCAR"]                                                   | Psuedo-potential file to be used for 02.fp. Note that the order of elements should correspond to the order in `type_map`. |
|**fp_incar** | String | "/sharedext4/../ch4/INCAR" | Input file for VASP
|**keywords** | String or list | "mn15/6-31g** nosymm scf(maxcyc=512)" | Keywords for Gaussian input.
|**multiplicity**| Integer or String | 1 | Spin multiplicity for Gaussian input. If set to `auto`, the spin multiplicity will be detected automatically. If set to `frag`, the "fragment=N" method will be used.
|**nproc** | Integer| 4 | The number of processors for Gaussian input.

## Test: Evaluating performances of model
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


We take Al as an example to show the parameter settings of `param.json`.
The first part is the fundamental setting for particular alloy system. 
```json
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
```json
    "_comment":"00.equi",
    "store_stable":true,
```
+ `store_stable`:(boolean) whether to store the stable energy and volume

```json
    "_comment": "01.eos",
    "vol_start":	12,
    "vol_end":		22,
    "vol_step":		0.5,
```
+ `vol_start`, `vol_end` and `vol_step` determine the volumetric range and accuracy of the **eos**.

```json
    "_comment": "02.elastic",
    "norm_deform":	2e-2,
    "shear_deform":	5e-2,
```
+ `norm_deform` and `shear_deform` are the scales of material deformation. 
This task uses the stress-strain relationship to calculate the elastic constant.

```json
    "_comment":"03.vacancy",
    "supercell":[3,3,3],
```
+ `supercell`:(list of integer) the supercell size used to generate vacancy defect and interstitial defect
```json
    "_comment":"04.interstitial",
    "insert_ele":["Al"],
    "reprod-opt":false,
```
+ `insert_ele`:(list of string) the elements used to generate point interstitial defect
+ `repord-opt`:(boolean) whether to reproduce trajectories of interstitial defect

```json
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




## Set up machine
When switching into a new machine, you may modifying the `MACHINE`, according to the actual circumstance. Once you have finished, the `MACHINE` can be re-used for any DP-GEN tasks without any extra efforts.

An example for `MACHINE` is:
```json
{
  "train": [
    {
      "machine": {
        "machine_type": "slurm",
        "hostname": "localhost",
        "port": 22,
        "username": "1600017784",
        "work_path": "/gpfs/share/home/1600017784/yuzhi/methane/work"
      },
      "resources": {
        "numb_node": 1,
        "numb_gpu": 1,
        "task_per_node": 4,
        "partition": "AdminGPU",
        "exclude_list": [],
        "source_list": [
          "/gpfs/share/home/1600017784/env/train_tf112_float.env"
        ],
        "module_list": [],
        "time_limit": "23:0:0",
        "qos": "bigdata"
      },
      "deepmd_path": "/gpfs/share/software/deepmd-kit/0.12.4/gpu/gcc/4.9.0/tf1120-lowprec"
    }
  ],
  "model_devi": [
    {
      "machine": {
        "machine_type": "slurm",
        "hostname": "localhost",
        "port": 22,
        "username": "1600017784",
        "work_path": "/gpfs/share/home/1600017784/yuzhi/methane/work"
      },
      "resources": {
        "numb_node": 1,
        "numb_gpu": 1,
        "task_per_node": 2,
        "partition": "AdminGPU",
        "exclude_list": [],
        "source_list": [
          "/gpfs/share/home/1600017784/env/lmp_tf112_float.env"
        ],
        "module_list": [],
        "time_limit": "23:0:0",
        "qos": "bigdata"
      },
      "command": "lmp_serial",
      "group_size": 1
    }
  ],
  "fp": [
    {
      "machine": {
        "machine_type": "slurm",
        "hostname": "localhost",
        "port": 22,
        "username": "1600017784",
        "work_path": "/gpfs/share/home/1600017784/yuzhi/methane/work"
      },
      "resources": {
        "cvasp": true,
        "task_per_node": 4,
        "numb_gpu": 1,
        "exclude_list": [],
        "with_mpi": false,
        "source_list": [],
        "module_list": [
          "mpich/3.2.1-intel-2017.1",
          "vasp/5.4.4-intel-2017.1",
          "cuda/10.1"
        ],
        "time_limit": "120:0:0",
        "partition": "AdminGPU",
        "_comment": "that's All"
      },
      "command": "vasp_gpu",
      "group_size": 1
    }
  ]
}
```
Following table illustrates which key is needed for three types of machine: `train`,`model_devi`  and `fp`. Each of them is a list of dicts. Each dict can be considered as an independent environmnet for calculation. 

 Key   | `train`          | `model_devi`                                                    | `fp`                                                     |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
| machine | NEED  | NEED | NEED 
| resources | NEED | NEED | NEED 
| deepmd_path | NEED | 
| command |  |NEED |  NEED | 
| group_size | | NEED | NEED | 

The following table gives explicit descriptions on keys in param.json.


 Key   | Type       | Example                                                  | Discription                                                     |
| :---------------- | :--------------------- | :-------------------------------------- | :-------------------------------------------------------------|
|deepmd_path | String |"/gpfs/share/software/deepmd-kit/0.12.4/gpu/gcc/4.9.0/tf1120-lowprec" | Installed directory of DeepMD-Kit 0.x, which should contain `bin lib include`.
| python_path | String | "/gpfs/home/tzhu/anaconda3/envs/python3.6/bin/python" | Python path for DeePMD-kit 1.x installed. This option should not be used with `deepmd_path` together.
| machine | Dict | | Settings of the machine for TASK.
| resources | Dict | | Resources needed for calculation.
| # Followings are keys in resources 
| numb_node | Integer | 1 | Node count required for the job 
| task_per_node | Integer | 4 | Number of CPU cores required 
| `numb_gpu` | Integer | 4 | Number of GPUs required 
| source_list | List of string | "....../vasp.env" | Environment needed for certain job. For example, if "env" is in the list, 'source env' will be written in the script.
| module_list | List of string | [ "Intel/2018", "Anaconda3"] | For example, If "Intel/2018" is in the list, "module load Intel/2018" will be written in the script.
| time_limit | String (time format) | 23:00:00 | Maximal time permitted for the job | 
mem_limit | Interger | 16 | Maximal memory permitted to apply for the job.
| with_mpi | Boolean | true | Deciding whether to use mpi for calculation. If it's true and machine type is Slurm, "srun" will be prefixed to `command` in the script.
| qos | "string"| "bigdata" | Deciding priority, dependent on particular settings of your HPC.
| # End of resources
| command | String | "lmp_serial" | Executable path of software, such as `lmp_serial`, `lmp_mpi` and `vasp_gpu`, `vasp_std`, etc.
| group_size | Integer | 5 | DP-GEN will put these jobs together in one submitting script.
| allow_failure | Boolean | false | Allow the command returns a non-zero exit code.

## Troubleshooting

## License 
The project dpgen is licensed under [GNU LGPLv3.0](./LICENSE).





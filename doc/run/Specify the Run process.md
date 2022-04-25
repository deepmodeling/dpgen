# Specify the Run process

In the run process of the DP-GEN, we may use different softwares or versions for exploration, labeling, and training. We can specify the task as we expect in param.json. We have provided different examples of param.json in dpgen/examples/run/. 

In this section, we give a description of the param.json, taking dpgen/examples/run/dp2.x-lammps-vasp/param_CH4_deepmd-kit-2.0.1.json as an example. Here, DeePMD-kit (v2.x), LAMMPS and VASP codes are used for training, exploration and labelling respectively.

## Specify the basics

The basics related keys in param.json is given as follows

```
  "type_map": [
    "H",
    "C"
  ],
  "mass_map": [
    1,
    12
  ],
```

The basics related keys specify the information about the system. This is a param.json for a gas-phase methane molecule.

- type_map: Atom types.
- mass_map: Standard atom weights.

## Specify the init data

The init data related keys in param.json is given as follows

```
  "init_data_prefix": "....../init/",
  "init_data_sys": [
    "CH4.POSCAR.01x01x01/02.md/sys-0004-0001/deepmd"
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
```

The init data related keys specify the data for traning initial DP models and structures used for model_devi calculations. 

- init_data_prefix: Prefix of initial data directories.
- init_data_sys: Directories of initial data. You may use either absolute or relative path here.
- sys_configs_prefix: Prefix of `sys_configs`.
- sys_configs: Containing directories of structures to be explored in iterations. Wildcard characters are supported here.

## Specify the training task

The training task related keys in param.json is given as follows

```
  "numb_models": 4,
  "train_param": "input.json",
  "default_training_param": {
     "model": {
        ...
     },
     "learning_rate": {
        ...
        },
     "loss": {
        ...
     },
     "training": {
     ...
     }
  },
```

- numb_models：Number of models to be trained in `00.train`.
- default_training_param：Training parameters for `deepmd-kit` in `00.train`. You can find instructions from here: (https://github.com/deepmodeling/deepmd-kit).

## Specify the exploration task

The exploration task related keys in param.json is given as follows

```
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
```

- model_devi_dt：Timestep for MD.
- model_devi_skip：Number of structures skipped for fp in each MD.
- model_devi_f_trust_lo：Lower bound of forces for the selection. If List, should be set for each index in `sys_configs`, respectively. 
- model_devi_f_trust_hi：Upper bound of forces for the selection. If List, should be set for each index in `sys_configs`, respectively.
- model_devi_clean_traj: If type of model_devi_clean_traj is boolean type then it denote whether to clean traj folders in MD since they are too large. If it is Int type, then the most recent n iterations of traj folders will be retained, others will be removed.
- model_devi_jobs: Settings for exploration in `01.model_devi`. Each dict in the list corresponds to one iteration. The index of `model_devi_jobs` exactly accord with index of iterations.
> - sys_idx: Systems to be selected as the initial structure of MD and be explored. The index corresponds exactly to the `sys_configs`.
> - temps: Temperature (**K**) in MD.
> - press: Pressure (**Bar**) in MD.
> - trj_freq: Frequecy of trajectory saved in MD.
> - nsteps: Running steps of MD.
> - ensemble: Determining which ensemble used in MD, **options** include “npt” and “nvt”.
> - _idx: The index of iteration.

## Specify the labeling task

The labeling task related keys in param.json is given as follows

```
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
- fp_style: Software for First Principles.
- shuffle_poscar: 
- fp_task_max: Maximum of  structures to be calculated in `02.fp` of each iteration.
- fp_task_min: Minimum of structures to calculate in `02.fp` of each iteration.
- fp_pp_path: Directory of psuedo-potential file to be used for 02.fp exists. 
- fp_pp_files: Psuedo-potential file to be used for 02.fp. Note that the order of elements should correspond to the order in `type_map`. 
- fp_incar: Input file for VASP. INCAR must specify KSPACING and KGAMMA.

All the keys of the DP-GEN are explained in detail in the section Parameters.

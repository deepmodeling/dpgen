# Example-of-param.json

We have provided different examples of param.json in dpgen/examples/run/. In this section, we give a description of the param.json, taking dpgen/examples/run/dp2.x-lammps-vasp/param_CH4_deepmd-kit-2.0.1.json as an example. This is a param.json for a gas-phase methane molecule. Here, DeePMD-kit (v2.x), LAMMPS and VASP codes are used for training, exploration and labeling respectively.

## basics

The basics related keys in param.json are given as follows

```json
  "type_map": [
    "H",
    "C"
  ],
  "mass_map": [
    1,
    12
  ],
```

The basics related keys specify the basic information about the system. {dargs:argument}`type_map <run_jdata/type_map>` gives the atom types, i.e. "H" and "C". {dargs:argument}`mass_map <run_jdata/mass_map>` gives the standard atom weights, i.e. "1" and "12". 

## data

The data related keys in param.json are given as follows

```json
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

The data related keys specify the init data for training initial DP models and structures used for model_devi calculations. {dargs:argument}`init_data_prefix <run_jdata/init_data_prefix>` and {dargs:argument}`init_data_sys <run_jdata/init_data_sys>` specify the location of the init data. {dargs:argument}`sys_configs_prefix <run_jdata/sys_configs_prefix>` and {dargs:argument}`sys_configs <run_jdata/sys_configs>` specify the location of the structures. 

Here, the init data is provided at "...... /init/CH4.POSCAR.01x01x01/02.md/sys-0004-0001/deepmd". These structures are divided into two groups and provided at "....../init/CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00000*/POSCAR" and "....../init/CH4.POSCAR.01x01x01/01.scale_pert/sys-0004-0001/scale*/00001*/POSCAR". 

## training

The training related keys in param.json are given as follows

```json
  "numb_models": 4,
  "default_training_param": {
  },
```
The training related keys specify the details of training tasks. {dargs:argument}`numb_models <run_jdata/numb_models>` specifies the number of models to be trained. "default_training_param" specifies the training parameters for `deepmd-kit`. 

Here, 4 DP models will be trained in `00.train`. A detailed explanation of training parameters can be found in DeePMD-kit’s documentation (https://docs.deepmodeling.com/projects/deepmd/en/master/).

## exploration

The exploration related keys in param.json are given as follows

```json
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
The exploration related keys specify the details of exploration tasks. {dargs:argument}`model_devi_dt <run_jdata[model_devi_engine=lammps]/model_devi_dt>` specifies timestep for MD simulation. {dargs:argument}`model_devi_skip <run_jdata[model_devi_engine=lammps]/model_devi_skip>` specifies the number of structures skipped for saving in each MD. {dargs:argument}`model_devi_f_trust_lo <run_jdata[model_devi_engine=lammps]/model_devi_f_trust_lo>` and {dargs:argument}`model_devi_f_trust_hi <run_jdata[model_devi_engine=lammps]/model_devi_f_trust_hi>` specify the lower and upper bound of model_devi of forces for the selection. {dargs:argument}`model_devi_clean_traj <run_jdata[model_devi_engine=lammps]/model_devi_clean_traj>` specifies whether to clean traj folders in MD. If type of model_devi_clean_traj is boolean type then it denote whether to clean traj folders in MD since they are too large. In {dargs:argument}`model_devi_jobs <run_jdata[model_devi_engine=lammps]/model_devi_jobs>`, {dargs:argument}`sys_idx <run_jdata[model_devi_engine=lammps]/model_devi_jobs/sys_idx>` specifies the group of structures used for model_devi calculations, {dargs:argument}`temps <run_jdata[model_devi_engine=lammps]/model_devi_jobs/temps>` specifies the temperature (K) in MD, {dargs:argument}`press <run_jdata[model_devi_engine=lammps]/model_devi_jobs/press>` specifies the pressure (Bar) in MD, {dargs:argument}`trj_freq <run_jdata[model_devi_engine=lammps]/model_devi_jobs/trj_freq>` specifies the frequency of trajectory saved in MD, {dargs:argument}`nsteps <run_jdata[model_devi_engine=lammps]/model_devi_jobs/trj_freq>` specifies the running steps of MD, {dargs:argument}`ensemble <run_jdata[model_devi_engine=lammps]/model_devi_jobs/ensemble>` specifies the ensemble used in MD, and "_idx" specifies the index of iteration.

Here, MD simulations are performed at the temperature of 100 K and the pressure of 1.0 Bar with an integrator time of 2 fs under the nvt ensemble. Two iterations are set in {dargs:argument}`model_devi_jobs <run_jdata[model_devi_engine=lammps]/model_devi_jobs>`. MD simulations are run for 300 and 3000 time steps with the first and second groups of structures in {dargs:argument}`sys_configs <run_jdata/sys_configs>` in 00 and 01 iterations. We choose to save all structures generated in MD simulations and have set {dargs:argument}`trj_freq <run_jdata[model_devi_engine=lammps]/model_devi_jobs/trj_freq>` as 10, so 30 and 300 structures are saved in 00 and 01 iterations. If the "max_devi_f" of saved structure falls between 0.05 and 0.15, DP-GEN will treat the structure as a candidate. We choose to clean traj folders in MD since they are too large. If you want to save the most recent n iterations of traj folders, you can set {dargs:argument}`model_devi_clean_traj <run_jdata[model_devi_engine=lammps]/model_devi_clean_traj>` to be an integer.

## labeling 

The labeling related keys in param.json are given as follows

```json
  "fp_style": "vasp",
  "shuffle_poscar": false,
  "fp_task_max": 20,
  "fp_task_min": 1,
  "fp_pp_path": "....../methane/",
  "fp_pp_files": [
    "POTCAR"
  ],
  "fp_incar": "....../INCAR_methane"
```

The labeling related keys specify the details of labeling tasks. {dargs:argument}`fp_style <run_jdata/fp_style>` specifies software for First Principles. {dargs:argument}`fp_task_max <run_jdata/fp_task_max>` and {dargs:argument}`fp_task_min <run_jdata/fp_task_min>` specify the minimum and maximum of structures to be calculated in `02.fp` of each iteration. {dargs:argument}`fp_pp_path <run_jdata[fp_style=vasp]/fp_pp_path>` and {dargs:argument}`fp_pp_files <run_jdata[fp_style=vasp]/fp_pp_files>` specify the location of the psuedo-potential file to be used for 02.fp. {dargs:argument}`run_jdata[fp_style=vasp]/fp_incar` specifies input file for VASP. INCAR must specify KSPACING and KGAMMA.

Here, a minimum of 1 and a maximum of 20 structures will be labeled using the VASP code with the INCAR provided at "....../INCAR_methane" and POTCAR provided at "....../methane/POTCAR" in each iteration. Note that the order of elements in POTCAR should correspond to the order in {dargs:argument}`type_map <run_jdata/type_map>`. 

All the keys of the DP-GEN are explained in detail in the section Parameters.

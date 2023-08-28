## Simplify
When you have a dataset containing lots of repeated data, this step will help you simplify your dataset. The workflow contains three stages: train, model_devi, and fp. The train stage and the fp stage are as the same as the run step, and the model_devi stage will calculate model deviations of the rest data that has not been confirmed accurate. Data with small model deviations will be confirmed accurate, while the program will pick data from those with large model deviations to the new dataset.

Use the following script to start the workflow:
```bash
dpgen simplify param.json machine.json
```

Here is an example of `param.json` for QM7 dataset:
```json
{
    "type_map": [
        "C",
        "H",
        "N",
        "O",
        "S"
    ],
    "mass_map": [
        12.011,
        1.008,
        14.007,
        15.999,
        32.065
    ],
    "pick_data": "/scratch/jz748/simplify/qm7",
    "init_data_prefix": "",
    "init_data_sys": [],
    "sys_batch_size": [
        "auto"
    ],
    "numb_models": 4,
    "default_training_param": {
        "model": {
            "type_map": [
                "C",
                "H",
                "N",
                "O",
                "S"
            ],
            "descriptor": {
                "type": "se_a",
                "sel": [
                    7,
                    16,
                    3,
                    3,
                    1
                ],
                "rcut_smth": 1.00,
                "rcut": 6.00,
                "neuron": [
                    25,
                    50,
                    100
                ],
                "resnet_dt": false,
                "axis_neuron": 12
            },
            "fitting_net": {
                "neuron": [
                    240,
                    240,
                    240
                ],
                "resnet_dt": true
            }
        },
        "learning_rate": {
            "type": "exp",
            "start_lr": 0.001,
            "stop_lr": 5e-8,
            "decay_rate": 0.99
        },
        "loss": {
            "start_pref_e": 0.02,
            "limit_pref_e": 1,
            "start_pref_f": 1000,
            "limit_pref_f": 1,
            "start_pref_v": 0,
            "limit_pref_v": 0,
            "start_pref_pf": 0,
            "limit_pref_pf": 0
        },
        "training": {
            "set_prefix": "set",
            "numb_steps": 10000,
            "disp_file": "lcurve.out",
            "disp_freq": 1000,
            "numb_test": 1,
            "save_freq": 1000,
            "disp_training": true,
            "time_training": true,
            "profiling": false,
            "profiling_file": "timeline.json"
        },
        "_comment": "that's all"
    },
    "fp_style": "gaussian",
    "shuffle_poscar": false,
    "fp_task_max": 1000,
    "fp_task_min": 10,
    "fp_pp_path": "/home/jzzeng/",
    "fp_pp_files": [],
    "fp_params": {
        "keywords": "mn15/6-31g** force nosymm scf(maxcyc=512)",
        "nproc": 28,
        "multiplicity": 1,
        "_comment": " that's all "
    },
    "init_pick_number":100,
    "iter_pick_number":100,
    "model_devi_f_trust_lo":0.25,
    "model_devi_f_trust_hi":0.45,
    "_comment": " that's all "
}
```

Here {dargs:argument}`pick_data <simplify_jdata/pick_data>` is the directory to data to simplify where the program recursively detects systems `System` with `deepmd/npy` format. {dargs:argument}`init_pick_number <simplify_jdata/init_pick_number>` and {dargs:argument}`iter_pick_number <simplify_jdata/iter_pick_number>` are the numbers of picked frames. {dargs:argument}`model_devi_f_trust_lo <simplify_jdata/model_devi_f_trust_lo>` and {dargs:argument}`model_devi_f_trust_hi <simplify_jdata/model_devi_f_trust_hi>` mean the range of the max deviation of atomic forces in a frame. {dargs:argument}`fp_style <simplify_jdata/fp_style>` can be either `gaussian` or `vasp` currently. Other parameters are as the same as those of generator.

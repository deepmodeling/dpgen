## Troubleshooting
1. The most common problem is whether two settings correspond with each other, including:
    - The order of elements in `type_map` and `mass_map` and **`fp_pp_files`**.
    - Size of `init_data_sys` and `init_batch_size`.
    - Size of `sys_configs` and `sys_batch_size`.
    - Size of `sel_a` and actual types of atoms in your system.
    - Index of `sys_configs` and `sys_idx`.

2. Please verify the directories of `sys_configs`. If there isn't any POSCAR for `01.model_devi` in one iteration, it may happen that you write the false path of `sys_configs`.
3. Correct format of JSON file.
4. The frames of one system should be larger than `batch_size` and `numb_test` in `default_training_param`. It happens that one iteration adds only a few structures and causes error in next iteration's training. In this condition, you may let `fp_task_min` be larger than `numb_test`.


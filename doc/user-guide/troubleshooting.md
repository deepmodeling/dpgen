## Troubleshooting
1. The most common problem is whether two settings correspond with each other, including:
    - The order of elements in `type_map` and `mass_map` and **`fp_pp_files`**.
    - Size of `init_data_sys` and `init_batch_size`.
    - Size of `sys_configs` and `sys_batch_size`.
    - Size of `sel_a` and actual types of atoms in your system.
    - Index of `sys_configs` and `sys_idx`.

2. Please verify the directories of `sys_configs`. If there isn't any POSCAR for `01.model_devi` in one iteration, it may happen that you write the false path of `sys_configs`. Note that `init_data_sys` is a list, while `sys_configs` should be a two-dimensional list. The first dimension corresponds to `sys_idx`, and the second level are some poscars under each group. Refer to the [sample file](github.com/deepmodeling/dpgen/blob/master/examples/run/dp2.x-lammps-vasp/param_CH4_deepmd-kit-2.0.1.json ). 

3. Correct format of JSON file.

4. The frames of one system should be larger than `batch_size` and `numb_test` in `default_training_param`. It happens that one iteration adds only a few structures and causes error in next iteration's training. In this condition, you may let `fp_task_min` be larger than `numb_test`.

5. If you found the dpgen with the same version on two machines behaves differently, you may have modified the code in one of them.


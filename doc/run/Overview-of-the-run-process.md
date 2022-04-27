# Overview of the Run process

The run process contains a series of successive iterations, succussively undertaken in order such as heating the system to certain temperature. Each iteration is composed of three steps: exploration, labeling, and training. Accordingly, there are three sub-folders: 00.train, 01.model_devi, and 02.fp in each iteration.

00.train: DP-GEN will train several (default 4) models based on initial and generated data. The only difference between these models is the random seed for neural network initialization.

01.model_devi : represent for model-deviation. DP-GEN will use models obtained from 00.train to run Molecular Dynamics(default LAMMPS). Larger deviation for structure properties (default is force of atoms) means less accuracy of the models. Using this criterion, a few fructures will be selected and put into next stage 02.fp for more accurate calculation based on First Principles.

02.fp : Selected structures will be calculated by first principles methods(default VASP). DP-GEN will obtain some new data and put them together with initial data and data generated in previous iterations. After that a new training will be set up and DP-GEN will enter next iteration!

In the run process of the DP-GEN, we need to specify the basic information about the system, the initial data, details of the training, exploration, and labeling tasks. In addition, we need to specify the software, machine environment, and computing resource and enable the process of job generation, submission, query, and collection automatically. We can perform the run process as we expect by specifying the keywords in param.json and machine.json, and they will be introduced in detail in the following sections. 

Here, we give a general description of the run process. We can execute the run process of DP-GEN easily by:

```
dpgen run param.json machine.json
```

The following files or folders will be created and upgraded by codes：

- iter.00000x contains the main results that DP-GEN generates in the first iteration.
- record.dpgen records the current stage of the run process.
- dpgen.log includes time and iteration information.

When the first iteration is completed, the folder structure of iter.000000 is like this:

```
$ ls iter.000000
00.train 01.model_devi 02.fp
```

In folder iter.000000/ 00.train:

- Folder 00x contains the input and output files of the DeePMD-kit, in which a model is trained.
- graph.00x.pb is the model DeePMD-kit generates. The only difference between these models is the random seed for neural network initialization.

In folder iter.000000/ 01.model_devi：

- Folder confs contains the initial configurations for LAMMPS MD converted from POSCAR you set in "sys_configs" of param.json. 
- Folder task.000.00000x contains the input and output files of the LAMMPS. In folder task.000.00000x, file model_devi.out records the model deviation of concerned labels, energy and force in MD. It serves as the criterion for selecting which structures and doing first-principle calculations.

In folder iter.000000/ 02.fp：

- candidate.shuffle.000.out records which structures will be selected from last step 01.model_devi.  There are always far more candidates than the maximum you expect to calculate at one time. In this condition, DP-GEN will randomly choose up to `"fp_task_max"` structures and form the folder task.*.
- rest_accurate.shuffle.000.out records the other structures where our model is accurate ("max_devi_f" is less than `"model_devi_f_trust_lo"`, no need to calculate any more), 
- rest_failed.shuffled.000.out records the other structures where our model is too inaccurate (lager than `"model_devi_f_trust_hi"`, there may be some error).
- data.000: After first-principle calculations, DP-GEN will collect these data and change them into the format DeePMD-kit needs. In the next iteration's 00.train, these data will be trained together as well as initial data.

DP-GEN identifies the stage of run process by a record file, record.dpgen, which will be created and upgraded by codes. Each line contains two number: the first is index of iteration, and the second ,ranging from 0 to 9 ,records which stage in each iteration is currently running.

| Index of iterations  | Stage in eachiteration      | Process          |
|:---------------------|:----------------------------|:-----------------|
| 0                    | 0                           | make_train       |
| 0                    | 1                           | run_train        |
| 0                    | 2                           | post_train       |
| 0                    | 3                           | make_model_devi  |
| 0                    | 4                           | run_model_devi   |
| 0                    | 5                           | post_model_devi  |
| 0                    | 6                           | make_fp          |
| 0                    | 7                           | run_fp           |
| 0                    | 8                           | post_fp          |

0,1,2 correspond to make_train, run_train, post_train. DP-GEN will write scripts in make_train, run the task by specific machine in run_train and collect result in post_train. The records for model_devi and fp stage follow similar rules.

If the process of DP-GEN stops for some reason, DP-GEN will automatically recover the main process by record.dpgen. You may also change it manually for your purpose, such as removing the last iterations and recovering from one checkpoint.

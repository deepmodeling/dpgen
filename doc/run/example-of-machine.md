# Example of machine.json

## DPDispatcher Update Note

DPDispatcher has updated and the api of machine.json is changed. DP-GEN will use the new DPDispatcher if the value of key "api_version" in machine.json is equal to or large than 1.0. And for now, DPDispatcher is maintained on a separate repo (https://github.com/deepmodeling/dpdispatcher). Please check the documents (https://deepmd.readthedocs.io/projects/dpdispatcher/en/latest/) for more information about the new DPDispatcher. 

DP-GEN will use the old DPDispatcher if the key "api_version" is not specified in machine.json or the "api_version" is smaller than 1.0. This gurantees that the old machine.json still works.

## New DPDispatcher

Each iteration in the run process of DP-GEN is composed of three steps: exploration, labeling, and training. Accordingly, machine.json is composed of three parts: train, model_devi, and fp. Each part is a list of dicts. Each dict can be considered as an independent environment for calculation. 

In this section, we will show you how to perform train task at a local workstation, model_devi task at a local Slurm cluster, and fp task at a remote PBS cluster using the new DPDispatcher. For each task, three types of keys are needed:
- Command: provides the command used to execute each step.
- Machine: specifies the machine environment (local workstation, local or remote cluster, or cloud server).
- Resources: specify the number of groups, nodes, CPU, and GPU; enable the virtual environment.

### Performing train tasks at a local workstation

In this example, we perform the `train` task on a local workstation.

```json
"train": [
    {
      "command": "dp",
      "machine": {
        "batch_type": "Shell",
        "context_type": "local",
        "local_root": "./",
        "remote_root": "/home/user1234/work_path"
      },
      "resources": {
        "number_node": 1,
        "cpu_per_node": 4,
        "gpu_per_node": 1,
        "group_size": 1,
        "source_list": ["/home/user1234/deepmd.env"]
      }
    }
  ],
```

The "command" for the train task in the DeePMD-kit is "dp".

In machine parameters, "batch_type" specifies the type of job scheduling system. If there is no job scheduling system, we can use the "Shell" to perform the task. "context_type" specifies the method of data transfer, and "local" means copying and moving data via local file storage systems (e.g. cp, mv, etc.). In DP-GEN, the paths of all tasks are automatically located and set by the software, and therefore "local_root" is always set to "./". The input file for each task will be sent to the "remote_root" and the task will be performed there, so we need to make sure that the path exists.

In the resources parameter, "number_node", "cpu_per_node", and "gpu_per_node" specify the number of nodes, the number of CPUs, and the number of GPUs required for a task respectively. "group_size", which needs to be highlighted, specifies how many tasks will be packed into a group. In the training tasks, we need to train 4 models. If we only have one GPU, we can set the "group_size" to 4. If "group_size" is set to 1, 4  models will be trained on one GPU at the same time, as there is no job scheduling system. Finally, the environment variables can be activated by "source_list". In this example, "source /home/user1234/deepmd.env" is executed before "dp" to load the environment variables necessary to perform the training task.

### Perform model_devi task at a local Slurm cluster

In this example, we perform the model_devi task at a local Slurm workstation.

```json
"model_devi": [
    {
      "command": "lmp",
      "machine": {
       "context_type": "local",
        "batch_type": "Slurm",
        "local_root": "./",
        "remote_root": "/home/user1234/work_path"
      },
      "resources": {
        "number_node": 1,
        "cpu_per_node": 4,
        "gpu_per_node": 1,
        "queue_name": "QueueGPU",
        "custom_flags" : ["#SBATCH --mem=32G"],
        "group_size": 10,
        "source_list": ["/home/user1234/lammps.env"]
      }
    }
],
```

The "command" for the model_devi task in the LAMMPS is "lmp".

In the machine parameter, we specify the type of job scheduling system by changing the "batch_type" to "Slurm".

In the resources parameter, we specify the name of the queue to which the task is submitted by adding "queue_name". We can add additional lines to the calculation script via the "custom_flags". In the model_devi steps, there are frequently many short tasks, so we usually pack multiple tasks (e.g. 10) into a group for submission. Other parameters are similar to that of the local workstation.

### Perform fp task in a remote PBS cluster

In this example, we perform the fp task at a remote PBS cluster that can be accessed via SSH.

```json
"fp": [
    {
      "command": "mpirun -n 32 vasp_std",
      "machine": {
       "context_type": "SSHContext",
        "batch_type": "PBS",
        "local_root": "./",
        "remote_root": "/home/user1234/work_path",
        "remote_profile": {
          "hostname": "39.xxx.xx.xx",
          "username": "user1234"
         }
      },
      "resources": {
        "number_node": 1,
        "cpu_per_node": 32,
        "gpu_per_node": 0,
        "queue_name": "QueueCPU",
        "group_size": 5,
        "source_list": ["/home/user1234/vasp.env"]
      }
    }
],
```

VASP code is used for fp task and mpi is used for parallel computing, so "mpirun -n 32" is added to specify the number of parallel threads.

In the machine parameter, "context_type" is modified to "SSHContext" and "batch_type" is modified to "PBS". It is worth noting that "remote_root" should be set to an accessible path on the remote PBS cluster. "remote_profile" is added to specify the information used to connect the remote cluster, including hostname, username,  port, etc. 

In the resources parameter, we set "gpu_per_node" to 0 since it is cost-effective to use the CPU for VASP calculations.

Explicit descriptions of keys in machine.json will be given in the following section.

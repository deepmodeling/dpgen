## Relaxation run

The work path of each task should be in the form like `confs/mp-*/relaxation` and all task is in the form like `confs/mp-*/relaxation/relax_task`. 

The `machine.json` file should be applied in this process and the machine parameters (eg. GPU or CPU) are determined according to the task type (VASP or LAMMPS). Then in each work path, the corresponding tasks would be submitted and the results would be sent back through [make_dispatcher](https://github.com/deepmodeling/dpgen/blob/devel/dpgen/dispatcher/Dispatcher.py).

Take `deepmd` run for example:
```bash
nohup dpgen autotest run relaxation.json machine-ali.json > run.result 2>&1 &
tree confs/std-fcc/relaxation/
```

the output would be:
```
confs/std-fcc/relaxation/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- jr.json
`-- relax_task
    |-- conf.lmp
    |-- dump.relax
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- log.lammps
    |-- outlog
    |-- POSCAR -> ../../POSCAR
    `-- task.json
```
`dump.relax` is the file storing configurations and `log.lammps` is the output file for lammps.

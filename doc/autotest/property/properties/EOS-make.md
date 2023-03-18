## EOS make

**Step 1.** Before `make` in EOS, the equilibrium configuration `CONTCAR` must be present in `confs/mp-*/relaxation`.

**Step 2.** For the input example in the previous section, when we do `make`, 40 tasks would be generated as `confs/mp-*/eos_00/task.000000, confs/mp-*/eos_00/task.000001, ... , confs/mp-*/eos_00/task.000039`. The suffix `00` is used  for possible `refine` later.

**Step 3.** If the task directory, for example `confs/mp-*/eos_00/task.000000` is not empty, the old input files in it including `INCAR`, `POSCAR`, `POTCAR`, `conf.lmp`, `in.lammps` would be deleted.

**Step 4.** In each task directory, `POSCAR.orig` would link to `confs/mp-*/relaxation/CONTCAR`. Then the `scale` parameter can be calculated as:

```txt
scale = (vol_current / vol_equi) ** (1. / 3.)
```

`vol_current` is the corresponding volume per atom of the current task and `vol_equi` is the volume per atom of the equilibrium configuration. Then the `poscar_scale` function in `dpgen.auto_test.lib.vasp` module would help to generate `POSCAR` file with `vol_current` in `confs/mp-*/eos_00/task.[0-9]*[0-9]`.

**Step 5.** According to the task type, the input file including `INCAR`, `POTCAR` or `conf.lmp`, `in.lammps` would be written in every `confs/mp-*/eos_00/task.[0-9]*[0-9]`.

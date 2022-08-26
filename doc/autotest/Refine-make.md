## Refine-make

```bash
dpgen autotest make refine.json
tree confs/std-fcc/vacancy_01/
```
the output will be:

```
confs/std-fcc/vacancy_01/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
`-- task.000000
    |-- conf.lmp
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- POSCAR -> ../../vacancy_00/task.000000/CONTCAR
    |-- supercell.json -> ../../vacancy_00/task.000000/supercell.json
    `-- task.json
```

an new directory `vacancy_01` would be established and the starting configuration links to previous results.

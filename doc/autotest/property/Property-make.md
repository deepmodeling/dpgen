## Property make

```bash
dpgen autotest make property.json
```
**EOS output:**

```
confs/std-fcc/eos_00/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- task.000000
|   |-- conf.lmp
|   |-- eos.json
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps
|   |-- inter.json
|   |-- POSCAR
|   |-- POSCAR.orig -> ../../relaxation/relax_task/CONTCAR
|   `-- task.json
|-- task.000001
|   |-- conf.lmp
|   |-- eos.json
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps
|   |-- inter.json
|   |-- POSCAR
|   |-- POSCAR.orig -> ../../relaxation/relax_task/CONTCAR
|   `-- task.json
...
`-- task.000019
    |-- conf.lmp
    |-- eos.json
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps
    |-- inter.json
    |-- POSCAR
    |-- POSCAR.orig -> ../../relaxation/relax_task/CONTCAR
    `-- task.json
```

`eos.json` records the `volume` and `scale` of the corresponding task.

**Elastic output:**

```
confs/std-fcc/elastic_00/
|-- equi.stress.json
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- POSCAR -> ../relaxation/relax_task/CONTCAR
|-- task.000000
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- POSCAR
|   |-- strain.json
|   `-- task.json
|-- task.000001
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- POSCAR
|   |-- strain.json
|   `-- task.json
...
`-- task.000023
    |-- conf.lmp
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- POSCAR
    |-- strain.json
    `-- task.json
```

`equi.stress.json` records the stress information of the equilibrium task and `strain.json` records the deformation information of the corresponding task.

**Vacancy output:**

```
confs/std-fcc/vacancy_00/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- POSCAR -> ../relaxation/relax_task/CONTCAR
`-- task.000000
    |-- conf.lmp
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- POSCAR
    |-- supercell.json
    `-- task.json
```
`supercell.json` records the supercell size information of the corresponding task.

**Interstitial output:**

```
confs/std-fcc/interstitial_00/
|-- element.out
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- POSCAR -> ../relaxation/relax_task/CONTCAR
|-- task.000000
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- POSCAR
|   |-- supercell.json
|   `-- task.json
`-- task.000001
    |-- conf.lmp
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- POSCAR
    |-- supercell.json
    `-- task.json
```

`element.out` records the inserted element type of each task and `supercell.json` records the supercell size information of the corresponding task.

**Surface output:**

```
confs/std-fcc/surface_00/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- POSCAR -> ../relaxation/relax_task/CONTCAR
|-- task.000000
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- miller.json
|   |-- POSCAR
|   |-- POSCAR.tmp
|   `-- task.json
|-- task.000001
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- miller.json
|   |-- POSCAR
|   |-- POSCAR.tmp
|   `-- task.json
...
`-- task.000008
    |-- conf.lmp
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- miller.json
    |-- POSCAR
    |-- POSCAR.tmp
    `-- task.json
```

`miller.json` records the miller index of the corresponding task.

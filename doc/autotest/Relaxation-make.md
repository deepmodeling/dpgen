## Relaxation-make

The list of the directories storing structures are `["confs/std-*"]` in the previous example. For single element system, if `POSCAR` doesn't exist in the directories: `std-fcc`, `std-hcp`, `std-dhcp`, `std-bcc`, `std-diamond`, and `std-sc`, the package will automatically generate the standard crystal structures **`fcc`**, **`hcp`**, **`dhcp`**, **`bcc`**, **`diamond`**, and **`sc`** in the corresponding directories, respectively. In other conditions and for multi-component system (more than 1), if `POSCAR` doesn't exist, the package will terminate and print the error **"no configuration for autotest"**.

### VASP relaxation
Take the input example of Al in the previous section, when we do `make` as follows:
```bash
dpgen autotest make relaxation.json
```
the following files would be generated:
```bash
tree confs/std-fcc/relaxation/
```

```
confs/std-fcc/relaxation/
|-- INCAR
|-- POTCAR
`-- relax_task
    |-- INCAR -> ../INCAR
    |-- inter.json
    |-- KPOINTS
    |-- POSCAR -> ../../POSCAR
    |-- POTCAR -> ../POTCAR
    `-- task.json
```
`inter.json` records the information in the `interaction` dictionary and `task.json` records the information in the `relaxation` dictionary.

### LAMMPS relaxation
```bash
dpgen autotest make relaxation.json
tree confs/std-fcc/
```
the output would be:
```
confs/std-fcc/
|-- POSCAR
`-- relaxation
    |-- frozen_model.pb -> ../../../frozen_model.pb
    |-- in.lammps
    `-- relax_task
        |-- conf.lmp
        |-- frozen_model.pb -> ../frozen_model.pb
        |-- in.lammps -> ../in.lammps
        |-- inter.json
        |-- POSCAR -> ../../POSCAR
        `-- task.json
```
the `conf.lmp` is the input configuration and `in.lammps` is the input command file for lammps.

**in.lammps**: the package would generate the file `confs/mp-*/relaxation/in.lammps` as follows and we refer the user to the further information of [fix box/relax](https://lammps.sandia.gov/doc/fix_box_relax.html) function in lammps:

```txt
clear
units 	          metal
dimension	  3
boundary	  p p p
atom_style	  atomic
box               tilt large
read_data         conf.lmp
mass              1 26.982
neigh_modify      every 1 delay 0 check no
pair_style deepmd frozen_model.pb
pair_coeff
compute           mype all pe
thermo            100
thermo_style      custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype
dump              1 all custom 100 dump.relax id type xs ys zs fx fy fz
min_style         cg
fix               1 all box/relax iso 0.0
minimize          1.000000e-12 1.000000e-06 5000 500000
fix               1 all box/relax aniso 0.0
minimize          1.000000e-12 1.000000e-06 5000 500000
variable          N equal count(all)
variable          V equal vol
variable          E equal "c_mype"
variable          tmplx equal lx
variable          tmply equal ly
variable          Pxx equal pxx
variable          Pyy equal pyy
variable          Pzz equal pzz
variable          Pxy equal pxy
variable          Pxz equal pxz
variable          Pyz equal pyz
variable          Epa equal ${E}/${N}
variable          Vpa equal ${V}/${N}
variable          AA equal (${tmplx}*${tmply})
print "All done"
print "Total number of atoms = ${N}"
print "Final energy per atoms = ${Epa}"
print "Final volume per atoms = ${Vpa}"
print "Final Base area = ${AA}"
print "Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}"
```

If user provides lammps input command file `in.lammps`, the `thermo_style` and `dump` commands should be the same as the above file.

**interatomic potential model**: the `frozen_model.pb` in `confs/mp-*/relaxation` would link to the `frozen_model.pb` file given in the input.





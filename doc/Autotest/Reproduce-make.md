
```bash
dpgen autotest make reproduce.json
tree confs/std-fcc/interstitial_reprod/
```
the output will be:

```
confs/std-fcc/interstitial_reprod/
|-- frozen_model.pb -> ../../../frozen_model.pb
|-- in.lammps
|-- task.000000
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- POSCAR
|   `-- task.json
|-- task.000001
|   |-- conf.lmp
|   |-- frozen_model.pb -> ../frozen_model.pb
|   |-- in.lammps -> ../in.lammps
|   |-- inter.json
|   |-- POSCAR
|   `-- task.json
...
`-- task.000038
    |-- conf.lmp
    |-- frozen_model.pb -> ../frozen_model.pb
    |-- in.lammps -> ../in.lammps
    |-- inter.json
    |-- POSCAR
    `-- task.json
```

every singe frame in the initial data is split into each task and the following `in.lammps` would help to do the `static` calculation:

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
run               0
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
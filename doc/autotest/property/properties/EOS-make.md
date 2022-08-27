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

For EOS calculations by VASP, if `change_box` is `True`, `ISIF` in VASP would be 4, else `ISIF` would be 2. The default value of `change_box` is `True`. For further information of the use of `ISIF` in VASP, we refer users to [ISIF command](https://www.vasp.at/wiki/index.php/ISIF).

For EOS calculations by LAMMPS, when `change_box` is `True`, an example of `in.lammps` for AlMg is given as below and the `scale` parameter in line 5 is calculated by the equation above.

```txt
clear
variable        GPa2bar	equal 1e4
variable        B0	equal 70
variable        bp	equal 0
variable	xx	equal scale
variable        yeta	equal 1.5*(${bp}-1)
variable        Px0	equal 3*${B0}*(1-${xx})/${xx}^2*exp(${yeta}*(1-${xx}))
variable        Px	equal ${Px0}*${GPa2bar}
units           metal
dimension       3
boundary	p p p
atom_style	atomic
box             tilt large
read_data       conf.lmp
mass            1 1
mass            2 1
neigh_modify    every 1 delay 0 check no
pair_style      deepmd frozen_model.pb
pair_coeff
compute         mype all pe
thermo          100
thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype
dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz
min_style       cg
fix             1 all box/relax iso ${Px}
minimize        1.000000e-12 1.000000e-06 5000 500000
fix             1 all box/relax aniso ${Px}
minimize        1.000000e-12 1.000000e-06 5000 500000
variable        N equal count(all)
variable        V equal vol
variable        E equal "c_mype"
variable        Pxx equal pxx
variable        Pyy equal pyy
variable        Pzz equal pzz
variable        Pxy equal pxy
variable        Pxz equal pxz
variable        Pyz equal pyz
variable        Epa equal ${E}/${N}
variable        Vpa equal ${V}/${N}
print "All done"
print "Total number of atoms  = ${N}"
print "Relax at Press         = ${Px} Bar"
print "Final energy per atoms = ${Epa} eV"
print "Final volume per atoms = ${Vpa} A^3"
print "Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}"
```

when `change_box` is `False`, an example of `in.lammps` for AlMg is given as:
```txt
clear
units 	metal
dimension	3
boundary	p	p    p
atom_style	atomic
box         tilt large
read_data   conf.lmp
mass            1 1
mass            2 1
neigh_modify    every 1 delay 0 check no
pair_style deepmd frozen_model.pb
pair_coeff
compute         mype all pe
thermo          100
thermo_style    custom step pe pxx pyy pzz pxy pxz pyz lx ly lz vol c_mype
dump            1 all custom 100 dump.relax id type xs ys zs fx fy fz
min_style       cg
minimize        1.000000e-12 1.000000e-06 5000 500000
variable        N equal count(all)
variable        V equal vol
variable        E equal "c_mype"
variable        tmplx equal lx
variable        tmply equal ly
variable        Pxx equal pxx
variable        Pyy equal pyy
variable        Pzz equal pzz
variable        Pxy equal pxy
variable        Pxz equal pxz
variable        Pyz equal pyz
variable        Epa equal ${E}/${N}
variable        Vpa equal ${V}/${N}
variable        AA equal (${tmplx}*${tmply})
print "All done"
print "Total number of atoms = ${N}"
print "Final energy per atoms = ${Epa}"
print "Final volume per atoms = ${Vpa}"
print "Final Base area = ${AA}"
print "Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}"
```


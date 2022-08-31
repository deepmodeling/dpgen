## Elastic make

**Step 1.** The `DeformedStructureSet` module in [pymatgen.analysis.elasticity.strain](https://pymatgen.org/pymatgen.analysis.elasticity.strain.html) is used to generate a set of independently deformed structures. `equi.stress.out` file is written to record the equilibrium stress in the Elastic directory. For the example in the previous section, `equi.stress.out` should be in `confs/mp-*/elastic_00`.

**Step 2.** If there are `init_from_suffix` and `output_suffix` parameter in the `properties` part, the [refine process](../../refine/Refine-get-started-and-input-examples) follows. Else, the deformed structure (`POSCAR`) and strain information (`strain.out`) are written in the task directory, for example, in `confs/mp-*/elastic_00/task.000000`.

**Step 3.** When doing `elastic` by VASP, `ISIF=2`. When doing by LAMMPS, the following `in.lammps` would be written.

```txt
units 	metal
dimension	3
boundary	p	p    p
atom_style	atomic
box             tilt large
read_data       conf.lmp
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
minimize        0 1.000000e-10 5000 500000
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
print "Total number of atoms = ${N}"
print "Final energy per atoms = ${Epa}"
print "Final volume per atoms = ${Vpa}"
print "Final Stress (xx yy zz xy xz yz) = ${Pxx} ${Pyy} ${Pzz} ${Pxy} ${Pxz} ${Pyz}"
```

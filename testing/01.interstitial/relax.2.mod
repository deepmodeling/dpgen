#----Relax system-----------------------------------------------
reset_timestep	0

thermo		${thermo_freq}
thermo_style	custom step pe lx ly lz press pxx pyy pzz c_pot

dump		1 all custom ${dump_freq} dump.relax.2 id type xs ys zs 

fix		2 fnull setforce 0 0 0
min_style	cg
minimize	${mini_etol} ${mini_ftol} 5000 5000
unfix		2


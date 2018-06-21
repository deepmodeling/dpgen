#----Eenergy Computation----------------------------------------
compute		pot all pe

#----Relax system-----------------------------------------------
reset_timestep	0

thermo		${thermo_freq}
thermo_style	custom step pe lx ly lz press pxx pyy pzz c_pot

dump		1 all custom ${dump_freq} dump.relax id type xs ys zs 

fix		1 all box/relax iso 0.0 
min_style	cg
minimize	${mini_etol} ${mini_ftol} 5000 5000

run		0
undump		1


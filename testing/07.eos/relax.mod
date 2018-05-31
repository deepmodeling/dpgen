#----Relax system-----------------------------------------------
variable	mini_etol	equal 1e-12
variable	mini_ftol	equal 1e-6
reset_timestep	0
fix		1 all box/relax iso 0.0 
min_style	cg
minimize	${mini_etol} ${mini_ftol} 5000 5000
unfix		1


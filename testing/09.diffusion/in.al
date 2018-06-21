# melt the fcc lattice first and then do stat
variable	test_temp	equal 933.47
variable	temp_damp	equal 0.5
variable	test_pres	equal 1.0
variable	pres_damp	equal 0.5
variable	liquid_temp	equal 1500
variable        rlx_step        equal 20000

variable        dumpfreq equal 100 # dump snapshot frequency
variable        thermofreq equal 100 # thermo output frequency
variable        dt equal 0.001
variable        nsteps equal 200000

dimension 	3
boundary	p p p
units		metal

atom_style	atomic
neighbor	1.0 bin
neigh_modify	every 10 delay 0 check no

variable	numb_x equal 8
variable	numb_y equal 8
variable	numb_z equal 8
variable	latt_a equal 4.04
variable	protct equal 0.001
variable	size_x equal ${numb_x}*${latt_a}-${protct}
variable	size_y equal ${numb_y}*${latt_a}-${protct}
variable	size_z equal ${numb_z}*${latt_a}-${protct}

lattice		fcc ${latt_a}
region		whole block -${protct} ${size_x} -${protct} ${size_y} -${protct} ${size_z} units box
create_box	1 whole
create_atoms	1 region whole
mass		1 27

include		potential.mod

velocity        all create ${liquid_temp} 12345 loop geom
timestep        ${dt}

fix             melt all npt temp ${liquid_temp} ${liquid_temp} ${temp_damp} iso ${test_pres} ${test_pres} ${pres_damp}
run ${rlx_step}

unfix melt
fix             equil all npt temp ${test_temp} ${test_temp} ${temp_damp} iso ${test_pres} ${test_pres} ${pres_damp}
run ${rlx_step}
unfix equil

reset_timestep 0
fix             do_stat all nve
#npt temp ${test_temp} ${test_temp} ${temp_damp} iso ${test_pres} ${test_pres} ${pres_damp}
compute         msd all msd com yes
variable        twopoint equal c_msd[4]/6/(step*dt+1.0e-6)*10
fix             9 all vector 10 c_msd[4]
variable        fitslope equal slope(f_9)/6/(dt)
compute         2 all vacf
fix             5 all vector 1 c_2[4]
variable        diff equal dt*trap(f_5)

thermo_style custom step time  temp ke etotal press c_msd[4] v_twopoint v_fitslope v_diff

thermo          ${thermofreq}
thermo_modify flush yes

run            ${nsteps}

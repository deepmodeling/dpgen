#----Input parameters-------------------------------------------
variable	ao		equal 4.04
variable	ncopy		equal 1
variable	thermo_freq	equal 5
variable	dump_freq	equal 50
variable	mini_etol	equal 1e-10
variable	mini_ftol	equal 1e-10

#----Create system----------------------------------------------
clear
units 		metal
dimension	3
boundary	p	p    p      
atom_style	atomic

variable	hnc equal ${ncopy}
variable	mhnc equal 0
lattice         fcc ${ao}
region		simbox block ${mhnc} ${hnc} ${mhnc} ${hnc} ${mhnc} ${hnc}
create_box	1 simbox
lattice 	fcc ${ao} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
create_atoms	1 region simbox
mass		1 27

group		fnull empty

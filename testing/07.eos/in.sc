units		metal
dimension	3
boundary	p p p
atom_style	atomic
variable	latt_a		equal LATT_A

lattice		sc ${latt_a}
region		whole block 0 1 0 1 0 1
create_box	1 whole
create_atoms	1 region whole
mass		1 27

include		potential.mod

thermo		1
thermo_style	custom step pe pxx pyy pzz lx ly lz

#include		relax.mod

run		0
variable	natoms equal count(all)
variable	petotal equal "pe"
variable	latt_a equal "lx"
variable	epa equal ${petotal}/${natoms}
variable	vpa equal ${latt_a}*${latt_a}*${latt_a}/${natoms}

print "ENER_PER_ATOM ${epa}"
print "VOLM_PER_ATOM ${vpa}"




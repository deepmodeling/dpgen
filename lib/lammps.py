#!/usr/bin/env python3

import random, os, sys
import subprocess as sp

def cvt_lammps_conf (fin, 
                     fout) :
    thisfile = os.path.abspath(__file__)
    thisdir = os.path.dirname(thisfile)
    cmd = os.path.join(thisdir, 'ovito_file_convert.py')
    sp.check_call([cmd, fin, fout])    

def make_lammps_input(ensemble, 
                      conf_file,
                      graphs,
                      nsteps, 
                      trj_freq, 
                      mass_map,
                      temp, 
                      pres = None, 
                      max_seed = 1000000) :
    ret = "variable        NSTEPS          equal %d\n" % nsteps
    ret+= "variable        THERMO_FREQ     equal %d\n" % trj_freq
    ret+= "variable        DUMP_FREQ       equal %d\n" % trj_freq
    ret+= "variable        TEMP            equal %f\n" % temp
    ret+= "variable        PRES            equal %f\n" % pres
    ret+= "\n"
    ret+= "units           metal\n"
    ret+= "boundary        p p p\n"
    ret+= "atom_style      atomic\n"
    ret+= "\n"
    ret+= "neighbor        1.0 bin\n"
    ret+= "neigh_modify    every 10\n"
    ret+= "\n"
    ret+= "read_data       %s\n" % conf_file    
    for jj in range(len(mass_map)) :
        ret+= "mass            %d %f\n" %(jj+1, mass_map[jj])
    graph_list = ""
    for ii in graphs :
        graph_list += ii + " "
    ret+= "pair_style      deepmd %s ${THERMO_FREQ} model_devi.out\n" % graph_list
    ret+= "pair_coeff      \n"
    ret+= "\n"
    ret+= "thermo_style    custom step temp pe ke etotal press vol lx ly lz xy xz yz\n"
    ret+= "thermo          ${THERMO_FREQ}\n"
    ret+= "dump            1 all custom ${DUMP_FREQ} traj/*.lammpstrj id type x y z\n"
    ret+= "\n"
    ret+= "velocity        all create ${TEMP} %d" % random.randrange(max_seed)
    ret+= "\n"
    if ensemble == "npt" :
        assert (pres is not None)
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} 0.5 tri ${PRES} ${PRES} 3.0\n"
    elif ensemble == "nvt" :
        ret+= "fix             1 all nvt temp ${TEMP} ${TEMP} 0.5\n"
    else :
        raise RuntimeError("unknown emsemble " + ensemble)
    ret+= "\n"
    ret+= "timestep        0.002\n"
    ret+= "run             ${NSTEPS}\n"
    return ret
        
# ret = make_lammps_input ("npt", "al.lmp", ['graph.000.pb', 'graph.001.pb'], 20000, 20, [27], 1000, pres = 1.0)
# print (ret)
# cvt_lammps_conf('POSCAR', 'tmp.lmp')


    
 

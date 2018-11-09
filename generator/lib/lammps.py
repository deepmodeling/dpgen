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
                      dt,
                      neidelay,
                      trj_freq, 
                      mass_map,
                      temp, 
                      tau_t = 0.1,
                      pres = None,
                      tau_p = 0.5,
                      max_seed = 1000000) :
    ret = "variable        NSTEPS          equal %d\n" % nsteps
    ret+= "variable        THERMO_FREQ     equal %d\n" % trj_freq
    ret+= "variable        DUMP_FREQ       equal %d\n" % trj_freq
    ret+= "variable        TEMP            equal %f\n" % temp
    ret+= "variable        PRES            equal %f\n" % pres
    ret+= "variable        TAU_T           equal %f\n" % tau_t
    ret+= "variable        TAU_P           equal %f\n" % tau_p
    ret+= "\n"
    ret+= "units           metal\n"
    ret+= "boundary        p p p\n"
    ret+= "atom_style      atomic\n"
    ret+= "\n"
    ret+= "neighbor        1.0 bin\n"
    if neidelay is not None :
        ret+= "neigh_modify    delay %d\n" % neidelay
    ret+= "\n"
    ret+= "box          tilt large\n"
    ret+= "read_data       %s\n" % conf_file
    ret+= "change_box   all triclinic\n"
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
    ret+= "velocity        all create ${TEMP} %d" % (random.randrange(max_seed-1)+1)
    ret+= "\n"
    if ensemble.split('-')[0] == 'npt' :
        assert (pres is not None)
    if ensemble == "npt" or ensemble == "npt-i" or ensemble == "npt-iso" :
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n"
    elif ensemble == 'npt-a' or ensemble == 'npt-aniso' : 
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} aniso ${PRES} ${PRES} ${TAU_P}\n"
    elif ensemble == 'npt-t' or ensemble == 'npt-tri' : 
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} tri ${PRES} ${PRES} ${TAU_P}\n"
    elif ensemble == "nvt" :
        ret+= "fix             1 all nvt temp ${TEMP} ${TEMP} ${TAU_T}\n"
    else :
        raise RuntimeError("unknown emsemble " + ensemble)
    ret+= "\n"
    ret+= "timestep        %f\n" % dt
    ret+= "run             ${NSTEPS}\n"
    return ret
        
# ret = make_lammps_input ("npt", "al.lmp", ['graph.000.pb', 'graph.001.pb'], 20000, 20, [27], 1000, pres = 1.0)
# print (ret)
# cvt_lammps_conf('POSCAR', 'tmp.lmp')


    
 

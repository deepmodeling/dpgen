#!/usr/bin/env python3

import random, os, sys, dpdata
import numpy as np
import subprocess as sp
import scipy.constants as pc
from distutils.version import LooseVersion

def _sample_sphere() :
    while True:
        vv = np.array([np.random.normal(), np.random.normal(), np.random.normal()])
        vn = np.linalg.norm(vv)
        if vn < 0.2 :
            continue
        return vv / vn

def make_lammps_input(ensemble, 
                      conf_file,
                      graphs,
                      nsteps,
                      dt,
                      neidelay,
                      trj_freq, 
                      mass_map,
                      temp, 
                      jdata,
                      tau_t = 0.1,
                      pres = None,
                      tau_p = 0.5,
                      pka_e = None,
                      ele_temp_f = None,
                      ele_temp_a = None,
                      max_seed = 1000000,
                      nopbc = False,
                      deepmd_version = '0.1') :
    if (ele_temp_f is not None or ele_temp_a is not None) and LooseVersion(deepmd_version) < LooseVersion('1'):
        raise RuntimeError('the electron temperature is only supported by deepmd-kit >= 1.0.0, please upgrade your deepmd-kit')
    if ele_temp_f is not None and ele_temp_a is not None:
        raise RuntimeError('the frame style ele_temp and atom style ele_temp should not be set at the same time')
    ret = "variable        NSTEPS          equal %d\n" % nsteps
    ret+= "variable        THERMO_FREQ     equal %d\n" % trj_freq
    ret+= "variable        DUMP_FREQ       equal %d\n" % trj_freq
    ret+= "variable        TEMP            equal %f\n" % temp
    if ele_temp_f is not None:
        ret+= "variable        ELE_TEMP        equal %f\n" % ele_temp_f
    if ele_temp_a is not None:
        ret+= "variable        ELE_TEMP        equal %f\n" % ele_temp_a
    ret+= "variable        PRES            equal %f\n" % pres
    ret+= "variable        TAU_T           equal %f\n" % tau_t
    ret+= "variable        TAU_P           equal %f\n" % tau_p
    ret+= "\n"
    ret+= "units           metal\n"
    if nopbc:
        ret+= "boundary        f f f\n"
    else:
        ret+= "boundary        p p p\n"
    ret+= "atom_style      atomic\n"
    ret+= "\n"
    ret+= "neighbor        1.0 bin\n"
    if neidelay is not None :
        ret+= "neigh_modify    delay %d\n" % neidelay
    ret+= "\n"
    ret+= "box          tilt large\n"
    ret+= "if \"${restart} > 0\" then \"read_restart dpgen.restart.*\" else \"read_data %s\"\n" % conf_file
    ret+= "change_box   all triclinic\n"
    for jj in range(len(mass_map)) :
        ret+= "mass            %d %f\n" %(jj+1, mass_map[jj])
    graph_list = ""
    for ii in graphs :
        graph_list += ii + " "
    if LooseVersion(deepmd_version) < LooseVersion('1'):
        # 0.x
        ret+= "pair_style      deepmd %s ${THERMO_FREQ} model_devi.out\n" % graph_list
    else:
        # 1.x
        keywords = ""
        if jdata.get('use_clusters', False):
            keywords += "atomic "
        if jdata.get('use_relative', False):
            keywords += "relative %s " % jdata['epsilon']
        if jdata.get('use_relative_v', False):
            keywords += "relative_v %s " % jdata['epsilon_v']
        if ele_temp_f is not None:
            keywords += "fparam ${ELE_TEMP}"
        if ele_temp_a is not None:
            keywords += "aparam ${ELE_TEMP}"
        ret+= "pair_style      deepmd %s out_freq ${THERMO_FREQ} out_file model_devi.out %s\n" % (graph_list, keywords)
    ret+= "pair_coeff      * *\n"
    ret+= "\n"
    ret+= "thermo_style    custom step temp pe ke etotal press vol lx ly lz xy xz yz\n"
    ret+= "thermo          ${THERMO_FREQ}\n"
    model_devi_merge_traj = jdata.get('model_devi_merge_traj', False)
    if(model_devi_merge_traj is True):
        ret+= "dump            1 all custom ${DUMP_FREQ} all.lammpstrj id type x y z fx fy fz\n"
    else:
        ret+= "dump            1 all custom ${DUMP_FREQ} traj/*.lammpstrj id type x y z fx fy fz\n"
    ret+= "restart         10000 dpgen.restart\n"
    ret+= "\n"
    if pka_e is None :
        ret+= "if \"${restart} == 0\" then \"velocity        all create ${TEMP} %d\"" % (random.randrange(max_seed-1)+1)
    else :
        sys = dpdata.System(conf_file, fmt = 'lammps/lmp')
        sys_data = sys.data
        pka_mass = mass_map[sys_data['atom_types'][0] - 1]
        pka_vn = pka_e * pc.electron_volt / \
                 (0.5 * pka_mass * 1e-3 / pc.Avogadro * (pc.angstrom / pc.pico) ** 2)
        pka_vn = np.sqrt(pka_vn)
        print(pka_vn)
        pka_vec = _sample_sphere()
        pka_vec *= pka_vn
        ret+= 'group           first id 1\n'
        ret+= 'if \"${restart} == 0\" then \"velocity        first set %f %f %f\"\n' % (pka_vec[0], pka_vec[1], pka_vec[2])
        ret+= 'fix	       2 all momentum 1 linear 1 1 1\n'
    ret+= "\n"
    if ensemble.split('-')[0] == 'npt' :
        assert (pres is not None)
        if nopbc:
            raise RuntimeError('ensemble %s is conflicting with nopbc' % ensemble)
    if ensemble == "npt" or ensemble == "npt-i" or ensemble == "npt-iso" :
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n"
    elif ensemble == 'npt-a' or ensemble == 'npt-aniso' : 
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} aniso ${PRES} ${PRES} ${TAU_P}\n"
    elif ensemble == 'npt-t' or ensemble == 'npt-tri' : 
        ret+= "fix             1 all npt temp ${TEMP} ${TEMP} ${TAU_T} tri ${PRES} ${PRES} ${TAU_P}\n"
    elif ensemble == "nvt" :
        ret+= "fix             1 all nvt temp ${TEMP} ${TEMP} ${TAU_T}\n"
    elif ensemble == 'nve' :
        ret+= "fix             1 all nve\n"
    else :
        raise RuntimeError("unknown emsemble " + ensemble)
    if nopbc:
        ret+= "velocity        all zero linear\n"
        ret+= "fix             fm all momentum 1 linear 1 1 1\n"
    ret+= "\n"
    ret+= "timestep        %f\n" % dt
    ret+= "run             ${NSTEPS} upto\n"
    return ret
        
# ret = make_lammps_input ("npt", "al.lmp", ['graph.000.pb', 'graph.001.pb'], 20000, 20, [27], 1000, pres = 1.0)
# print (ret)
# cvt_lammps_conf('POSCAR', 'tmp.lmp')


def get_dumped_forces(
        file_name):
    with open(file_name) as fp:        
        lines = fp.read().split('\n')
    natoms = None
    for idx,ii in enumerate(lines):
        if 'ITEM: NUMBER OF ATOMS' in ii:
            natoms = int(lines[idx+1])
            break
    if natoms is None:
        raise RuntimeError('wrong dump file format, cannot find number of atoms', file_name)
    idfx = None
    for idx,ii in enumerate(lines):
        if 'ITEM: ATOMS' in ii:
            keys = ii
            keys = keys.replace('ITEM: ATOMS', '')
            keys = keys.split()
            idfx = keys.index('fx')
            idfy = keys.index('fy')
            idfz = keys.index('fz')
            break
    if idfx is None:
        raise RuntimeError('wrong dump file format, cannot find dump keys', file_name)
    ret = []
    for ii in range(idx+1, idx+natoms+1):
        words = lines[ii].split()
        ret.append([ float(words[ii]) for ii in [idfx, idfy, idfz] ])
    ret = np.array(ret)
    return ret

def get_all_dumped_forces(
        file_name):
    with open(file_name) as fp:        
        lines = fp.read().split('\n')

    ret = []
    exist_natoms = False
    exist_atoms = False

    for idx,ii in enumerate(lines):

        if 'ITEM: NUMBER OF ATOMS' in ii:
            natoms = int(lines[idx+1])
            exist_natoms = True

        if 'ITEM: ATOMS' in ii:
            keys = ii
            keys = keys.replace('ITEM: ATOMS', '')
            keys = keys.split()
            idfx = keys.index('fx')
            idfy = keys.index('fy')
            idfz = keys.index('fz')
            exist_atoms = True
        
            single_traj = []
            for jj in range(idx+1, idx+natoms+1):
                words = lines[jj].split()
                single_traj.append([ float(words[jj]) for jj in [idfx, idfy, idfz] ])
            single_traj = np.array(single_traj)
            ret.append(single_traj)
    
    if exist_natoms is False:
        raise RuntimeError('wrong dump file format, cannot find number of atoms', file_name)
    if exist_atoms is False:
        raise RuntimeError('wrong dump file format, cannot find dump keys', file_name)
    return ret

if __name__ == '__main__':
    ret = get_dumped_forces('40.lammpstrj')
    print(ret)

#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import lib.vasp as vasp
import lib.lammps as lammps
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import yaml
import phonopy
#from phonolammps import Phonolammps



global_equi_name = '00.equi'
global_task_name = '06.phonon'

'''
link poscar
link potcar
make incar
'''
def cmpt_vasp(jdata, conf_dir,opt) :
    fp_params = jdata['vasp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']
    supercell_matrix=jdata['supercell_matrix']
    
    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    cwd = os.getcwd()
    os.chdir(task_path)
    if os.path.isfile('vasprun.xml'):
        os.system('phonopy --fc vasprun.xml')
        os.system('phonopy --dim="%d %d %d" -c POSCAR-unitcell band.conf'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
    else:
        print('vasprun.xml No such file')
    if opt=='Y':
        ph = phonopy.load(supercell_matrix=supercell_matrix,primitive_matrix='auto',unitcell_filename="POSCAR-unitcell",force_constants_filename='FORCE_CONSTANTS')
        ph.auto_band_structure(plot=True).show()
    
    
def cmpt_deepmd_lammps(jdata, conf_dir,opt) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    supercell_matrix=jdata['supercell_matrix']

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'deepmd')

    os.chdir(task_path)
    if os.path.isfile('FORCE_CONSTANTS'):
        os.system('phonopy --dim="%d %d %d" -c POSCAR-unitcell band.conf'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
    else:
        print('FORCE_CONSTANTS No such file')
    if opt=='Y':
        ph = phonopy.load(supercell_matrix=supercell_matrix,primitive_matrix='auto',unitcell_filename="POSCAR-unitcell",force_constants_filename='FORCE_CONSTANTS')
        ph.auto_band_structure(plot=True).show()


def _main() :
    parser = argparse.ArgumentParser(
        description="cmpt 06.phonon")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('OPT', type=str,
                        help='show the band structue or not [Y/N]')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        cmpt_vasp(jdata, args.CONF,args.OPT)               
    elif args.TASK == 'deepmd' :
        cmpt_deepmd_lammps(jdata, args.CONF,args.OPT)
    elif args.TASK == 'meam' :
        cmpt_meam_lammps(jdata, args.CONF)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()


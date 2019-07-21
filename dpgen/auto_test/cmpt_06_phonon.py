#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

import yaml
import phonopy




global_equi_name = '00.equi'
global_task_name = '06.phonon'

'''
link poscar
link potcar
make incar
'''

def get_force_from_dump(cell):
    na=len(cell)
    with open("dump.relax","r") as fp:
        lines = fp.readlines()
        flag=False
        index=[]
        forces=[]
        for line in lines:
            if flag:
                data=line.split()
                index.append(data[0])
                forces.append(data[5:9])
            if len(forces)==na:
                break
            if 'fx fy fz' in line:
                flag=True
        index = np.asarray(index, int)
        indexing = np.argsort(index)
        if len(forces)==na:
            forces=np.asarray(np.reshape(forces,(na,3)),float)[indexing, :]
        else:
            raise RuntimeError('Incomplete result: dump.relax')
    return forces

def cmpt_vasp(jdata, conf_dir) :
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
        os.system('phonopy-bandplot --gnuplot band.yaml > band.dat')
    else:
        print('vasprun.xml No such file')
    
    
def cmpt_deepmd_lammps(jdata, conf_dir) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    supercell_matrix=jdata['supercell_matrix']

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'deepmd')
    task_poscar = os.path.join(task_path, 'POSCAR')

    os.chdir(task_path)
    if os.path.isfile('FORCE_CONSTANTS'):
        os.system('phonopy --dim="%d %d %d" -c POSCAR band.conf'%(supercell_matrix[0],supercell_matrix[1],supercell_matrix[2]))
        os.system('phonopy-bandplot --gnuplot band.yaml > band.dat')
    else:
        print('FORCE_CONSTANTS No such file')


def _main() :
    parser = argparse.ArgumentParser(
        description="cmpt 06.phonon")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        cmpt_vasp(jdata, args.CONF)               
    elif args.TASK == 'deepmd' :
        cmpt_deepmd_lammps(jdata, args.CONF)
    elif args.TASK == 'meam' :
        cmpt_meam_lammps(jdata, args.CONF)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()


#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob, sys
import subprocess as sp
import numpy as np
import lib.vasp as vasp
import lib.lammps as lammps
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.analysis.elasticity.stress import Stress

global_equi_name = '00.equi'
global_task_name = '02.elastic'

def cmpt_vasp(jdata, conf_dir) :
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    equi_stress = Stress(np.loadtxt(os.path.join(task_path, 'equi.stress.out')))

    lst_dfm_path = glob.glob(os.path.join(task_path, 'dfm-*'))
    lst_strain = []
    lst_stress = []
    for ii in lst_dfm_path :
        strain = np.loadtxt(os.path.join(ii, 'strain.out'))
        stress = vasp.get_stress(os.path.join(ii, 'OUTCAR'))
        # convert from pressure in kB to stress
        stress *= -1000
        lst_strain.append(Strain(strain))
        lst_stress.append(Stress(stress))
    et = ElasticTensor.from_independent_strains(lst_strain, lst_stress, eq_stress = equi_stress, vasp = False)
    # et = ElasticTensor.from_independent_strains(lst_strain, lst_stress, eq_stress = None)
    # bar to GPa
    # et = -et / 1e4 
    for ii in range(6) :
        for jj in range(6) :
            sys.stdout.write ("%7.2f " % (et.voigt[ii][jj] / 1e4))
        sys.stdout.write('\n')

def cmpt_deepmd_lammps(jdata, conf_dir) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'lmp')
    equi_stress = Stress(np.loadtxt(os.path.join(task_path, 'equi.stress.out')))

    lst_dfm_path = glob.glob(os.path.join(task_path, 'dfm-*'))
    lst_strain = []
    lst_stress = []
    for ii in lst_dfm_path :
        strain = np.loadtxt(os.path.join(ii, 'strain.out'))
        stress = lammps.get_stress(os.path.join(ii, 'log.lammps'))
        # convert from pressure to stress
        stress = -stress
        lst_strain.append(Strain(strain))
        lst_stress.append(Stress(stress))
    et = ElasticTensor.from_independent_strains(lst_strain, lst_stress, eq_stress = equi_stress, vasp = False)
    # et = ElasticTensor.from_independent_strains(lst_strain, lst_stress, eq_stress = None)
    # bar to GPa
    # et = -et / 1e4 
    for ii in range(6) :
        for jj in range(6) :
            sys.stdout.write ("%7.2f " % (et.voigt[ii][jj] / 1e4))
        sys.stdout.write('\n')

def _main() :
    parser = argparse.ArgumentParser(
        description="gen 01.eos")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

    print('# generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        cmpt_vasp(jdata, args.CONF)               
    elif args.TASK == 'lammps' :
        cmpt_deepmd_lammps(jdata, args.CONF)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

    

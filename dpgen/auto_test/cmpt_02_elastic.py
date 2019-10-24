#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob, sys
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.analysis.elasticity.stress import Stress

global_equi_name = '00.equi'
global_task_name = '02.elastic'

def result_et(et,conf_dir,task_path):
    with open(os.path.join(task_path,'result'),'w') as fp:
        fp.write('conf_dir:%s\n'% (conf_dir))
        for ii in range(6) :
            for jj in range(6) :
                fp.write ("%7.2f " % (et.voigt[ii][jj] / 1e4))
            fp.write('\n')
        BV = et.k_voigt / 1e4
        GV = et.g_voigt / 1e4
        EV = 9*BV*GV/(3*BV+GV)
        uV = 0.5*(3*BV-2*GV)/(3*BV+GV)
        fp.write("# Bulk   Modulus BV = %.2f GPa\n" % (BV))
        fp.write("# Shear  Modulus GV = %.2f GPa\n" % (GV))
        fp.write("# Youngs Modulus EV = %.2f GPa\n" % (EV))
        fp.write("# Poission Ratio uV = %.2f \n" % (uV))

def print_et (et): 
    for ii in range(6) :
        for jj in range(6) :
            sys.stdout.write ("%7.2f " % (et.voigt[ii][jj] / 1e4))
        sys.stdout.write('\n')
    BV = et.k_voigt / 1e4
    GV = et.g_voigt / 1e4
    EV = 9*BV*GV/(3*BV+GV)
    uV = 0.5*(3*BV-2*GV)/(3*BV+GV)
    print("# Bulk   Modulus BV = %.2f GPa" % (BV))
    print("# Shear  Modulus GV = %.2f GPa" % (GV))
    print("# Youngs Modulus EV = %.2f GPa" % (EV))
    print("# Poission Ratio uV = %.2f " % (uV))

def cmpt_vasp(jdata, conf_dir) :
    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        fp_params = jdata['vasp_params']
        kspacing = fp_params['kspacing']
        vasp_str='vasp-k%.2f' % kspacing 
    task_path = os.path.join(task_path, vasp_str)

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
    print_et(et)
    result_et(et,conf_dir,task_path)

def cmpt_deepmd_lammps(jdata, conf_dir, task_name) :
    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_name)
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
    print_et(et)
    result_et(et,conf_dir,task_path)

def _main() :
    parser = argparse.ArgumentParser(
        description="cmpt 02.elastic")
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
    elif args.TASK == 'deepmd' :
        cmpt_deepmd_lammps(jdata, args.CONF, args.TASK)
    elif args.TASK == 'meam' :
        cmpt_deepmd_lammps(jdata, args.CONF, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

    

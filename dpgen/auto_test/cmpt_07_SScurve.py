#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob, sys
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps


global_equi_name = '00.equi'
global_task_name = '07.SScurve'

def cmpt_vasp(jdata, conf_dir) :
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']
    strain_direct=jdata['strain_direct']
    a,b=strain_direct[0],strain_direct[1]    

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)

    lst_dfm_path = glob.glob(os.path.join(task_path, 'dfm-*'))
    lst_strain = []
    lst_stress = []
    print('conf_dir:',conf_dir)
    print('strain\t stress(kB)')
    for ii in lst_dfm_path :
        strain = np.loadtxt(os.path.join(ii, 'strain.out'))
        stress = vasp.get_stress(os.path.join(ii, 'OUTCAR'))
        # convert from pressure in kB to stress
        stress = -stress
        lst_strain.append(strain[a,b])
        lst_stress.append(stress[a,b])
    index=np.argsort(lst_strain)
    for ii in range(len(lst_strain)):
        print('%7.4f %7.4f'%(lst_strain[index[ii]],lst_stress[index[ii]]))
        

def cmpt_lammps(jdata, conf_dir, task_name) :
    strain_direct=jdata['strain_direct']
    a,b=strain_direct[0],strain_direct[1]

    conf_path = os.path.abspath(conf_dir)
    conf_poscar = os.path.join(conf_path, 'POSCAR')
    task_path = re.sub('confs', global_task_name, conf_path)
    task_path = os.path.join(task_path, task_name)
    

    lst_dfm_path = glob.glob(os.path.join(task_path, 'dfm-*'))
    lst_strain = []
    lst_stress = []
    print('conf_dir:',conf_dir)
    print('strain\t stress(kB)')
    for ii in lst_dfm_path :
        strain = np.loadtxt(os.path.join(ii, 'strain.out'))
        stress = lammps.get_stress(os.path.join(ii, 'log.lammps'))
        # convert from pressure to stress
        stress = -stress/1000
        lst_strain.append(strain[a,b])
        lst_stress.append(stress[a,b])
    index=np.argsort(lst_strain)
    for ii in range(len(lst_strain)):
        print('%7.4f %7.4f'%(lst_strain[index[ii]],lst_stress[index[ii]]))

def _main() :
    parser = argparse.ArgumentParser(
        description="cmpt 07.SScurve")
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
        cmpt_lammps(jdata, args.CONF, args.TASK)
    elif args.TASK == 'meam' :
        cmpt_lammps(jdata, args.CONF, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

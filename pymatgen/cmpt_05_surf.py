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
global_task_name = '05.surf'

def cmpt_vasp(jdata, conf_dir, supercell) :
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']

    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    equi_path = os.path.abspath(equi_path)
    equi_outcar = os.path.join(equi_path, 'OUTCAR')
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    task_path = os.path.abspath(task_path)

    equi_natoms, equi_epa, equi_vpa = vasp.get_nev(equi_outcar)

    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_path_widecard = os.path.join(task_path, 'struct-%s-*' % (copy_str))
    struct_path_list = glob.glob(struct_path_widecard)
    struct_path_list.sort()
    if len(struct_path_list) == 0:
        print("# cannot find results for conf %s supercell %s" % (conf_dir, supercell))
    for ii in struct_path_list :
        structure_dir = os.path.basename(ii)
        outcar = os.path.join(ii, 'OUTCAR')
        natoms, epa, vpa = vasp.get_nev(outcar)
        evac = epa * natoms - equi_epa * natoms
        sys.stdout.write ("%s: %7.3f \n" % (structure_dir, evac))    

def cmpt_deepmd_lammps(jdata, conf_dir) :
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, 'lmp')
    equi_path = os.path.abspath(equi_path)
    equi_log = os.path.join(equi_path, 'log.lammps')
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.join(task_path, 'lmp')
    task_path = os.path.abspath(task_path)

    equi_natoms, equi_epa, equi_vpa = lammps.get_nev(equi_log)

    struct_path_widecard = os.path.join(task_path, 'struct-m*m')
    struct_path_list = glob.glob(struct_path_widecard)
    struct_path_list.sort()
    if len(struct_path_list) == 0:
        print("# cannot find results for conf %s" % (conf_dir))
    for ii in struct_path_list :
        structure_dir = os.path.basename(ii)
        lmp_log = os.path.join(ii, 'log.lammps')
        natoms, epa, vpa = lammps.get_nev(lmp_log)
        AA = lammps.get_base_area(lmp_log)
        Cf = 1.60217657e-16 / (1e-20 * 2) * 0.001
        evac = (epa * natoms - equi_epa * natoms) / AA * Cf 
        sys.stdout.write ("%s: \t%7.3f  \n" % (structure_dir, evac))

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

    

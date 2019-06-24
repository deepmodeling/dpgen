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
global_task_name = '03.vacancy'

def comput_e_shift(poscar, task_name) :
    a_types = vasp.get_poscar_types(poscar)
    a_natoms = vasp.get_poscar_natoms(poscar)
    ener_shift = 0
    if not os.path.isdir('stables') :
        raise RuntimeError('no dir "stable". Stable energy and volume of components should be computed before calculating formation energy of an alloy')
    for ii in range(len(a_types)) :
        ref_e_file = a_types[ii] + '.' + task_name + '.e'
        ref_e_file = os.path.join('stables', ref_e_file)
        ener = float(open(ref_e_file).read())
        ener_shift += a_natoms[ii] * ener
    return ener_shift 

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
    print("# ", task_path)

    equi_natoms, equi_epa, equi_vpa = vasp.get_nev(equi_outcar)

    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_path_widecard = os.path.join(task_path, 'struct-%s-*' % (copy_str))
    struct_path_list = glob.glob(struct_path_widecard)
    struct_path_list.sort()
    if len(struct_path_list) == 0:
        print("# cannot find results for conf %s supercell %s" % (conf_dir, supercell))
    sys.stdout.write ("Structure: \tVac_E(eV)  E(eV) equi_E(eV)\n")
    for ii in struct_path_list :
        struct_poscar = os.path.join(ii, 'POSCAR')
        energy_shift = comput_e_shift(struct_poscar, 'vasp-k%.2f' % kspacing)
        structure_dir = os.path.basename(ii)
        outcar = os.path.join(ii, 'OUTCAR')
        natoms, epa, vpa = vasp.get_nev(outcar)
        evac = epa * natoms - equi_epa * natoms
        sys.stdout.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa*natoms))
        # evac = epa * natoms - energy_shift
        # sys.stdout.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, energy_shift))
        # sys.stdout.write ("%s: %7.3f \n" % (structure_dir, evac))    

def cmpt_deepmd_lammps(jdata, conf_dir, supercell, task_name) :
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, task_name)
    equi_path = os.path.abspath(equi_path)
    equi_log = os.path.join(equi_path, 'log.lammps')
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.join(task_path, task_name)
    task_path = os.path.abspath(task_path)
    print("# ", task_path)

    equi_natoms, equi_epa, equi_vpa = lammps.get_nev(equi_log)

    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_path_widecard = os.path.join(task_path, 'struct-%s-*' % (copy_str))
    struct_path_list = glob.glob(struct_path_widecard)
    struct_path_list.sort()
    if len(struct_path_list) == 0:
        print("# cannot find results for conf %s supercell %s" % (conf_dir, supercell))
    sys.stdout.write ("Structure: \tVac_E(eV)  E(eV) equi_E(eV)\n")
    for ii in struct_path_list :
        struct_poscar = os.path.join(ii, 'POSCAR')
        energy_shift = comput_e_shift(struct_poscar, task_name)
        structure_dir = os.path.basename(ii)
        lmp_log = os.path.join(ii, 'log.lammps')
        natoms, epa, vpa = lammps.get_nev(lmp_log)
        evac = epa * natoms - equi_epa * natoms
        sys.stdout.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa * natoms))
        # evac = epa * natoms - energy_shift
        # sys.stdout.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, energy_shift))
        # sys.stdout.write ("%s: %7.3f\n" % (structure_dir, evac))

def _main() :
    parser = argparse.ArgumentParser(
        description="cmpt 03.vacancy")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('COPY', type=int, nargs = 3,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('# generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        cmpt_vasp(jdata, args.CONF, args.COPY)               
    elif args.TASK == 'deepmd' :
        cmpt_deepmd_lammps(jdata, args.CONF, args.COPY, args.TASK)
    elif args.TASK == 'meam' :
        cmpt_deepmd_lammps(jdata, args.CONF, args.COPY, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()

    

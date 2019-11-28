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
global_task_name = '04.interstitial'

def cmpt_vasp(jdata, conf_dir, supercell, insert_ele) :
    for ii in insert_ele:
        _cmpt_vasp(jdata, conf_dir, supercell, ii)

def _cmpt_vasp(jdata, conf_dir, supercell, insert_ele) :

    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing

    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, vasp_str)
    equi_path = os.path.abspath(equi_path)
    equi_outcar = os.path.join(equi_path, 'OUTCAR')
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.join(task_path, vasp_str)
    task_path = os.path.abspath(task_path)

    equi_natoms, equi_epa, equi_vpa = vasp.get_nev(equi_outcar)

    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_path_widecard = os.path.join(task_path, 'struct-%s-%s-*' % (insert_ele, copy_str))
    print(struct_path_widecard)
    struct_path_list = glob.glob(struct_path_widecard)
    struct_path_list.sort()
    if len(struct_path_list) == 0:
        print("# cannot find results for conf %s supercell %s" % (conf_dir, supercell))
    sys.stdout.write ("Insert_ele-Struct: Inter_E(eV)  E(eV) equi_E(eV)\n")
    result = os.path.join(task_path,'result')
    with open(result,'w') as fp:
        fp.write('conf_dir:%s\n'% (conf_dir))
        fp.write ("Insert_ele-Struct: Inter_E(eV)  E(eV) equi_E(eV)\n")
        for ii in struct_path_list :
            structure_dir = os.path.basename(ii)
            outcar = os.path.join(ii, 'OUTCAR')
            natoms, epa, vpa = vasp.get_nev(outcar)
            evac = epa * natoms - equi_epa * natoms
            sys.stdout.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa * natoms))
            fp.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa * natoms))
    fp.close()

def cmpt_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name) :
    for ii in insert_ele:
        _cmpt_deepmd_reprod_traj(jdata, conf_dir, supercell, ii, task_name)

def _cmpt_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name) :
    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    vasp_path = os.path.join(task_path, vasp_str)
    lmps_path = os.path.join(task_path, task_name + vasp_str.replace('vasp','-reprod'))
    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_widecard = os.path.join(vasp_path, 'struct-%s-%s-*' % (insert_ele,copy_str))
    vasp_struct = glob.glob(struct_widecard)
    vasp_struct.sort()
    cwd=os.getcwd()

    for vs in vasp_struct :
        # compute vasp
        outcar = os.path.join(vs, 'OUTCAR')
        vasp_ener = np.array(vasp.get_energies(outcar))
        vasp_ener_file = os.path.join(vs, 'ener.vasp.out')
        # compute reprod
        struct_basename  = os.path.basename(vs)
        ls = os.path.join(lmps_path, struct_basename)
        frame_widecard = os.path.join(ls, 'frame.*')
        frames = glob.glob(frame_widecard)
        frames.sort()
        lmp_ener = []
        for ii in frames :
            log_lmp = os.path.join(ii, 'log.lammps')
            natoms, epa, vpa = lammps.get_nev(log_lmp)
            lmp_ener.append(epa)
        lmp_ener = np.array(lmp_ener)
        lmp_ener = np.reshape(lmp_ener,  [-1,1])
        lmp_ener_file = os.path.join(ls, 'ener.lmp.out')
        vasp_ener = np.reshape(vasp_ener, [-1,1]) / natoms
        error_start = 1
        lmp_ener -= lmp_ener[-1] - vasp_ener[-1]
        diff = lmp_ener - vasp_ener
        diff = diff[error_start:]
        error = np.linalg.norm(diff) / np.sqrt(np.size(lmp_ener))
        np.savetxt(vasp_ener_file, vasp_ener[error_start:])
        np.savetxt(lmp_ener_file,  lmp_ener[error_start:])
        print(os.path.basename(ls), 'EpA_std_err=',error)

def cmpt_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_name) :
    for ii in insert_ele:
        _cmpt_deepmd_lammps(jdata, conf_dir, supercell, ii, task_name)

def _cmpt_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_name) :
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    equi_path = os.path.join(equi_path, task_name.split('-')[0])
    equi_path = os.path.abspath(equi_path)
    equi_log = os.path.join(equi_path, 'log.lammps')
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.join(task_path, task_name)
    task_path = os.path.abspath(task_path)

    equi_natoms, equi_epa, equi_vpa = lammps.get_nev(equi_log)

    copy_str = "%sx%sx%s" % (supercell[0], supercell[1], supercell[2])
    struct_path_widecard = os.path.join(task_path, 'struct-%s-%s-*' % (insert_ele,copy_str))
    struct_path_list = glob.glob(struct_path_widecard)
    print(struct_path_widecard)
    struct_path_list.sort()
    if len(struct_path_list) == 0:
        print("# cannot find results for conf %s supercell %s" % (conf_dir, supercell))
    sys.stdout.write ("Insert_ele-Struct: Inter_E(eV)  E(eV) equi_E(eV)\n")
    result = os.path.join(task_path,'result')
    with open(result,'w') as fp:
        fp.write('conf_dir:%s\n'% (conf_dir))
        fp.write ("Insert_ele-Struct: Inter_E(eV)  E(eV) equi_E(eV)\n")
        for ii in struct_path_list :
            structure_dir = os.path.basename(ii)
            lmp_log = os.path.join(ii, 'log.lammps')
            natoms, epa, vpa = lammps.get_nev(lmp_log)
            evac = epa * natoms - equi_epa * natoms
            sys.stdout.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa * natoms))
            fp.write ("%s: %7.3f  %7.3f %7.3f \n" % (structure_dir, evac, epa * natoms, equi_epa * natoms))
    fp.close()

def _main() :
    parser = argparse.ArgumentParser(
        description="cmpt 04.interstitial")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('COPY', type=int, nargs = 3,
                        help='define the supercell')
    parser.add_argument('ELEMENT', type=str, nargs = '+',
                        help='the inserted element')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

#    print('# generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        cmpt_vasp(jdata, args.CONF, args.COPY, args.ELEMENT)
    elif args.TASK == 'deepmd' :
        cmpt_deepmd_lammps(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'deepmd-reprod' :
        cmpt_deepmd_reprod_traj(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'meam' :
        cmpt_deepmd_lammps(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    elif args.TASK == 'meam-reprod' :
        cmpt_deepmd_reprod_traj(jdata, args.CONF, args.COPY, args.ELEMENT, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)

if __name__ == '__main__' :
    _main()

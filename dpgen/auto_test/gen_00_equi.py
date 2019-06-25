#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps

global_task_name = '00.equi'

'''
link poscar
link potcar
make incar
'''
def make_vasp(jdata, conf_dir) :
    fp_params = jdata['vasp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    kgamma = fp_params['kgamma']

    conf_path = os.path.abspath(conf_dir)
    equi_path = re.sub('confs', global_task_name, conf_path)
    os.makedirs(equi_path, exist_ok = True)
    cwd = os.getcwd()
    from_poscar = os.path.join(conf_path, 'POSCAR')
    to_poscar = os.path.join(equi_path, 'POSCAR')
    if os.path.exists(to_poscar) :
        assert(filecmp.cmp(from_poscar, to_poscar))
    else :
        os.chdir(equi_path)
        os.symlink(os.path.relpath(from_poscar), 'POSCAR')
        os.chdir(cwd)
    is_alloy = \
               os.path.exists(
                   os.path.join(
                       os.path.join(conf_path, '..'),
                       'alloy'
                   )
               )    
    vasp_path = os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    os.makedirs(vasp_path, exist_ok = True)
    os.chdir(vasp_path)
    print(vasp_path)
    # gen incar
    if is_alloy :
        fc = vasp.make_vasp_relax_incar(ecut, ediff, True,  True, True, npar, kpar, kspacing, kgamma)
    else :
        fc = vasp.make_vasp_relax_incar(ecut, ediff, False, True, True, npar, kpar, kspacing, kgamma)
    with open('INCAR', 'w') as fp :
        fp.write(fc)
    # gen poscar
    if os.path.exists('POSCAR') :
        os.remove('POSCAR')
    os.symlink(os.path.relpath(to_poscar), 'POSCAR')
    # gen potcar
    with open('POSCAR','r') as fp :
        lines = fp.read().split('\n')
        ele_list = lines[5].split()
    potcar_map = jdata['potcar_map']
    potcar_list = []
    for ii in ele_list :
        assert(os.path.exists(potcar_map[ii]))
        potcar_list.append(potcar_map[ii])
    with open('POTCAR', 'w') as outfile:
        for fname in potcar_list:
            with open(fname) as infile:
                outfile.write(infile.read())
    os.chdir(cwd)

def make_deepmd_lammps (jdata, conf_dir) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))    
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]

    conf_path = os.path.abspath(conf_dir)
    equi_path = re.sub('confs', global_task_name, conf_path)
    os.makedirs(equi_path, exist_ok = True)
    cwd = os.getcwd()
    from_poscar = os.path.join(conf_path, 'POSCAR')
    to_poscar = os.path.join(equi_path, 'POSCAR')
    if os.path.exists(to_poscar) :
        assert(filecmp.cmp(from_poscar, to_poscar))
    else :
        os.chdir(equi_path)
        os.symlink(os.path.relpath(from_poscar), 'POSCAR')
        os.chdir(cwd)
    # lmp path
    lmp_path = os.path.join(equi_path, 'deepmd')
    os.makedirs(lmp_path, exist_ok = True)    
    print(lmp_path)
    # lmp conf
    conf_file = os.path.join(lmp_path, 'conf.lmp')
    lammps.cvt_lammps_conf(to_poscar, os.path.relpath(conf_file))
    ptypes = vasp.get_poscar_types(to_poscar)
    lammps.apply_type_map(conf_file, deepmd_type_map, ptypes)    
    # lmp input
    fc = lammps.make_lammps_equi(os.path.basename(conf_file), 
                                 ntypes, 
                                 lammps.inter_deepmd, 
                                 deepmd_models_name)
    with open(os.path.join(lmp_path, 'lammps.in'), 'w') as fp :
        fp.write(fc)
    # link models
    os.chdir(lmp_path)
    for ii in deepmd_models_name :
        if os.path.exists(ii) :
            os.remove(ii)
    for (ii,jj) in zip(deepmd_models, deepmd_models_name) :
        os.symlink(os.path.relpath(ii), jj)
    os.chdir(cwd)

def make_meam_lammps (jdata, conf_dir) :
    meam_potfile_dir = jdata['meam_potfile_dir']
    meam_potfile_dir = os.path.abspath(meam_potfile_dir)
    meam_potfile = jdata['meam_potfile']
    meam_potfile = [os.path.join(meam_potfile_dir,ii) for ii in meam_potfile]
    meam_potfile_name = jdata['meam_potfile']
    type_map = jdata['meam_type_map']
    ntypes = len(type_map)
    meam_param = {'meam_potfile' :      jdata['meam_potfile'],
                  'meam_type':          jdata['meam_param_type']}

    conf_path = os.path.abspath(conf_dir)
    equi_path = re.sub('confs', global_task_name, conf_path)
    os.makedirs(equi_path, exist_ok = True)
    cwd = os.getcwd()
    from_poscar = os.path.join(conf_path, 'POSCAR')
    to_poscar = os.path.join(equi_path, 'POSCAR')
    if os.path.exists(to_poscar) :
        assert(filecmp.cmp(from_poscar, to_poscar))
    else :
        os.chdir(equi_path)
        os.symlink(os.path.relpath(from_poscar), 'POSCAR')
        os.chdir(cwd)
    # lmp path
    lmp_path = os.path.join(equi_path, 'meam')
    os.makedirs(lmp_path, exist_ok = True)    
    print(lmp_path)
    # lmp conf
    conf_file = os.path.join(lmp_path, 'conf.lmp')
    lammps.cvt_lammps_conf(to_poscar, os.path.relpath(conf_file))
    ptypes = vasp.get_poscar_types(to_poscar)
    lammps.apply_type_map(conf_file, type_map, ptypes)    
    # lmp input
    fc = lammps.make_lammps_equi(os.path.basename(conf_file), 
                                 ntypes, 
                                 lammps.inter_meam, 
                                 meam_param)
    with open(os.path.join(lmp_path, 'lammps.in'), 'w') as fp :
        fp.write(fc)
    # link models
    os.chdir(lmp_path)
    for ii in meam_potfile_name :
        if os.path.exists(ii) :
            os.remove(ii)
    for (ii,jj) in zip(meam_potfile, meam_potfile_name) :
        os.symlink(os.path.relpath(ii), jj)
    os.chdir(cwd)

def _main() :
    parser = argparse.ArgumentParser(
        description="gen 00.equi")
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
        make_vasp(jdata, args.CONF)               
    elif args.TASK == 'deepmd' :
        make_deepmd_lammps(jdata, args.CONF)
    elif args.TASK == 'meam' :
        make_meam_lammps(jdata, args.CONF)
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()


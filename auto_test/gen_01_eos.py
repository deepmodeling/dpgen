#!/usr/bin/env python3

import os, re, argparse, filecmp, json, glob
import subprocess as sp
import numpy as np
import lib.vasp as vasp
import lib.lammps as lammps
from pymatgen.core.structure import Structure

global_equi_name = '00.equi'
global_task_name = '01.eos'

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
    vol_start = jdata['vol_start']
    vol_end = jdata['vol_end']
    vol_step = jdata['vol_step']

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', global_task_name, conf_path)
    os.makedirs(task_path, exist_ok = True)
    cwd = os.getcwd()
    from_poscar = os.path.join(conf_path, 'POSCAR')
    to_poscar = os.path.join(task_path, 'POSCAR')
    if os.path.exists(to_poscar) :
        assert(filecmp.cmp(from_poscar, to_poscar))
    else :
        os.chdir(task_path)
        os.symlink(os.path.relpath(from_poscar), 'POSCAR')
        os.chdir(cwd)
    vol_to_poscar = vasp.poscar_vol(to_poscar) / vasp.poscar_natoms(to_poscar)
    # print(to_poscar, vol_to_poscar)
    is_alloy = \
               os.path.exists(
                   os.path.join(
                       os.path.join(conf_path, '..'),
                       'alloy'
                   )
               )    
    vasp_path = os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    os.makedirs(vasp_path, exist_ok = True)
    os.chdir(vasp_path)
    # gen incar
    if is_alloy :
        fc = vasp.make_vasp_relax_incar(ecut, ediff, True,  True, False, npar, kpar, kspacing, kgamma)
    else :
        fc = vasp.make_vasp_relax_incar(ecut, ediff, False, True, False, npar, kpar, kspacing, kgamma)
    with open('INCAR', 'w') as fp :
        fp.write(fc)
    # gen potcar
    with open(to_poscar,'r') as fp :
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
    # loop over volumes
    for vol in np.arange(vol_start, vol_end, vol_step) :
        vol_path = os.path.join(vasp_path, 'vol-%.2f' % vol)        
        os.makedirs(vol_path, exist_ok = True)
        os.chdir(vol_path)
        print(vol_path)
        for ii in ['INCAR', 'POTCAR', 'POSCAR.orig', 'POSCAR'] :
            if os.path.exists(ii) :
                os.remove(ii)
        # link incar, potcar
        os.symlink(os.path.relpath(os.path.join(vasp_path, 'INCAR')), 'INCAR')
        os.symlink(os.path.relpath(os.path.join(vasp_path, 'POTCAR')), 'POTCAR')
        # gen poscar
        os.symlink(os.path.relpath(to_poscar), 'POSCAR.orig')
        scale = (vol / vol_to_poscar) ** (1./3.)
        # print(scale)
        vasp.poscar_scale('POSCAR.orig', 'POSCAR', scale)
        # print(vol_path, vasp.poscar_vol('POSCAR') / vasp.poscar_natoms('POSCAR'))
        os.chdir(cwd)

def make_deepmd_lammps (jdata, conf_dir) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]
    vol_start = jdata['vol_start']
    vol_end = jdata['vol_end']
    vol_step = jdata['vol_step']

    # # get equi props
    # equi_path = re.sub('confs', global_equi_name, conf_path)
    # equi_path = os.path.join(equi_path, 'lmp')
    # equi_log = os.path.join(equi_path, 'log.lammps')
    # if not os.path.isfile(equi_log) :
    #     raise RuntimeError("the system should be equilibriated first")
    # natoms, epa, vpa = lammps.get_nev(equi_log)
    # task path
    task_path = re.sub('confs', global_task_name, conf_dir)
    task_path = os.path.abspath(task_path)
    os.makedirs(task_path, exist_ok = True)
    cwd = os.getcwd()
    conf_path = os.path.abspath(conf_dir)
    from_poscar = os.path.join(conf_path, 'POSCAR')
    to_poscar = os.path.join(task_path, 'POSCAR')
    if os.path.exists(to_poscar) :
        assert(filecmp.cmp(from_poscar, to_poscar))
    else :
        os.chdir(task_path)
        os.symlink(os.path.relpath(from_poscar), 'POSCAR')
        os.chdir(cwd)
    volume = vasp.poscar_vol(to_poscar)
    natoms = vasp.poscar_natoms(to_poscar)
    vpa = volume / natoms
    # structrure
    ss = Structure.from_file(to_poscar)
    # lmp path
    lmp_path = os.path.join(task_path, 'deepmd')
    os.makedirs(lmp_path, exist_ok = True)
    # # lmp conf
    # conf_file = os.path.join(lmp_path, 'conf.lmp')
    # lammps.cvt_lammps_conf(to_poscar, conf_file)
    # ptypes = vasp.get_poscar_types(to_poscar)
    # lammps.apply_type_map(conf_file, deepmd_type_map, ptypes)
    for vol in np.arange(vol_start, vol_end, vol_step) :
        vol_path = os.path.join(lmp_path, 'vol-%.2f' % vol)        
        print('# generate %s' % (vol_path))
        os.makedirs(vol_path, exist_ok = True)
        os.chdir(vol_path)
        #print(vol_path)
        for ii in ['conf.lmp', 'conf.lmp'] + deepmd_models_name :
            if os.path.exists(ii) :
                os.remove(ii)                
        # # link conf
        # os.symlink(os.path.relpath(conf_file), 'conf.lmp')
        # make conf
        scale_ss = ss.copy()
        scale_ss.scale_lattice(vol * natoms)
        scale_ss.to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', deepmd_type_map, ptypes)
        # link models
        for (ii,jj) in zip(deepmd_models, deepmd_models_name) :
            os.symlink(os.path.relpath(ii), jj)
        # make lammps input
        scale = (vol / vpa) ** (1./3.)
        fc = lammps.make_lammps_press_relax('conf.lmp', ntypes, scale,lammps.inter_deepmd, deepmd_models_name)
        with open(os.path.join(vol_path, 'lammps.in'), 'w') as fp :
            fp.write(fc)
        os.chdir(cwd)

def make_deepmd_lammps_fixv (jdata, conf_dir) :
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_type_map = jdata['deepmd_type_map']
    ntypes = len(deepmd_type_map)    
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    deepmd_models = glob.glob(os.path.join(deepmd_model_dir, '*pb'))
    deepmd_models_name = [os.path.basename(ii) for ii in deepmd_models]
    vol_start = jdata['vol_start']
    vol_end = jdata['vol_end']
    vol_step = jdata['vol_step']

    # get equi props
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    task_path = re.sub('confs', global_task_name, conf_dir)
    equi_path = os.path.join(equi_path, 'deepmd')
    task_path = os.path.join(task_path, 'deepmd')
    equi_path = os.path.abspath(equi_path)
    task_path = os.path.abspath(task_path)
    equi_log = os.path.join(equi_path, 'log.lammps')
    equi_dump = os.path.join(equi_path, 'dump.relax')
    os.makedirs(task_path, exist_ok = True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    lammps.poscar_from_last_dump(equi_dump, task_poscar, deepmd_type_map)

    cwd = os.getcwd()
    volume = vasp.poscar_vol(task_poscar)
    natoms = vasp.poscar_natoms(task_poscar)
    vpa = volume / natoms
    # structrure
    ss = Structure.from_file(task_poscar)
    # make lammps.in
    fc = lammps.make_lammps_equi('conf.lmp', 
                                 ntypes, 
                                 lammps.inter_deepmd, 
                                 deepmd_models_name, 
                                 change_box = False)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    # make vols
    for vol in np.arange(vol_start, vol_end, vol_step) :
        vol_path = os.path.join(task_path, 'vol-%.2f' % vol)        
        print('# generate %s' % (vol_path))
        os.makedirs(vol_path, exist_ok = True)
        os.chdir(vol_path)
        for ii in ['conf.lmp', 'conf.lmp', 'lammps.in'] + deepmd_models_name :
            if os.path.exists(ii) :
                os.remove(ii)                
        # make conf
        scale_ss = ss.copy()
        scale_ss.scale_lattice(vol * natoms)
        scale_ss.to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', deepmd_type_map, ptypes)
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(deepmd_models, deepmd_models_name) :
            os.symlink(os.path.relpath(ii), jj)
        # make lammps input
        os.chdir(cwd)

def make_meam_lammps_fixv (jdata, conf_dir) :
    meam_potfile_dir = jdata['meam_potfile_dir']
    meam_potfile_dir = os.path.abspath(meam_potfile_dir)
    meam_potfile = jdata['meam_potfile']
    meam_potfile = [os.path.join(meam_potfile_dir,ii) for ii in meam_potfile]
    meam_potfile_name = jdata['meam_potfile']
    type_map = jdata['meam_type_map']
    ntypes = len(type_map)
    meam_param = {'meam_potfile' :      jdata['meam_potfile'],
                  'meam_type':          jdata['meam_param_type']}

    vol_start = jdata['vol_start']
    vol_end = jdata['vol_end']
    vol_step = jdata['vol_step']

    # get equi props
    equi_path = re.sub('confs', global_equi_name, conf_dir)
    task_path = re.sub('confs', global_task_name, conf_dir)
    equi_path = os.path.join(equi_path, 'meam')
    task_path = os.path.join(task_path, 'meam')
    equi_path = os.path.abspath(equi_path)
    task_path = os.path.abspath(task_path)
    equi_log = os.path.join(equi_path, 'log.lammps')
    equi_dump = os.path.join(equi_path, 'dump.relax')
    os.makedirs(task_path, exist_ok = True)
    task_poscar = os.path.join(task_path, 'POSCAR')
    lammps.poscar_from_last_dump(equi_dump, task_poscar, type_map)

    cwd = os.getcwd()
    volume = vasp.poscar_vol(task_poscar)
    natoms = vasp.poscar_natoms(task_poscar)
    vpa = volume / natoms
    # structrure
    ss = Structure.from_file(task_poscar)
    # make lammps.in
    fc = lammps.make_lammps_equi('conf.lmp', 
                                 ntypes, 
                                 lammps.inter_meam, 
                                 meam_param, 
                                 change_box = False)
    f_lammps_in = os.path.join(task_path, 'lammps.in')
    with open(f_lammps_in, 'w') as fp :
        fp.write(fc)
    # make vols
    for vol in np.arange(vol_start, vol_end, vol_step) :
        vol_path = os.path.join(task_path, 'vol-%.2f' % vol)        
        print('# generate %s' % (vol_path))
        os.makedirs(vol_path, exist_ok = True)
        os.chdir(vol_path)
        for ii in ['conf.lmp', 'conf.lmp', 'lammps.in'] + meam_potfile_name :
            if os.path.exists(ii) :
                os.remove(ii)                
        # make conf
        scale_ss = ss.copy()
        scale_ss.scale_lattice(vol * natoms)
        scale_ss.to('POSCAR', 'POSCAR')
        lammps.cvt_lammps_conf('POSCAR', 'conf.lmp')
        ptypes = vasp.get_poscar_types('POSCAR')
        lammps.apply_type_map('conf.lmp', type_map, ptypes)
        # link lammps.in
        os.symlink(os.path.relpath(f_lammps_in), 'lammps.in')
        # link models
        for (ii,jj) in zip(meam_potfile, meam_potfile_name) :
            os.symlink(os.path.relpath(ii), jj)
        # make lammps input
        os.chdir(cwd)

def _main() :
    parser = argparse.ArgumentParser(
        description="gen 01.eos")
    parser.add_argument('TASK', type=str,
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    parser.add_argument('-f', '--fix-shape', action = 'store_true',
                        help='fix shape of box')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

    # print('generate %s task with conf %s' % (args.TASK, args.CONF))
    if args.TASK == 'vasp':
        make_vasp(jdata, args.CONF)               
    elif args.TASK == 'deepmd' :
        if args.fix_shape :
            make_deepmd_lammps_fixv(jdata, args.CONF)
        else :
            make_deepmd_lammps(jdata, args.CONF)        
    elif args.TASK == 'meam' :
        if args.fix_shape:
            make_meam_lammps_fixv(jdata, args.CONF)
        else :
            raise RuntimeError("not implemented ", args.TASK)            
    else :
        raise RuntimeError("unknow task ", args.TASK)
    
if __name__ == '__main__' :
    _main()


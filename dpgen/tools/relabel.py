#!/usr/bin/env python3

import os,sys,json,glob,argparse,shutil
import numpy as np
import subprocess as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from dpgen.generator.lib.pwscf import make_pwscf_input
from dpgen.generator.lib.siesta import make_siesta_input
from dpgen.generator.run import make_vasp_incar, update_mass_map
import dpdata

def get_lmp_info(input_file) :
    with open(input_file) as fp:
        lines = [line.rstrip('\n') for line in fp]
    for ii in lines :
        words = ii.split()
        if len(words) >= 4 and words[0] == 'variable' :
            if words[1] == 'TEMP' :
                temp = float(words[3])
            elif words[1] == 'PRES' :
                pres = float(words[3])
        if len(words) >= 4 and words[0] == 'fix' :
            if words[3] == 'nvt' :
                ens = 'nvt'
            elif words[3] == 'npt' :
                ens = 'npt'
    return ens, temp, pres


def link_pp_files(tdir, fp_pp_path, fp_pp_files) :
    cwd = os.getcwd()
    os.chdir(tdir)
    for ii in fp_pp_files :
        if os.path.lexists(ii) :
            os.remove(ii)
        os.symlink(os.path.join(fp_pp_path, ii), ii)    
    os.chdir(cwd)


def copy_pp_files(tdir, fp_pp_path, fp_pp_files) :
    cwd = os.getcwd()
    os.chdir(tdir)
    for ii in fp_pp_files :
        if os.path.lexists(ii) :
            os.remove(ii)
        if os.path.exists(ii) :
            os.remove(ii)
        shutil.copyfile(os.path.join(fp_pp_path, ii), ii)    
    os.chdir(cwd)
    
        
def make_vasp(tdir, fp_params) :
    cwd = os.getcwd()
    os.chdir(tdir)
    incar = make_vasp_incar(fp_params)
    with open('INCAR', 'w') as fp:
        fp.write(incar)
    os.chdir(cwd)
        
def make_vasp_incar(tdir, fp_incar) :
    cwd = os.getcwd()
    os.chdir(tdir)
    shutil.copyfile(fp_incar, 'INCAR')
    os.chdir(cwd)        

def make_pwscf(tdir, fp_params, mass_map, fp_pp_path, fp_pp_files, user_input) : 
    cwd = os.getcwd()
    os.chdir(tdir)
    sys_data = dpdata.System('POSCAR').data
    sys_data['atom_masses'] = mass_map
    ret = make_pwscf_input(sys_data, fp_pp_files, fp_params)
    open('input', 'w').write(ret)        
    os.chdir(cwd)

def make_siesta(tdir, fp_params, fp_pp_path, fp_pp_files) :
    cwd = os.getcwd()
    os.chdir(tdir)
    sys_data = dpdata.System('POSCAR').data
    ret = make_siesta_input(sys_data, fp_pp_files, fp_params)
    open('input', 'w').write(ret)
    os.chdir(cwd)

def create_init_tasks(target_folder, param_file, output, fp_json, verbose = True) :
    target_folder = os.path.abspath(target_folder)
    output = os.path.abspath(output)
    tool_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'template')    
    jdata = json.load(open(os.path.join(target_folder, param_file)))
    update_mass_map(jdata)
    fp_jdata = json.load(open(fp_json))
    # fp settings
    mass_map = jdata['mass_map']
    type_map = jdata['type_map']
    fp_style = fp_jdata['fp_style']
    fp_pp_path = fp_jdata['fp_pp_path']
    fp_pp_files = fp_jdata['fp_pp_files']
    cwd_ = os.getcwd()
    os.chdir(target_folder)
    fp_pp_path = os.path.abspath(fp_pp_path)
    os.chdir(cwd_)
    # init data sys
    init_data_prefix = jdata['init_data_prefix']
    init_data_sys = jdata['init_data_sys']
    for idx,ii in enumerate(init_data_sys):        
        sys = dpdata.LabeledSystem(os.path.join(init_data_prefix, ii), fmt = 'deepmd/npy', type_map = type_map)
        nframes = sys.get_nframes()
        sys_dir = os.path.join(output, 'init_system.%03d' % idx)
        os.makedirs(sys_dir, exist_ok = True)
        if verbose :
            print('# working on ' + sys_dir)
        with open(os.path.join(sys_dir,'record'), 'w') as fp:
            fp.write(os.path.join(init_data_prefix, ii) + '\n')
        for ff in range(nframes) :
            task_dir = os.path.join(sys_dir, 'task.%06d' % ff)
            os.makedirs(task_dir, exist_ok = True)
            sys.to_vasp_poscar(os.path.join(task_dir, 'POSCAR'), frame_idx=ff)
            # make fp
            cwd_ = os.getcwd()
            os.chdir(task_dir)
            for pp in fp_pp_files :
                if os.path.lexists(pp) :
                    os.remove(pp)
                os.symlink(os.path.relpath(os.path.join(output, pp)), pp)
            if fp_style == 'vasp':
                if os.path.lexists('INCAR') :
                    os.remove('INCAR')
                os.symlink(os.path.relpath(os.path.join(output, 'INCAR')), 'INCAR')
            elif fp_style == 'pwscf':
                try:                    
                    fp_params = fp_jdata['user_fp_params']
                    user_input = True
                except Exception:
                    fp_params = fp_jdata['fp_params']
                    user_input = False
                make_pwscf('.', fp_params, mass_map, fp_pp_files, fp_pp_files, user_input)
            elif fp_style == 'siesta':
                make_siesta('.', fp_params, fp_pp_files, fp_pp_files)
            os.chdir(cwd_)            
    

def create_tasks(target_folder, param_file, output, fp_json, verbose = True, numb_iter = -1) :
    target_folder = os.path.abspath(target_folder)
    output = os.path.abspath(output)
    tool_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'template')    
    jdata = json.load(open(os.path.join(target_folder, param_file)))
    update_mass_map(jdata)
    fp_jdata = json.load(open(fp_json))
    # goto input        
    cwd = os.getcwd()
    os.chdir(target_folder)
    sys = jdata['sys_configs']    
    # fp settings
    mass_map = jdata['mass_map']
    fp_style = fp_jdata['fp_style']
    fp_pp_path = fp_jdata['fp_pp_path']
    fp_pp_files = fp_jdata['fp_pp_files']
    cwd_ = os.getcwd()
    os.chdir(target_folder)
    fp_pp_path = os.path.abspath(fp_pp_path)
    os.chdir(cwd_)
    # fp_params = fp_jdata['fp_params']
    # collect tasks from iter dirs
    sys_tasks = [[] for ii in sys]
    sys_tasks_record = [[] for ii in sys]
    sys_tasks_cc = [0 for ii in sys]
    numb_sys = len(sys)
    iters = glob.glob('iter.[0-9]*[0-9]')
    iters.sort()
    # iters = iters[:2]
    for ii in iters[:numb_iter] :
        iter_tasks = glob.glob(os.path.join(ii, '02.fp', 'task.[0-9]*[0-9].[0-9]*[0-9]'))
        iter_tasks.sort()
        if verbose :
            print('# check iter ' + ii + ' with %6d tasks' % len(iter_tasks))
        for jj in iter_tasks :
            sys_idx = int(os.path.basename(jj).split('.')[-2])
            sys_tasks[sys_idx].append(jj)
            if os.path.islink(os.path.join(jj, 'conf.lmp')):
                linked_file = os.path.realpath(os.path.join(jj, 'conf.lmp'))
            elif os.path.islink(os.path.join(jj, 'conf.dump')):
                linked_file = os.path.realpath(os.path.join(jj, 'conf.dump'))
            else:
                raise RuntimeError('cannot file linked conf file')
            linked_keys = linked_file.split('/')
            task_record = linked_keys[-5] + '.' + linked_keys[-3] + '.' + linked_keys[-1].split('.')[0]
            task_record_keys = task_record.split('.')
            ens, temp, pres = get_lmp_info(os.path.join(ii, '01.model_devi', linked_keys[-3], 'input.lammps'))
            human_record = 'iter: %s   system: %s   model_devi_task: %s   frame: %6d   fp_task: %s   ens: %s   temp: %10.2f   pres: %10.2f' \
                           % (task_record_keys[1],
                              task_record_keys[3],
                              linked_keys[-3],
                              int(task_record_keys[-1]),
                              os.path.basename(jj),
                              ens, temp, pres
                           )
            # print(human_record)
            sys_tasks_record[sys_idx].append(human_record)            
    # for ii in range(numb_sys) :
    #     for jj in range(len(sys_tasks[ii])) :
    #         print(sys_tasks[ii][jj], sys_tasks_record[ii][jj])

    # mk output
    os.makedirs(output, exist_ok = True)
    if fp_style == 'vasp':
        copy_pp_files(output, fp_pp_path, fp_pp_files)
        make_vasp_incar(fp_params, output)
    if fp_style == 'pwscf' :
        copy_pp_files(output, fp_pp_path, fp_pp_files)
    if fp_style == 'siesta' :
        copy_pp_files(output, fp_pp_path, fp_pp_files)
    for si in range(numb_sys) :
        sys_dir = os.path.join(output, 'system.%03d' % si)
        if verbose :
            print('# working on ' + sys_dir)
        for tt,rr in zip(sys_tasks[si], sys_tasks_record[si]) :            
            # copy poscar
            source_path = os.path.join(('iter.%s/02.fp' % rr.split()[1]), rr.split()[9])
            source_file = os.path.join(source_path, 'POSCAR')
            target_path = os.path.join(sys_dir, 'task.%06d'%sys_tasks_cc[si])
            sys_tasks_cc[si] += 1
            os.makedirs(target_path, exist_ok = True)
            target_file = os.path.join(target_path, 'POSCAR')
            target_recd = os.path.join(target_path, 'record')
            if os.path.exists(target_file) :
                os.remove(target_file)
            if os.path.exists(target_recd) :
                os.remove(target_recd)
            shutil.copyfile(source_file, target_file)
            with open(target_recd, 'w') as fp:
                fp.write('\n'.join([target_folder, rr, '']))
            # make fp
            cwd_ = os.getcwd()
            os.chdir(target_path)
            for pp in fp_pp_files :
                if os.path.lexists(pp) :
                    os.remove(pp)
                os.symlink(os.path.relpath(os.path.join(output, pp)), pp)
            if fp_style == 'vasp':
                if os.path.lexists('INCAR') :
                    os.remove('INCAR')
                os.symlink(os.path.relpath(os.path.join(output, 'INCAR')), 'INCAR')
            elif fp_style == 'pwscf':
                try:                    
                    fp_params = fp_jdata['user_fp_params']
                    user_input = True
                except Exception:
                    fp_params = fp_jdata['fp_params']
                    user_input = False
                make_pwscf('.', fp_params, mass_map, fp_pp_files, fp_pp_files, user_input)
            elif fp_style == 'siesta':
                make_siesta('.', fp_params, mass_map, fp_pp_files, fp_pp_files)
            os.chdir(cwd_)
    os.chdir(cwd)


def _main()   :
    parser = argparse.ArgumentParser(description='Create tasks for relabeling from a DP-GEN job')
    parser.add_argument("JOB_DIR", type=str, 
                        help="the directory of the DP-GEN job")
    parser.add_argument("PARAM", type=str, default = 'fp.json',
                        help="the json file defines vasp tasks")
    parser.add_argument("OUTPUT", type=str, 
                        help="the output directory of relabel tasks")
    parser.add_argument('-p',"--parameter", type=str, default = 'param.json',
                        help="the json file provides DP-GEN paramters, should be located in JOB_DIR")
    parser.add_argument('-n',"--numb-iter", type=int, default = -1,
                        help="number of iterations to relabel")
    parser.add_argument('-v',"--verbose", action = 'store_true',
                        help="being loud")
    args = parser.parse_args()
            
    create_tasks(args.JOB_DIR, args.parameter, args.OUTPUT, args.PARAM, numb_iter = args.numb_iter, verbose = args.verbose)
    create_init_tasks(args.JOB_DIR, args.parameter, args.OUTPUT, args.PARAM, verbose = args.verbose)

if __name__ == '__main__':
    _main()
    

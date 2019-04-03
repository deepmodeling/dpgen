#!/usr/bin/env python3

import os,sys,json,glob,argparse,shutil
import numpy as np
import subprocess as sp
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from lib.vasp import make_vasp_incar
from lib.vasp import system_from_poscar


def get_lmp_info(input_file) :
    lines = [line.rstrip('\n') for line in open(input_file)]
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

        
def make_vasp(tdir, fp_params) :
    cwd = os.getcwd()
    os.chdir(tdir)
    incar = make_vasp_incar(fp_params)
    with open('INCAR', 'w') as fp:
        fp.write(incar)
    os.chdir(cwd)
        

def make_pwscf(tdir, fp_params, mass_map, fp_pp_path, fp_pp_files) :        
    cwd = os.getcwd()
    os.chdir(tdir)
    sys_data = system_from_poscar('POSCAR')
    sys_data['atom_masses'] = jdata['mass_map']
    ret = make_pwscf_input(sys_data, fp_pp_files, fp_params)
    open('input', 'w').write(ret)        
    os.chdir(cwd)
    

def create_tasks(target_folder, param_file, output, fp_json, verbose = True) :
    target_folder = os.path.abspath(target_folder)
    output = os.path.abspath(output)
    tool_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'template')    
    jdata = json.load(open(os.path.join(target_folder, param_file)))
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
    fp_params = fp_jdata['fp_params']
    # collect tasks from iter dirs
    sys_tasks = [[] for ii in sys]
    sys_tasks_record = [[] for ii in sys]
    sys_tasks_cc = [0 for ii in sys]
    numb_sys = len(sys)
    iters = glob.glob('iter.[0-9]*[0-9]')
    iters.sort()
    # iters = iters[:2]
    for ii in iters :
        iter_tasks = glob.glob(os.path.join(ii, '02.fp', 'task.[0-9]*[0-9].[0-9]*[0-9]'))
        iter_tasks.sort()
        if verbose :
            print('# check iter ' + ii + ' with %6d tasks' % len(iter_tasks))
        for jj in iter_tasks :
            sys_idx = int(os.path.basename(jj).split('.')[-2])
            sys_tasks[sys_idx].append(jj)
            linked_file = os.path.realpath(os.path.join(jj, 'conf.lmp'))
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
        link_pp_files(output, fp_pp_path, fp_pp_files)
        make_vasp(output, fp_params)
    for si in range(numb_sys) :
        sys_dir = os.path.join(output, 'system.%03d' % si)
        for tt,rr in zip(sys_tasks[si], sys_tasks_record[si]) :            
            # copy poscar
            source_path = os.path.join(('iter.%s/02.fp' % rr.split()[1]), rr.split()[9])
            source_file = os.path.join(source_path, 'POSCAR')
            target_path = os.path.join(sys_dir, '%06d'%sys_tasks_cc[si])
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
                make_pwscf('.', fp_params, mass_map, fp_pp_files, fp_pp_files)
            os.chdir(cwd_)


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
    parser.add_argument('-v',"--verbose", action = 'store_true',
                        help="being loud")
    args = parser.parse_args()
            
    create_tasks(args.JOB_DIR, args.parameter, args.OUTPUT, args.PARAM)

if __name__ == '__main__':
    _main()
    

#!/usr/bin/env python3

import os,sys,json,glob,argparse,shutil
import numpy as np
import subprocess as sp
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from relabel import get_lmp_info

def ascii_hist(count) :
    np = (count-1) // 5 + 1
    ret = " |"
    for ii in range(np):
        ret += '='
    return ret    

def stat_tasks(target_folder, param_file, verbose = True) :
    target_folder = os.path.abspath(target_folder)
    tool_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'template')    
    jdata = json.load(open(os.path.join(target_folder, param_file)))
    # goto input        
    cwd = os.getcwd()
    os.chdir(target_folder)
    sys = jdata['sys_configs']    
    numb_sys = len(sys)
    sys_tasks_count = [0 for ii in sys]
    sys_tasks_trait = [[] for ii in sys]
    sys_tasks_trait_count = [[] for ii in sys]
    # collect tasks from iter dirs
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
            sys_tasks_count[sys_idx] += 1
            linked_file = os.path.realpath(os.path.join(jj, 'conf.lmp'))
            linked_keys = linked_file.split('/')
            task_record = linked_keys[-5] + '.' + linked_keys[-3] + '.' + linked_keys[-1].split('.')[0]
            task_record_keys = task_record.split('.')
            ens, temp, pres = get_lmp_info(os.path.join(ii, '01.model_devi', linked_keys[-3], 'input.lammps'))
            trait = [ens, temp, pres]
            if not trait in sys_tasks_trait[sys_idx] :
                sys_tasks_trait[sys_idx].append(trait)
                sys_tasks_trait_count[sys_idx].append(0)
            t_idx = sys_tasks_trait[sys_idx].index(trait)
            sys_tasks_trait_count[sys_idx][t_idx] += 1
    sys_tasks_all = []
    for ii in range(numb_sys) :
        # print(sys[ii], sys_tasks_count[ii])
        tmp_all = []
        for jj in range(len(sys_tasks_trait[ii])) :
            tmp_all.append(sys_tasks_trait[ii][jj] + [sys_tasks_trait_count[ii][jj]])
        sys_tasks_all.append(tmp_all)
    for ii in sys_tasks_all:
        ii.sort()
    max_str_len = max([len(str(ii)) for ii in sys])
    sys_fmt = '%%%ds   %%6d' % (max_str_len+1)
    blank = max_str_len - 50
    str_blk = ""
    for ii in range(blank):
        str_blk += " "
    trait_fmt = str_blk + 'ens: %s   T: %10.2f   P: %12.2f   count:   %6d'
    for ii in range(numb_sys):
        print(sys_fmt % (str(sys[ii]), sys_tasks_count[ii]))
        for jj in range(len(sys_tasks_all[ii])):
            hist_str = ascii_hist(sys_tasks_all[ii][jj][3])
            print((trait_fmt + hist_str) % (sys_tasks_all[ii][jj][0],
                                            sys_tasks_all[ii][jj][1],
                                            sys_tasks_all[ii][jj][2],
                                            sys_tasks_all[ii][jj][3]))
            

def _main()   :
    parser = argparse.ArgumentParser(description='Some data statistics of DP-GEN iterations')
    parser.add_argument("JOB_DIR", type=str, 
                        help="the directory of the DP-GEN job")
    parser.add_argument('-p',"--parameter", type=str, default = 'param.json',
                        help="the json file provides DP-GEN paramters, should be located in JOB_DIR")
    parser.add_argument('-v',"--verbose", action = 'store_true',
                        help="being loud")
    args = parser.parse_args()
    
    stat_tasks(args.JOB_DIR, args.parameter, args.verbose)


if __name__ == '__main__':
    _main()

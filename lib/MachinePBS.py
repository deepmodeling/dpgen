#!/usr/bin/env python3

import os
import subprocess as sp
import multiprocessing as mp

# global vars
env_source_list = []
b_vcores = False
encoding_sys = 'utf-8'

def add_source_file (source_file) :    
    if os.path.isfile (source_file) :
        sf = os.path.abspath (source_file)
        global env_source_list
        env_source_list += [sf]
    else :
        raise RuntimeError ("no file " + source_file)

def has_virtual_cores (yes_or_no) :
    global b_vcores
    b_vcores = yes_or_no

def get_node_list () :
    node_file = os.getenv('PBS_NODEFILE')
    with open(node_file, 'r') as fp:
        nodes = fp.read().rstrip().split('\n')
    a = list(set(nodes))
    a.sort()
    # a = ['cu03', 'cu04']
    return a

def get_core_per_node () :
    ncpu = os.cpu_count()
    global b_vcores
    if b_vcores :
        ncpu = ncpu // 2
    return ncpu

def cmd_source () :
    global env_source_list
    cmd = ""
    for ii in env_source_list :
        cmd += "source " + ii + ";"
    return cmd
    
def exec_batch (cmd,
                cmd_dir_,
                task_batch,
                args_batch,
                work_threads) :
    host_list = get_node_list()
    nnode = len(host_list)
    ntasks = len(task_batch)

    cwd = os.getcwd()
    cmd_dir = os.path.abspath(cmd_dir_)
    os.chdir(cmd_dir)
    ret = []
    for ii in range(ntasks) :
        os.chdir(task_batch[ii])
        exec_host = host_list[ii%nnode]
        with open("nodes", "w") as fp:
            fp.write(exec_host)
        sph = sp.Popen("mpirun -n 1 -hostfile nodes %s" % cmd, shell = True)
        ret.append(sph)
        os.chdir(cmd_dir)
    os.chdir(cwd)
    return ret
    
    
# print(get_node_list()    )
# print(get_core_per_node()   )
# exec_batch("pwd", "/home/wang_han/study/deepgen.test/", ['iter.000000/00.train/000', 'iter.000000/00.train/001', 'iter.000000/00.train/002', 'iter.000000/00.train/003'], ["", "", "", ""], 7)

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
    hostname = sp.check_output ("hostname", shell = True).rstrip()
    hostname = hostname.decode(encoding_sys) 
    return [hostname]

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

def exec_cmd (node,
              cmd,
              cmd_dir_,
              cmd_args) :
    cmd_dir = os.path.abspath(cmd_dir_)
    run_cmd = ""
    run_cmd += cmd_source()
    run_cmd += "cd " + cmd_dir + ";"
    run_cmd += cmd + " " + cmd_args
    return sp.Popen(run_cmd, shell = True)

def exec_batch (cmd,
                cmd_dir_,
                task_batch,
                args_batch,
                work_thread) :
    cmd_dir = os.path.abspath(cmd_dir_)
    cwd = os.getcwd()
    os.chdir(cmd_dir)
    host_list = get_node_list()
    nnode = len(host_list)
    ps = []
    for jj in range(len(task_batch)) :
        work_path = task_batch[jj]
        work_args = args_batch[jj]
        host = host_list[jj % nnode]
        work_name = os.path.basename (work_path)
        ps.append(exec_cmd(host, cmd, work_path, work_args))    
    os.chdir(cwd)
    return ps

# # nlist = get_node_list()
# # print (nlist)
# cpnode = get_core_per_node()
# print (cpnode)
# # add_source_file ('a')
# # add_source_file ('b')
# print (env_source_list)
# cmd = cmd_source ()
# print (cmd)
# sp.check_call (cmd, shell = True)

# p = exec_cmd("localhost", "ls", "~/study", "-lt")
# p.wait()

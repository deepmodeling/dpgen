#!/usr/bin/env python3

import os
import subprocess as sp
import multiprocessing as mp

# global vars
env_source_list = []
b_vcores = True
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
    node_dict = os.getenv('SLURM_JOB_NODELIST')
    ret = sp.check_output ("scontrol show hostnames " + node_dict, shell = True).rstrip()
    ret = ret.decode(encoding_sys)
    return ret.split('\n')

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
    ssh_run = 'ssh %s "%s" 2>/dev/null' % (node, run_cmd)
    return sp.Popen(ssh_run, shell = True)

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
    ntpnode = get_core_per_node()
    nsource = ntpnode * nnode    
    numb_jobs = (nsource // work_thread)
    fp = open ("mpmd.conf", "w")
    # print ("in exec batch")
    # print (len(task_batch))
    # print (work_thread)
    # print (nsource)
    iend = numb_jobs
    if len(task_batch) < numb_jobs :
        iend = len(task_batch)
    for ii in range(iend) :
        thread_s = ii 
        thread_e = (ii+1) - 1
        fp.write ("%d-%d %s/run.mpmd.sh %s %s\n" % (thread_s, thread_e, cmd_dir, task_batch[ii], args_batch[ii]))
    fp.close()
    fp = open ("run.mpmd.sh", "w")
    fp.write ("#!/bin/bash\n")
    fp.write (cmd_source())
    fp.write ("\n")
    fp.write ("cd $1\n")
    fp.write ("shift 1\n")
    fp.write ("%s $*\n" % cmd)
    fp.close()
    os.chmod("run.mpmd.sh", 0o755)
    # print ("srun --tasks-per-node=%d --multi-prog ./mpmd.conf"% (ntpnode//work_thread))
    sph = sp.Popen ("srun --tasks-per-node=%d --multi-prog ./mpmd.conf" % (ntpnode//work_thread),
                    shell = True)
    os.chdir(cwd)
    return [sph]
    

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

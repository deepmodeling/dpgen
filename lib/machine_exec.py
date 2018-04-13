#!/usr/bin/env python3

import os
import numpy as np
import subprocess as sp
import logging
import time

def run_node_tasks (max_thread,
                    work_thread,
                    all_task,
                    run_cmd) :
    numb_jobs = max_thread // work_thread
    task_chunks = [all_task[i:i + numb_jobs] for i in range(0, len(all_task), numb_jobs)]
    base_path = os.getcwd() + "/"
    count_batch = 0
    for task_batch in task_chunks :
        ps = []
        for work_path in task_batch :
            work_name = os.path.basename (work_path)
            log_task ("%03d %s: %s" % (count_batch, work_name, run_cmd))
            os.chdir(work_path)    
            ps.append(sp.Popen(run_cmd, shell = True))
            os.chdir(base_path)
        while True :
            if not(any(p.wait() for p in ps)) :
                break
            time.sleep(1)
        count_batch += 1    


def exec_hosts (machine_env,
                cmd, 
                work_thread,
                task_dirs,
                task_args = None) :
    ntasks = len(task_dirs)
    if task_args != None :
        assert ntasks == len(task_args) or len(task_args) == 1
        if len(task_args) == 1 :
            tmp_arg = task_args[0]
            task_args = [tmp_arg for ii in range(ntasks)]
    else :
        task_args = ["" for ii in range(ntasks)]
    assert ntasks == len(task_args)    

    host_list = machine_env.get_node_list()
    nnode = len(host_list)
    ntpnode = machine_env.get_core_per_node()
    nsource = ntpnode * nnode
    numb_jobs = (nsource // work_thread)
    task_chunks = [task_dirs[i:i + numb_jobs] for i in range(0, ntasks, numb_jobs)]
    args_chunks = [task_args[i:i + numb_jobs] for i in range(0, ntasks, numb_jobs)]
    nbatch = len(task_chunks)
    assert nbatch == len(args_chunks)

    base_path = os.getcwd() + "/"
    for ii in range(nbatch) :
        task_batch = task_chunks[ii]
        args_batch = args_chunks[ii]
        ps = []
        for jj in range(len(task_batch)) :
            work_path = task_batch[jj]
            work_args = args_batch[jj]
            host = host_list[jj % nnode]
            work_name = os.path.basename (work_path)
            logging.info(("%s %03d %s: %s %s" % (host, ii, work_name, cmd, work_args)))
            ps.append(machine_env.exec_cmd(host, cmd, work_path, work_args))
        while True :
            if not(any(p.wait() for p in ps)) :
                break
            time.sleep(1)    

def exec_hosts_batch (machine_env,
                      cmd, 
                      work_thread,
                      task_dirs,
                      task_args = None) :
    ntasks = len(task_dirs)
    if task_args != None :
        assert ntasks == len(task_args) or len(task_args) == 1
        if len(task_args) == 1 :
            tmp_arg = task_args[0]
            task_args = [tmp_arg for ii in range(ntasks)]
    else :
        task_args = ["" for ii in range(ntasks)]
    assert ntasks == len(task_args)    

    host_list = machine_env.get_node_list()
    nnode = len(host_list)
    ntpnode = machine_env.get_core_per_node()
    nsource = ntpnode * nnode
    numb_jobs = (nsource // work_thread)
    task_chunks = [task_dirs[i:i + numb_jobs] for i in range(0, ntasks, numb_jobs)]
    args_chunks = [task_args[i:i + numb_jobs] for i in range(0, ntasks, numb_jobs)]
    nbatch = len(task_chunks)
    assert nbatch == len(args_chunks)

    base_path = os.getcwd() + "/"
    for ii in range(nbatch) :
        task_batch = task_chunks[ii]
        args_batch = args_chunks[ii]
        logging.info(("%s %03d : %s with %d jobs %s" % (host_list, ii, cmd, len(task_batch), task_batch)))
        ps = machine_env.exec_batch(cmd, ".", task_batch, task_args, work_thread)
        while True :
            if not(any(p.wait() for p in ps)) :
                break
            time.sleep(1)    

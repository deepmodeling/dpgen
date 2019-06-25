#!/usr/bin/env python3

import os, time
from subprocess import Popen, PIPE
import subprocess as sp

from lib.BatchJob import JobStatus
from lib.BatchJob import BatchJob
from lib.SlurmJob import SlurmJob

def make_slurm_script (cmd,
                       numb_node = 1,
                       work_thread = 1,
                       numb_gpu = 0, 
                       task_args = None, 
                       time_limit = "0:30:0",
                       mem_limit = 32,
                       modules = None,
                       sources = None,
                       fin_tag = 'tag_finished') :
    ret = ""
    ret += "#!/bin/bash -l\n"
    ret += "#SBATCH -N %d\n" % numb_node
    ret += "#SBATCH --exclude tiger-i23g1\n"
    ret += "#SBATCH -t %s\n" % time_limit
    ret += "#SBATCH --mem %dG \n" % mem_limit
    ret += "#SBATCH --ntasks-per-node %d\n" % work_thread
    if numb_gpu > 0 :
        ret += "#SBATCH --gres=gpu:%d\n" % numb_gpu
    ret += "\n"
    for ii in modules :
        ret += "module load %s\n" % ii
    ret += "\n"
    for ii in sources :
        ret += "source %s\n" %ii
    ret += "\n"        
    if task_args is not None :
        ret += cmd + task_args + "\n"
    else :
        ret += cmd + "\n"
    ret += "if test $? -eq 0; then\n"
    ret += "    touch %s\n" % fin_tag
    ret += "fi\n"
    ret += "sleep 1\n"
    return ret

def make_slurm_script_group (cmd,
                             task_dir,
                             numb_node = 1,
                             work_thread = 1,
                             numb_gpu = 0,
                             task_args = None, 
                             time_limit = "0:30:0",
                             mem_limit = 32,
                             modules = None,
                             sources = None,
                             fin_tag = 'tag_finished') :
    if task_args is not None :
        assert(len(task_dir) == len(task_args))

    ret = ""
    ret += "#!/bin/bash -l\n"
    ret += "#SBATCH -N %d\n" % numb_node
    ret += "#SBATCH -t %s\n" % time_limit
    ret += "#SBATCH --mem %dG \n" % mem_limit
    ret += "#SBATCH --ntasks-per-node %d\n" % work_thread
    if numb_gpu > 0 :
        ret += "#SBATCH --gres=gpu:%d\n" % numb_gpu
    ret += "\n"
#    ret += "set -euo pipefail\n"
    for ii in modules :
        ret += "module load %s\n" % ii
    ret += "\n"
    for ii in sources :
        ret += "source %s\n" %ii
    ret += "\n"        
    ret += "cwd=`pwd`"
    ret += "\n"        
    for ii in range(len(task_dir)) :
        ret += "cd " + str(task_dir[ii]) + "\n"
        if task_args is not None :
            ret += cmd + task_args[ii] + "\n"
        else :
            ret += cmd + "\n"
        ret += "if test $? -ne 0; then exit ; fi\n"
        ret += "cd $cwd\n"
        ret += "\n"        
    ret += "if test $? -eq 0; then\n"
    ret += "    touch %s\n" % fin_tag
    ret += "fi\n"
    ret += "sleep 1\n"
    return ret

def exec_batch (cmd,
                numb_node,
                work_thread,
                numb_gpu,
                task_dirs,
                task_args = None, 
                time_limit = "24:0:0", 
                mem_limit = 32,
                modules = None,
                sources = None) :
    cwd = os.getcwd()
    job_list = []
    fin_tag = 'tag_finished'
    for ii,mydir in enumerate(task_dirs) :
        os.chdir(mydir)
        myarg = None
        if task_args is not None :
            myarg = task_args[ii]
        with open('_sub', 'w') as fp :
            fp.write(make_slurm_script(cmd,
                                       numb_node, work_thread, numb_gpu,
                                       myarg,
                                       time_limit,
                                       mem_limit,
                                       modules,
                                       sources,
                                       fin_tag))
        job = SlurmJob(os.getcwd(), '_sub', job_finish_tag = fin_tag)
        job_list.append (job)
        os.chdir(cwd)

    for ii in job_list:
        ii.submit()
#        time.sleep(1)

    while True :
        find_unfinish = False
        for job in job_list :
            stat = job.check_status ()
            if stat == JobStatus.terminated :
                raise RuntimeError("find terminated job")
                old_job_id = job.get_job_id()
                new_job_id = job.submit ()
                find_unfinish = True
            if stat != JobStatus.finished :
                find_unfinish = True
        if find_unfinish == False :
            return
        else :
            time.sleep (10)

def exec_batch_group (cmd,
                      numb_node,
                      work_thread,
                      numb_gpu,
                      task_dirs_,
                      group_size = 10,
                      task_args = None, 
                      time_limit = "24:0:0", 
                      mem_limit = 32,
                      modules = None,
                      sources = None) :    
    cwd = os.getcwd()
    job_list = []
    fin_tag = 'tag_finished'
    
    os.chdir(task_dirs_[0])
    os.chdir('..')
    working_dir = os.getcwd()
    os.chdir(cwd)

    task_dirs = []
    for ii in task_dirs_ :
        task_dirs.append(os.path.abspath(ii))
    if task_args is not None :
        assert(len(task_dirs) == len(task_args))
    if task_args is None :
        task_args = []
        for ii in task_dirs :
            task_args.append("")            

    ntasks = len(task_dirs)
    task_chunks = [task_dirs[i:i + group_size] for i in range(0, ntasks, group_size)]
    args_chunks = [task_args[i:i + group_size] for i in range(0, ntasks, group_size)]

    os.chdir(working_dir)
    for ii in range(len(task_chunks)):
        group_dir = "group.%06d" % ii
        if not os.path.isdir(group_dir) :
            os.mkdir(group_dir)
        os.chdir(group_dir)
        with open('_sub', 'w') as fp:
            fp.write(make_slurm_script_group(cmd,
                                             task_chunks[ii],
                                             numb_node, work_thread, numb_gpu,
                                             args_chunks[ii],
                                             time_limit,
                                             mem_limit,
                                             modules,
                                             sources,
                                             fin_tag))
            job = SlurmJob(os.getcwd(), '_sub', job_finish_tag = fin_tag)
        job_list.append (job)
        os.chdir(working_dir)
    os.chdir(cwd)

    # for ii,mydir in enumerate(task_dirs) :
    #     os.chdir(mydir)
    #     myarg = None
    #     if task_args is not None :
    #         myarg = task_args[ii]
    #     with open('_sub', 'w') as fp :
    #         fp.write(make_slurm_script(cmd, work_thread, numb_gpu, myarg, time_limit, mem_limit, modules, sources, fin_tag))
    #     job = SlurmJob(os.getcwd(), '_sub', job_finish_tag = fin_tag)
    #     job_list.append (job)
    #     os.chdir(cwd)

    for ii in job_list:
        ii.submit()
#        time.sleep(1)

    while True :
        find_unfinish = False
        for job in job_list :
            stat = job.check_status ()
            if stat == JobStatus.terminated :
                raise RuntimeError("find terminated job")
                old_job_id = job.get_job_id()
                new_job_id = job.submit ()
                find_unfinish = True
            if stat != JobStatus.finished :
                find_unfinish = True
        if find_unfinish == False :
            return
        else :
            time.sleep (10)


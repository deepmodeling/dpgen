#!/usr/bin/evn python3

"""
init: data
iter:
        00.train
        01.mode_devi
        02.vasp
        03.data
"""

import os
import sys
import argparse
import glob
import json
import random
import logging
import numpy as np
from lib.utils import make_iter_name
from lib.utils import create_path
from lib.utils import copy_file_list
from lib.utils import replace
from lib.utils import cmd_append_log
from lib.lammps import cvt_lammps_conf
from lib.lammps import make_lammps_input
import lib.MachineLocal as MachineLocal
import lib.MachineSlurm as MachineSlurm
from lib.machine_exec import exec_hosts
from lib.machine_exec import exec_hosts_batch

template_name = 'template'
train_name = '00.train'
train_task_fmt = '%03d'
train_files = ['input.json']
train_param = 'input.json'
train_tmpl_path = os.path.join(template_name, train_name)
model_devi_name = '01.model_devi'
model_devi_task_fmt = '%03d.%06d'
model_devi_conf_fmt = '%03d.%04d'

def get_job_names(jdata) :
    jobkeys = []
    for ii in jdata.keys() :
        if ii.split('_')[0] == "job" :
            jobkeys.append(ii)
    jobkeys.sort()
    return jobkeys

def make_model_devi_task_name (sys_idx, task_idx) :
    return "task." + model_devi_task_fmt % (sys_idx, task_idx)

def make_model_devi_conf_name (sys_idx, conf_idx) :
    return model_devi_conf_fmt % (sys_idx, conf_idx)


def make_train (iter_index, 
               jdata) :    
    # load json param
    numb_models = jdata['numb_models']
    init_data_sys_ = jdata['init_data_sys']    
    init_data_sys = []
    for ii in init_data_sys_ :
        init_data_sys.append(os.path.abspath(ii))
    for ii in init_data_sys :
        if not os.path.isdir(ii) :
            raise RuntimeError ("data sys %s does not exists" % ii)
    # establish work path
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    create_path(work_path)
    copy_file_list(train_files, train_tmpl_path, work_path)
    # establish tasks
    jinput = json.load(open(os.path.join(work_path, train_param), 'r'))
    jinput['systems'] = init_data_sys    
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(task_path)
        jinput['seed'] = random.randrange(sys.maxsize)
        with open(os.path.join(task_path, train_param), 'w') as outfile:
            json.dump(jinput, outfile, indent = 4)

def run_train (iter_index,
               jdata, 
               exec_machine) :    
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = jdata['deepmd_path']
    train_nthreads = jdata['train_nthreads']
    # paths
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    # make tasks
    all_task = []
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        all_task.append(task_path)
    command = os.path.join(deepmd_path, 'bin/dp_train') + ' ' + train_param
    command = cmd_append_log (command, 'train.log')
    # train models
    exec_hosts_batch(exec_machine, command, train_nthreads, all_task, None)

def post_train (iter_index,
                jdata) :
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = jdata['deepmd_path']
    # paths
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    all_task = []
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        all_task.append(task_path)
    command = os.path.join(deepmd_path, 'bin/dp_frz')
    command = cmd_append_log(command, 'freeze.log')
    # frz models
    exec_hosts(MachineLocal, command, 1, all_task, None)
    # symlink models
    for ii in range(numb_models) :
        task_file = os.path.join(train_task_fmt % ii, 'frozen_model.pb')
        ofile = os.path.join(work_path, 'graph.%03d.pb' % ii)
        os.remove(ofile)
        os.symlink(task_file, ofile)

def make_model_devi (iter_index, 
                     jdata, 
                     FPMananger) :
    model_devi_jobs = jdata['model_devi_jobs']
    ensemble = model_devi_jobs['ensemble']
    nsteps = model_devi_jobs['nsteps']
    trj_freq = model_devi_jobs['trj_freq']
    job_names = get_job_names (model_devi_jobs)
    assert (iter_index < len(job_names)) 
    cur_job_name = job_names[iter_index]    
    cur_job = model_devi_jobs[cur_job_name]
    if 'ensemble' in cur_job.keys() :
        ensemble = cur_job['ensemble']
    if 'nsteps' in cur_job.keys() :
        nsteps = cur_job['nsteps']
    if 'trj_freq' in cur_job.keys() :
        trj_freq = cur_job['trj_freq']

    conf_systems_glob = cur_job['systems']
    conf_systems = []
    for ss in conf_systems_glob :
        cur_systems = []
        for ii in ss :
            cur_systems += glob.glob(ii)
        cur_systems.sort()
        conf_systems.append (cur_systems)
    mass_map = cur_job['mass_map']
    temps = cur_job['temps']
    press = cur_job['press']

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = glob.glob(os.path.join(train_path, "graph*pb"))    
    task_model_list = [] 
    for ii in models: 
        task_model_list.append(os.path.join('..', os.path.basename(ii)))
    work_path = os.path.join(iter_name, model_devi_name)
    create_path(work_path)
    for mm in models :
        model_name = os.path.basename(mm)
        os.symlink(mm, os.path.join(work_path, model_name))

    all_task = []
    task_param = []
    conf_path = os.path.join(work_path, 'confs')
    create_path(conf_path)
    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        for cc in ss :            
            conf_name = make_model_devi_conf_name(sys_counter, conf_counter)
            poscar_name = conf_name + '.poscar'
            lmp_name = conf_name + '.lmp'
            os.symlink(cc, os.path.join(conf_path, poscar_name))
            all_task.append(conf_path)
            task_param.append(' ' + poscar_name + ' ' + lmp_name)
            conf_counter += 1
        sys_counter += 1
    exec_hosts(MachineLocal, os.path.join(os.getcwd(), "lib/ovito_file_convert.py"), 1, all_task, task_param)

    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        task_counter = 0
        for cc in ss :            
            for tt in temps:
                for pp in press:
                    task_name = make_model_devi_task_name(sys_counter, task_counter)
                    conf_name = make_model_devi_conf_name(sys_counter, conf_counter) + '.lmp'
                    task_path = os.path.join(work_path, task_name)
                    # print(task_path)
                    create_path(task_path)
                    create_path(os.path.join(task_path, 'traj'))
                    loc_conf_name = 'conf.lmp'
                    os.symlink(os.path.join(os.path.join('..','confs'), conf_name), 
                               os.path.join(task_path, loc_conf_name) )
                    file_c = make_lammps_input(ensemble,
                                               loc_conf_name,
                                               task_model_list,
                                               nsteps,
                                               trj_freq,
                                               mass_map,
                                               tt,
                                               pres = pp)
                    with open(os.path.join(task_path, 'input.lammps'), 'w') as fp :
                        fp.write(file_c)
                    # cvt_lammps_conf(cc, 'conf.lmp')
                    task_counter += 1
            conf_counter += 1
        sys_counter += 1

def run_model_devi (iter_index, 
                    jdata, 
                    exec_machine) :
    model_devi_np = jdata['model_devi_np']
    lmp_path = jdata['lmp_path']
    lmp_exec = os.path.join(lmp_path, 'lmp_mpi')

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = lmp_exec + " -i input.lammps"
    command = cmd_append_log(command, "model_devi.log")
    
    exec_hosts_batch(exec_machine, command, model_devi_np, all_task)

def post_model_devi (iter_index, 
                     jdata) :
    model_devi_trust = jdata['model_devi_trust']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = "awk '{if ($2 > %e) {print $1,$2}}' model_devi.out > sel.out" % model_devi_trust

    exec_hosts (MachineLocal, command, 1, all_task)

    # for ii in all_task :
    #     task_path = ii
    #     model_devi = np.loadtxt(os.path.join(task_path, "model_devi.out"))
    #     all_idx = model_devi[:,0]
    #     # all_idx = all_idx.astype(int)
    #     all_idx = all_idx.astype(int)
    #     sel_idx = []
    #     sel = model_devi[:,1] > model_devi_trust
    #     for kk in range(len(sel)) :
    #         if sel[kk]:
    #             sel_idx.append(all_idx[kk])
    #     print(sel)
    #     print(model_devi)
    #     print(all_idx)
    #     print (sel_idx)
    #     exit(0)
    
def run_iter (json_file, exec_machine) :
    prev_model = init_model
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_iter = jdata["numb_iter"]
    numb_task = 8
    record = "record.dpgen"

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        logging.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    for ii in range (numb_iter) :
        if ii > 0 :
            prev_model = glob.glob (make_iter_name(ii-1) + "/" + train_name + "/*pb")
        for jj in range (numb_task) :
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] : 
                continue
            if   jj == 0 :
                log_iter ("make_train", ii, jj)
                make_train (ii, jdata, prev_model) 
            elif jj == 1 :
                log_iter ("run_train", ii, jj)
                run_train  (ii, jdata, exec_machine)
            # else :
            #     raise RuntimeError ("unknow task %d, something wrong" % jj)

            record_iter (record, ii, jj)



def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("JSON", type=str, 
                        help="The json parameter")
    parser.add_argument("--machine", type=str, 
                        help="The machine settings")        
    args = parser.parse_args()

    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')
    # logging.basicConfig (filename="compute_string.log", filemode="a", level=logging.INFO, format='%(asctime)s %(message)s')

    machine_type = "local"    
    gmxrc = None
    vcores = None
    if args.machine != None :
        fp = open (args.machine, 'r')
        jdata = json.load (fp)
        machine_type = jdata["machine_type"]
        gmxrc = jdata["gmxrc"]
        vcores = jdata["virtual_cores"]

    if   machine_type == "local" :
        exec_machine = MachineLocal
    elif machine_type == "slurm" :
        exec_machine = MachineSlurm

    if vcores != None:
        exec_machine.has_virtual_cores(vcores)
    if gmxrc != None:
        exec_machine.add_source_file(gmxrc)

    logging.info ("start running")
    run_iter (args.JSON, exec_machine)
    logging.info ("finished!")

if __name__ == '__main__':
    # _main()
    fp = open ('param.json', 'r')
    jdata = json.load (fp)
    # post_train(0, jdata)
    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')
    # make_model_devi(0, jdata, None)
    # run_model_devi(0, jdata, MachineLocal)
    post_model_devi(0, jdata)

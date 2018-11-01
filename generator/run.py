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
import warnings
import shutil
import time
import numpy as np
import subprocess as sp
from lib.utils import make_iter_name
from lib.utils import create_path
from lib.utils import copy_file_list
from lib.utils import replace
from lib.utils import cmd_append_log
from lib.utils import log_iter
from lib.utils import record_iter
from lib.utils import log_task
from lib.lammps import cvt_lammps_conf
from lib.lammps import make_lammps_input
from lib.vasp import make_vasp_incar
from lib.vasp import system_from_poscar
from lib.pwscf import make_pwscf_input
import lib.MachineLocal as MachineLocal
import lib.MachineLocalGPU as MachineLocalGPU
import lib.MachineSlurm as MachineSlurm
import lib.MachinePBS as MachinePBS
from lib.machine_exec import exec_hosts
from lib.machine_exec import exec_hosts_batch
from lib.batch_exec import exec_batch
from lib.batch_exec import exec_batch_group
from lib.RemoteJob import SSHSession, JobStatus, SlurmJob, CloudMachineJob

template_name = 'template'
train_name = '00.train'
train_task_fmt = '%03d'
train_tmpl_path = os.path.join(template_name, train_name)
data_system_fmt = '%03d'
model_devi_name = '01.model_devi'
model_devi_task_fmt = data_system_fmt + '.%06d'
model_devi_conf_fmt = data_system_fmt + '.%04d'
fp_name = '02.fp'
fp_task_fmt = data_system_fmt + '.%06d'

import requests
from hashlib import sha1

def _verfy_ac(private_key, params):
    items= sorted(params.items())
    
    params_data = "";
    for key, value in items:
        params_data = params_data + str(key) + str(value)
    params_data = params_data + private_key
    sign = sha1()
    sign.update(params_data.encode())
    signature = sign.hexdigest()
    return signature

def _ucloud_remove_machine(machine, UHostId):
    ucloud_url = machine['url']
    ucloud_stop_param = machine['ucloud_param']
    ucloud_stop_param['Action'] = "StopUHostInstance"
    ucloud_stop_param['UHostId'] = UHostId
    ucloud_stop_param['Signature'] = _verfy_ac(machine['Private'], ucloud_stop_param)
    req = requests.get(ucloud_url, ucloud_stop_param)
    if req.json()['RetCode'] != 0 :
        raise RuntimeError ("failed to stop ucloud machine")

    ucloud_delete_param = machine['ucloud_param']
    ucloud_delete_param['Action'] = "TerminateUHostInstance"
    ucloud_delete_param['UHostId'] = UHostId
    ucloud_delete_param['Signature'] = _verfy_ac(machine['Private'], ucloud_delete_param)
    req = requests.get(ucloud_url, ucloud_delete_param)
    
    if req.json()['RetCode'] != 0 :
        raise RuntimeError ("failed to terminate ucloud machine")

def _ucloud_submit_jobs(machine,
                        command,
                        work_path,
                        tasks,
                        group_size,
                        forward_common_files,
                        forward_task_files,
                        backward_task_files) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    assert machine['machine_type'] == 'ucloud'
    ucloud_url = machine['url']
    ucloud_start_param = machine['ucloud_param']
    ucloud_start_param['Action'] = "CreateUHostInstance"
    ucloud_start_param['Name'] = "train"
    ucloud_start_param['Signature'] = _verfy_ac(machine['Private'], ucloud_start_param)

    njob = len(task_chunks)
    ucloud_machines = []
    ucloud_hostids = []
    for ii in range(njob) :
        req = requests.get(ucloud_url, ucloud_start_param)
        if req.json()['RetCode'] != 0 :
            print(json.dumps(req.json(),indent=2, sort_keys=True))
            raise RuntimeError ("failed to start ucloud machine")
        ucloud_machines.append(str(req.json()["IPs"][0]))
        ucloud_hostids.append(str(req.json()["UHostIds"][0]))

    machine_fin = [False for ii in ucloud_machines]
    total_machine_num = len(ucloud_machines)
    fin_machine_num = 0
    while not all(machine_fin):
        for idx,mac in enumerate(ucloud_machines):
            if not machine_fin[idx]:
                ucloud_check_param = {}
                ucloud_check_param['Action'] = "GetUHostInstanceVncInfo"
                ucloud_check_param['Region'] = machine['ucloud_param']['Region']
                ucloud_check_param['UHostId'] = ucloud_hostids[idx]
                ucloud_check_param['PublicKey'] = machine['ucloud_param']['PublicKey']
                ucloud_check_param['Signature'] = _verfy_ac(machine['Private'], ucloud_check_param)
                req = requests.get(ucloud_url, ucloud_check_param)
                print("the UHostId is", ucloud_hostids[idx])
                print(json.dumps(req.json(),indent=2, sort_keys=True))
                if req.json()['RetCode'] == 0 :
                    machine_fin[idx] = True
                    fin_machine_num = fin_machine_num + 1
        print("Current finish",fin_machine_num,"/", total_machine_num)
        time.sleep(10)
    
    ssh_sess = []
    ssh_param = {}
    ssh_param['port'] = 22
    ssh_param['username'] = 'root'
    ssh_param['work_path'] = machine['work_path']
    for ii in ucloud_machines :
        ssh_param['hostname'] = ii
        ssh_sess.append(SSHSession(ssh_param))

    job_list = []
    for ii in range(njob) :
        chunk = task_chunks[ii]
        rjob = CloudMachineJob(ssh_sess[ii], work_path)
        rjob.upload('.',  forward_common_files)
        rjob.upload(chunk, forward_task_files)
        rjob.submit(chunk, command)
        job_list.append(rjob)
    
    job_fin = [False for ii in job_list]
    while not all(job_fin) :
        for idx,rjob in enumerate(job_list) :
            if not job_fin[idx] :
                status = rjob.check_status()
                if status == JobStatus.terminated :
                    raise RuntimeError("find unsuccessfully terminated job on machine" % ucloud_machines[idx])
                elif status == JobStatus.finished :
                    rjob.download(task_chunks[idx], backward_task_files)
                    rjob.clean()
                    _ucloud_remove_machine(machine, ucloud_hostids[idx])
                    job_fin[idx] = True
        time.sleep(10)


def _group_submit_jobs(ssh_sess,
                       resources,
                       command,
                       work_path,
                       tasks,
                       group_size,
                       forward_common_files,
                       forward_task_files,
                       backward_task_files) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    job_list = []
    for chunk in task_chunks :
        rjob = SlurmJob(ssh_sess, work_path)
        rjob.upload('.',  forward_common_files)
        rjob.upload(chunk, forward_task_files)
        rjob.submit(resources, chunk, command)
        job_list.append(rjob)

    job_fin = [False for ii in job_list]
    while not all(job_fin) :
        for idx,rjob in enumerate(job_list) :
            if not job_fin[idx] :
                status = rjob.check_status()
                if status == JobStatus.terminated :
                    raise RuntimeError("find unsuccessfully terminated job in %s" % rjob.get_job_root())
                elif status == JobStatus.finished :
                    rjob.download(task_chunks[idx], backward_task_files)
                    rjob.clean()
                    job_fin[idx] = True
        time.sleep(10)

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

def make_fp_task_name(sys_idx, counter) : 
    return 'task.' + fp_task_fmt % (sys_idx, counter)

def get_sys_index(task) :
    task.sort()
    system_index = []
    for ii in task :
        system_index.append(os.path.basename(ii).split('.')[1])
    set_tmp = set(system_index)
    system_index = list(set_tmp)
    system_index.sort()
    return system_index
    
def check_empty_iter(iter_index, max_v = 0) :
    fp_path = os.path.join(make_iter_name(iter_index), fp_name)
    fp_tasks = glob.glob(os.path.join(fp_path, "task.*"))
    sys_index = get_sys_index(fp_tasks)
    empty_sys = []
    for ii in sys_index:
        sys_tasks = glob.glob(os.path.join(fp_path, "task." + ii + ".*"))
        empty_sys.append(len(sys_tasks) < max_v)
    return all(empty_sys)

def copy_model(numb_model, prv_iter_index, cur_iter_index) :
    cwd=os.getcwd()
    prv_train_path = os.path.join(make_iter_name(prv_iter_index), train_name)
    cur_train_path = os.path.join(make_iter_name(cur_iter_index), train_name)
    prv_train_path = os.path.abspath(prv_train_path)
    cur_train_path = os.path.abspath(cur_train_path)
    create_path(cur_train_path)
    for ii in range(numb_model) :
        prv_train_task = os.path.join(prv_train_path, train_task_fmt%ii)
        os.chdir(cur_train_path)
        os.symlink(os.path.relpath(prv_train_task), train_task_fmt%ii)
        os.symlink(os.path.join(train_task_fmt%ii, 'frozen_model.pb'), 'graph.%03d.pb' % ii)
        os.chdir(cwd)
    with open(os.path.join(cur_train_path, "copied"), 'w') as fp:
        None 

def poscar_natoms(lines) :
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    return numb_atoms

def poscar_shuffle(poscar_in, poscar_out) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    numb_atoms = poscar_natoms(lines)
    idx = np.arange(8, 8+numb_atoms)
    np.random.shuffle(idx)
    out_lines = lines[0:8]
    for ii in range(numb_atoms) :
        out_lines.append(lines[idx[ii]])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(out_lines))

def expand_idx (in_list) :
    ret = []
    for ii in in_list :
        if type(ii) == int :
            ret.append(ii)
        elif type(ii) == str:
            step_str = ii.split(':')
            if len(step_str) > 1 :
                step = int(step_str[1])
            else :
                step = 1
            range_str = step_str[0].split('-')
            assert(len(range_str)) == 2
            ret += range(int(range_str[0]), int(range_str[1]), step)
    return ret
        
def make_train (iter_index, 
               jdata) :    
    # load json param
    train_param = jdata['train_param']
    if iter_index > 0 :
        stop_batch = jdata['res_stop_batch']
        start_lr = jdata['res_start_lr']
        decay_steps = jdata['res_decay_steps']
        decay_rate = jdata['res_decay_rate']
    numb_models = jdata['numb_models']
    init_data_sys_ = jdata['init_data_sys']    
    fp_task_min = jdata['fp_task_min']    
    
    if iter_index > 0 and check_empty_iter(iter_index-1, fp_task_min) :
        log_task('prev data is empty, copy prev model')
        copy_model(numb_models, iter_index-1, iter_index)
        return
    else :
        iter_name = make_iter_name(iter_index)
        work_path = os.path.join(iter_name, train_name)
        copy_flag = os.path.join(work_path, 'copied')
        if os.path.isfile(copy_flag) :
            os.remove(copy_flag)

    init_data_sys = []
    init_batch_size = list(jdata['init_batch_size'])
    sys_batch_size = jdata['sys_batch_size']
    for ii in init_data_sys_ :
        init_data_sys.append(os.path.abspath(ii))
    if iter_index > 0 :
        for ii in range(iter_index) :
            fp_path = os.path.join(make_iter_name(ii), fp_name)
            fp_data_sys = glob.glob(os.path.join(fp_path, "data.*"))            
            for jj in fp_data_sys :
                tmp_box = np.loadtxt(os.path.join(jj, 'box.raw'))
                tmp_box = np.reshape(tmp_box, [-1,9])
                nframes = tmp_box.shape[0]
                if nframes < fp_task_min :
                    log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                    continue
                init_data_sys.append(os.path.abspath(jj))
                sys_idx = int(jj.split('.')[-1])
                init_batch_size.append(sys_batch_size[sys_idx])                
    for ii in init_data_sys :
        if not os.path.isdir(ii) :
            raise RuntimeError ("data sys %s does not exists, cwd is %s" % (ii, os.getcwd()))
    # establish work path
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    create_path(work_path)
    # establish tasks
    jinput = jdata['default_training_param']
    jinput['systems'] = init_data_sys    
    jinput['batch_size'] = init_batch_size
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(task_path)
        jinput['seed'] = random.randrange(sys.maxsize)
        with open(os.path.join(task_path, train_param), 'w') as outfile:
            json.dump(jinput, outfile, indent = 4)
    # link old models
    if iter_index > 0 :
        prev_iter_name = make_iter_name(iter_index-1)
        prev_work_path = os.path.join(prev_iter_name, train_name)
        for ii in range(numb_models) :
            prev_task_path =  os.path.join(prev_work_path, train_task_fmt%ii)
            old_model_files = glob.glob(
                os.path.join(prev_task_path, "model.ckpt*"))
            task_path = os.path.join(work_path, train_task_fmt % ii)
            task_old_path = os.path.join(task_path, 'old')
            create_path(task_old_path)
            cwd = os.getcwd()
            for jj in old_model_files:
                absjj = os.path.abspath(jj)
                basejj = os.path.basename(jj)
                os.chdir(task_old_path)
                os.symlink(os.path.relpath(absjj), basejj)
                os.chdir(cwd)            

def run_train (iter_index,
               jdata, 
               absmachine) :    
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = jdata['deepmd_path']
    train_param = jdata['train_param']
    train_resources = jdata['train_resources']

    # paths
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    # check if is copied
    copy_flag = os.path.join(work_path, 'copied')
    if os.path.isfile(copy_flag) :
        log_task('copied model, do not train')
        return
    # make tasks
    all_task = []
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        all_task.append(task_path)
    command =  os.path.join(deepmd_path, 'bin/dp_train')
    command += ' %s && ' % train_param
    command += os.path.join(deepmd_path, 'bin/dp_frz')

    run_tasks = [os.path.basename(ii) for ii in all_task]
    forward_files = [train_param]
    backward_files = ['frozen_model.pb', 'lcurve.out']

    if (type(absmachine) == dict) and \
       ('machine_type' in absmachine) and  \
       (absmachine['machine_type'] == 'ucloud') :
        _ucloud_submit_jobs(absmachine,
                            command, 
                            work_path,
                            run_tasks,
                            1,
                            [],
                            forward_files,
                            backward_files)
    else :
        _group_submit_jobs(absmachine,
                           train_resources,
                           command,
                           work_path,
                           run_tasks,
                           1,
                           [],
                           forward_files,
                           backward_files)

    # exec_batch(command,
    #            1,
    #            train_nppn,
    #            train_ngpu,
    #            all_task,
    #            time_limit = train_tlimit,
    #            modules = train_modules,
    #            sources = train_sources)

def post_train (iter_index,
                jdata) :
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = jdata['deepmd_path']
    # paths
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    # check if is copied
    copy_flag = os.path.join(work_path, 'copied')
    if os.path.isfile(copy_flag) :
        log_task('copied model, do not post train')
        return
    all_task = []
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        all_task.append(task_path)
    command = os.path.join(deepmd_path, 'bin/dp_frz')
    command = cmd_append_log(command, 'freeze.log')
    command = 'CUDA_VISIBLE_DEVICES="" ' + command
    # frz models
    # exec_hosts(MachineLocal, command, 1, all_task, None)
    # symlink models
    for ii in range(numb_models) :
        task_file = os.path.join(train_task_fmt % ii, 'frozen_model.pb')
        ofile = os.path.join(work_path, 'graph.%03d.pb' % ii)
        if os.path.isfile(ofile) :
            os.remove(ofile)
        os.symlink(task_file, ofile)    

def make_model_devi (iter_index, 
                     jdata) :
    model_devi_dt = jdata['model_devi_dt']
    model_devi_jobs = jdata['model_devi_jobs']
    assert (iter_index < len(model_devi_jobs))
    cur_job = model_devi_jobs[iter_index]
    # ensemble = model_devi_jobs['ensemble']
    # nsteps = model_devi_jobs['nsteps']
    # trj_freq = model_devi_jobs['trj_freq']
    # job_names = get_job_names (model_devi_jobs)
    # assert (iter_index < len(job_names)) 
    # cur_job_name = job_names[iter_index]    
    # cur_job = model_devi_jobs[cur_job_name]
    ensemble = cur_job['ensemble']
    nsteps = cur_job['nsteps']
    trj_freq = cur_job['trj_freq']
    sys_configs = jdata['sys_configs']
    shuffle_poscar = jdata['shuffle_poscar']

    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")
    conf_systems = []
    for idx in sys_idx :
        cur_systems = []
        ss = sys_configs[idx]
        for ii in ss :
            cur_systems += glob.glob(ii)
        cur_systems.sort()
        conf_systems.append (cur_systems)
    mass_map = jdata['mass_map']
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
    with open(os.path.join(work_path, 'cur_job.json'), 'w') as outfile:
        json.dump(cur_job, outfile, indent = 4)

    all_task = []
    task_param = []
    conf_path = os.path.join(work_path, 'confs')
    create_path(conf_path)
    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        for cc in ss :            
            conf_name = make_model_devi_conf_name(sys_idx[sys_counter], conf_counter)
            orig_poscar_name = conf_name + '.orig.poscar'
            poscar_name = conf_name + '.poscar'
            lmp_name = conf_name + '.lmp'
            if shuffle_poscar :
                os.symlink(cc, os.path.join(conf_path, orig_poscar_name))
                poscar_shuffle(os.path.join(conf_path, orig_poscar_name),
                               os.path.join(conf_path, poscar_name))
            else :
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
                    task_name = make_model_devi_task_name(sys_idx[sys_counter], task_counter)
                    conf_name = make_model_devi_conf_name(sys_idx[sys_counter], conf_counter) + '.lmp'
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
                                               model_devi_dt,
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
                    absmachine) :
    lmp_exec = jdata['lmp_command']
    model_devi_group_size = jdata['model_devi_group_size']
    model_devi_resources = jdata['model_devi_resources']

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = lmp_exec + " -i input.lammps"
    command = cmd_append_log(command, "model_devi.log")

    fp = open (os.path.join(work_path, 'cur_job.json'), 'r')
    cur_job = json.load (fp)
    traj_freq = cur_job['trj_freq']
    nsteps = cur_job['nsteps']
    nframes = nsteps // traj_freq + 1
    
    run_tasks_ = []
    for ii in all_task:
        fres = os.path.join(ii, 'model_devi.out')
        if os.path.isfile(fres) :
            nlines = np.loadtxt(fres).shape[0]
            if nframes != nlines :
                run_tasks_.append(ii)
        else :
            run_tasks_.append(ii)

    run_tasks = [os.path.basename(ii) for ii in run_tasks_]
    all_models = glob.glob(os.path.join(work_path, 'graph*pb'))
    model_names = [os.path.basename(ii) for ii in all_models]
    forward_files = ['conf.lmp', 'input.lammps', 'traj']
    backward_files = ['model_devi.out', 'model_devi.log', 'traj']

    if (type(absmachine) == dict) and \
       ('machine_type' in absmachine) and  \
       (absmachine['machine_type'] == 'ucloud') :
        _ucloud_submit_jobs(absmachine,
                            command,
                            work_path,
                            run_tasks,
                            model_devi_group_size,
                            model_names,
                            forward_files,
                            backward_files)
    else :
        _group_submit_jobs(absmachine,
                           model_devi_resources,
                           command,
                           work_path,
                           run_tasks,
                           model_devi_group_size,
                           model_names,
                           forward_files,
                           backward_files)

    # exec_hosts_batch(exec_machine, command, model_devi_np, run_tasks, None, verbose = True, gpu = True)
    # exec_batch_group(command,
    #                  model_devi_nn, model_devi_nppn, model_devi_ngpu,
    #                  run_tasks,
    #                  group_size = model_devi_group_size,
    #                  time_limit = model_devi_tlimit,
    #                  modules = model_devi_modules,
    #                  sources = model_devi_sources)

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

def _make_fp_vasp_inner (modd_path,
                         work_path,
                         model_devi_skip,
                         e_trust_lo,
                         e_trust_hi,
                         f_trust_lo,
                         f_trust_hi,
                         fp_task_min,
                         fp_task_max,
                         fp_link_files,
                         fp_params):
    """
    modd_path           string          path of model devi
    work_path           string          path of fp
    fp_task_max         int             max number of tasks
    fp_link_files       [string]        linked files for fp, POTCAR for example
    fp_params           map             parameters for fp
    """
    modd_task = glob.glob(os.path.join(modd_path, "task.*"))
    modd_task.sort()
    system_index = []
    for ii in modd_task :
        system_index.append(os.path.basename(ii).split('.')[1])
    set_tmp = set(system_index)
    system_index = list(set_tmp)
    system_index.sort()

    fp_tasks = []
    for ss in system_index :
        fp_candidate = []
        fp_rest = []
        modd_system_glob = os.path.join(modd_path, 'task.' + ss + '.*')
        modd_system_task = glob.glob(modd_system_glob)
        modd_system_task.sort()
        cc = 0
        for tt in modd_system_task :
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                all_conf = np.loadtxt(os.path.join(tt, 'model_devi.out'))
                sel_conf = []
                res_conf = []
                for ii in range(all_conf.shape[0]) :
                    if (all_conf[ii][1] < e_trust_hi and all_conf[ii][1] > e_trust_lo) or \
                       (all_conf[ii][4] < f_trust_hi and all_conf[ii][4] > f_trust_lo) and \
                       ii >= model_devi_skip :
                        sel_conf.append(int(all_conf[ii][0]))
                    else :
                        res_conf.append(int(all_conf[ii][0]))
                for ii in sel_conf:
                    fp_candidate.append([tt, ii])
                for ii in res_conf:
                    fp_rest.append([tt, ii])
        random.shuffle(fp_candidate)
        with open(os.path.join(work_path,'candidate.shuffled.%s.out'%ss), 'w') as fp:
            for ii in fp_candidate:
                fp.write(str(ii[0]) + " " + str(ii[1]) + "\n")
        random.shuffle(fp_rest)
        with open(os.path.join(work_path,'rest.shuffled.%s.out'%ss), 'w') as fp:
            for ii in fp_rest:
                fp.write(str(ii[0]) + " " + str(ii[1]) + "\n")
        numb_task = min(fp_task_max, len(fp_candidate))
        for cc in range(numb_task) :
            tt = fp_candidate[cc][0]
            ii = fp_candidate[cc][1]
            ss = os.path.basename(tt).split('.')[1]
            conf_name = os.path.join(tt, "traj")
            conf_name = os.path.join(conf_name, str(ii) + '.lammpstrj')
            conf_name = os.path.abspath(conf_name)
            fp_task_name = make_fp_task_name(int(ss), cc)
            fp_task_path = os.path.join(work_path, fp_task_name)
            create_path(fp_task_path)
            fp_tasks.append(fp_task_path)
            cwd = os.getcwd()
            os.chdir(fp_task_path)
            os.symlink(os.path.relpath(conf_name), 'conf.lmp')
            for pair in fp_link_files :
                os.symlink(pair[0], pair[1])
            os.chdir(cwd)
        if numb_task < fp_task_min:
            for cc in range(fp_task_min - numb_task) :
                tt = fp_rest[cc][0]
                ii = fp_rest[cc][1]
                ss = os.path.basename(tt).split('.')[1]
                conf_name = os.path.join(tt, "traj")
                conf_name = os.path.join(conf_name, str(ii) + '.lammpstrj')
                conf_name = os.path.abspath(conf_name)
                fp_task_name = make_fp_task_name(int(ss), cc + numb_task)
                fp_task_path = os.path.join(work_path, fp_task_name)
                create_path(fp_task_path)
                fp_tasks.append(fp_task_path)
                cwd = os.getcwd()
                os.chdir(fp_task_path)
                os.symlink(os.path.relpath(conf_name), 'conf.lmp')
                shutil.copyfile(os.path.relpath(conf_name), 'conf.lmp.bk')
                for pair in fp_link_files :
                    os.symlink(pair[0], pair[1])
                os.chdir(cwd)            
    return fp_tasks

def _link_fp_vasp_incar (iter_index,
                         jdata,
                         incar = 'INCAR') :
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)    
    incar_file = os.path.join(work_path, incar)
    incar_file = os.path.abspath(incar_file)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        os.symlink(os.path.relpath(incar_file), incar)
        os.chdir(cwd)

def _link_fp_vasp_pp (iter_index,
                      jdata) :
    fp_pp_path = jdata['fp_pp_path']
    fp_pp_files = jdata['fp_pp_files']
    assert(os.path.exists(fp_pp_path))
    fp_pp_path = os.path.abspath(fp_pp_path)
    
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        for jj in fp_pp_files:
            pp_file = os.path.join(fp_pp_path, jj)
            os.symlink(pp_file, jj)
        os.chdir(cwd)

def _make_fp_vasp_configs(iter_index, 
                          jdata):
    fp_task_max = jdata['fp_task_max']
    fp_params = jdata['fp_params']
    model_devi_skip = jdata['model_devi_skip']
    e_trust_lo = jdata['model_devi_e_trust_lo']
    e_trust_hi = jdata['model_devi_e_trust_hi']
    f_trust_lo = jdata['model_devi_f_trust_lo']
    f_trust_hi = jdata['model_devi_f_trust_hi']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    modd_path = os.path.join(iter_name, model_devi_name)
    cur_job = json.load(open(os.path.join(modd_path, 'cur_job.json'), 'r'))
    task_min = -1
    if 'task_min' in cur_job :
        task_min = cur_job['task_min']
    # make configs
    fp_tasks = _make_fp_vasp_inner(modd_path, work_path,
                                   model_devi_skip,
                                   e_trust_lo, e_trust_hi,
                                   f_trust_lo, f_trust_hi,
                                   task_min, fp_task_max,
                                   [],
                                   fp_params)
    return fp_tasks

def _fix_poscar_type (jdata, task_dirs) :
    type_map = jdata['type_map']
    for ii in task_dirs :
        poscar_file = os.path.join(ii, 'POSCAR')
        for idx,jj in enumerate(type_map):
            old_str = 'Type_%d' % (idx+1)
            new_str = jj
            replace(poscar_file, old_str, new_str)
    
def make_fp_vasp (iter_index, 
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return        
    #convert configs
    command = os.path.join(os.getcwd(), "lib/ovito_file_convert.py")
    command += " conf.lmp POSCAR"
    exec_hosts(MachineLocal, command, 1, fp_tasks, verbose = True)
    _fix_poscar_type(jdata, fp_tasks)
    # create incar
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    incar = make_vasp_incar(jdata['fp_params'])
    incar_file = os.path.join(work_path, 'INCAR')
    incar_file = os.path.abspath(incar_file)
    open(incar_file, 'w').write(incar)
    _link_fp_vasp_incar(iter_index, jdata)
    # create potcar
    _link_fp_vasp_pp(iter_index, jdata)
    # clean traj
    clean_traj = True
    if 'model_devi_clean_traj' in jdata :
        clean_traj = jdata['model_devi_clean_traj']
    if clean_traj:
        md_trajs = glob.glob(os.path.join(modd_path, 'task*/traj'))
        for ii in md_trajs :
            shutil.rmtree(ii)

            
def make_fp_pwscf(iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return        
    #convert configs
    command = os.path.join(os.getcwd(), "lib/ovito_file_convert.py")
    command += " conf.lmp POSCAR"
    exec_hosts(MachineLocal, command, 1, fp_tasks, verbose = True)
    _fix_poscar_type(jdata, fp_tasks)
    # make pwscf input
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_pp_files = jdata['fp_pp_files']
    fp_params = jdata['fp_params']
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    kspacing = fp_params['kspacing']
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        sys_data = system_from_poscar('POSCAR')
        sys_data['atom_masses'] = jdata['mass_map']
        ret = make_pwscf_input(sys_data, ecut, ediff, fp_pp_files, kspacing)
        open('input', 'w').write(ret)
        os.chdir(cwd)
    # link pp files
    _link_fp_vasp_pp(iter_index, jdata)
    # clean traj
    clean_traj = True
    if 'model_devi_clean_traj' in jdata :
        clean_traj = jdata['model_devi_clean_traj']
    if clean_traj:
        modd_path = os.path.join(iter_name, model_devi_name)
        md_trajs = glob.glob(os.path.join(modd_path, 'task*/traj'))
        for ii in md_trajs :
            shutil.rmtree(ii)
            
def make_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        make_fp_vasp(iter_index, jdata) 
    elif fp_style == "pwscf" :
        make_fp_pwscf(iter_index, jdata) 
    else :
        raise RuntimeError ("unsupported fp style")

def _vasp_check_fin (ii) :
    if os.path.isfile(os.path.join(ii, 'OUTCAR')) :
        with open(os.path.join(ii, 'OUTCAR'), 'r') as fp :
            content = fp.read()
            count = content.count('Elapse')
            if count != 1 :
                return False
    else :
        return False
    return True

def _qe_check_fin(ii) :
    if os.path.isfile(os.path.join(ii, 'output')) :
        with open(os.path.join(ii, 'output'), 'r') as fp :
            content = fp.read()
            count = content.count('JOB DONE')
            if count != 1 :
                return False
    else :
        return False
    return True
    
        
def run_fp_vasp (iter_index,
                 jdata,
                 absmachine,
                 check_fin,
                 log_file = "log") :
    fp_command = jdata['fp_command']
    fp_group_size = jdata['fp_group_size']
    fp_resources = jdata['fp_resources']
    # fp_command = ("OMP_NUM_THREADS=1 mpirun -n %d " % fp_np) + fp_command
    # cpu task in parallel
    if ('numb_gpu' not in fp_resources) or (fp_resources['numb_gpu'] == 0):
        fp_command = "srun " + fp_command
    fp_command = cmd_append_log(fp_command, log_file)

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    fp_run_tasks = []
    for ii in fp_tasks :
        if not check_fin(ii) :
            fp_run_tasks.append(ii)

    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]
    forward_files = ['POSCAR', 'INCAR', 'POTCAR']
    backward_files = ['OUTCAR']

    if (type(absmachine) == dict) and \
       ('machine_type' in absmachine) and  \
       (absmachine['machine_type'] == 'ucloud') :
        _ucloud_submit_jobs(absmachine,
                            fp_command,
                            work_path,
                            run_tasks,
                            fp_group_size,
                            [],
                            forward_files,
                            backward_files)
    else :
        _group_submit_jobs(absmachine,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files)

def run_fp_pwscf (iter_index,
                 jdata,
                 absmachine,
                 check_fin,
                 log_file = "log") :
    fp_command = jdata['fp_command']
    fp_group_size = jdata['fp_group_size']
    fp_resources = jdata['fp_resources']
    fp_pp_files = jdata['fp_pp_files']
    # fp_command = ("OMP_NUM_THREADS=1 mpirun -n %d " % fp_np) + fp_command
    # cpu task in parallel
    if ('numb_gpu' not in fp_resources) or (fp_resources['numb_gpu'] == 0):
        fp_command = "srun " + fp_command
    fp_command = cmd_append_log(fp_command, log_file)

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    fp_run_tasks = []
    for ii in fp_tasks :
        if not check_fin(ii) :
            fp_run_tasks.append(ii)

    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]
    forward_files = ['input'] + fp_pp_files
    backward_files = ['output']

    if (type(absmachine) == dict) and \
       ('machine_type' in absmachine) and  \
       (absmachine['machine_type'] == 'ucloud') :
        _ucloud_submit_jobs(absmachine,
                            fp_resources,
                            fp_command,
                            work_path,
                            run_tasks,
                            fp_group_size,
                            [],
                            forward_files,
                            backward_files)
    else :
        _group_submit_jobs(absmachine,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files)

#    exec_hosts_batch(exec_machine, fp_command, fp_np, fp_run_tasks, verbose = True, mpi = False, gpu=True)
    # exec_batch_group(fp_command,
    #                  fp_nn, fp_nppn, fp_ngpu,
    #                  fp_run_tasks,
    #                  group_size = fp_group_size,
    #                  time_limit = fp_tlimit,
    #                  modules = fp_modules,
    #                  sources = fp_sources)
        
def run_fp (iter_index,
            jdata,
            exec_machine) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        run_fp_vasp(iter_index, jdata, exec_machine, _vasp_check_fin) 
    elif fp_style == "pwscf" :
        run_fp_pwscf(iter_index, jdata, exec_machine, _qe_check_fin, log_file = 'output') 
    else :
        raise RuntimeError ("unsupported fp style") 


def post_fp_vasp (iter_index,
                  jdata,
                  to_config = 'template/tools.vasp/cessp2force_lin.py'):
    model_devi_jobs = jdata['model_devi_jobs']
    assert (iter_index < len(model_devi_jobs)) 
    cur_job = model_devi_jobs[iter_index]
    # job_names = get_job_names (model_devi_jobs)
    # assert (iter_index < len(job_names)) 
    # cur_job_name = job_names[iter_index]    
    # cur_job = model_devi_jobs[cur_job_name]

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    system_index = []
    for ii in fp_tasks :        
        system_index.append(os.path.basename(ii).split('.')[1])
    system_index.sort()
    set_tmp = set(system_index)
    system_index = list(set_tmp)
    system_index.sort()
    
    to_config = os.path.abspath(to_config)
    cmd_to_config = to_config + " OUTCAR "
    cmd_to_config = cmd_append_log(cmd_to_config, "to_config.log")
    exec_hosts(MachineLocal, cmd_to_config, 1, fp_tasks)

    convert_to_raw = 'template/tools.vasp/convert2raw.py'
    shuffle_raw = 'template/tools.raw/shuffle_raw.py'
    copy_raw = 'template/tools.raw/copy_raw.py'
    raw_to_set = 'template/tools.raw/raw_to_set.sh'
    convert_to_raw = os.path.abspath(convert_to_raw)
    shuffle_raw = os.path.abspath(shuffle_raw)
    copy_raw = os.path.abspath(copy_raw)
    raw_to_set = os.path.abspath(raw_to_set)
    cwd = os.getcwd()
    for ss in system_index :
        sys_config_data = glob.glob(os.path.join(work_path, "task.%s.*/test.configs"%ss))
        sys_data_path = os.path.join(work_path, 'data.%s/orig'%ss)
        create_path(sys_data_path)
        with open(os.path.join(sys_data_path, 'data.configs'), 'wb') as wfd:
            for f in sys_config_data :
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)
        os.chdir(sys_data_path)
        sp.check_call(convert_to_raw + ' data.configs', shell = True)
        os.chdir('..')
        sp.check_call(shuffle_raw + ' orig/ .', shell = True)
        if os.path.isfile('type.raw') :
            os.remove('type.raw')
        os.symlink('orig/type.raw', 'type.raw')
        sp.check_call(raw_to_set, shell = True)
        os.chdir(cwd)

def post_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        post_fp_vasp(iter_index, jdata) 
    elif fp_style == "pwscf" :
        post_fp_vasp(iter_index, jdata,
                     to_config = 'template/tools.pwscf/pwscf1frame.py') 
    else :
        raise RuntimeError ("unsupported fp style")            
    
def run_iter (json_file, exec_machine) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_iter = jdata["numb_iter"]
    max_tasks = 10000
    numb_task = 9
    record = "record.dpgen"

    train_machine = jdata['train_machine']    
    if ('machine_type' in train_machine) and  \
       (train_machine['machine_type'] == 'ucloud'):
        train_absmachine = train_machine
    else :
        train_absmachine = SSHSession(train_machine)

    model_devi_machine = jdata['model_devi_machine']    
    if ('machine_type' in model_devi_machine) and  \
       (model_devi_machine['machine_type'] == 'ucloud'):
        model_devi_absmachine = model_devi_machine
    else :
        model_devi_absmachine = SSHSession(model_devi_machine)

    fp_machine = jdata['fp_machine']    
    if ('machine_type' in fp_machine) and  \
       (fp_machine['machine_type'] == 'ucloud'):
        fp_absmachine = fp_machine
    else :
        fp_absmachine = SSHSession(fp_machine)

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        logging.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    for ii in range (numb_iter) :
        for jj in range (numb_task) :
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] : 
                continue
            if   jj == 0 :
                log_iter ("make_train", ii, jj)
                make_train (ii, jdata) 
            elif jj == 1 :
                log_iter ("run_train", ii, jj)
                run_train  (ii, jdata, train_absmachine)
            elif jj == 2 :
                log_iter ("post_train", ii, jj)
                post_train  (ii, jdata)
            elif jj == 3 :
                log_iter ("make_model_devi", ii, jj)
                make_model_devi  (ii, jdata)
            elif jj == 4 :
                log_iter ("run_model_devi", ii, jj)
                run_model_devi  (ii, jdata, model_devi_absmachine)
            elif jj == 5 :
                log_iter ("post_model_devi", ii, jj)
                post_model_devi  (ii, jdata)
            elif jj == 6 :
                log_iter ("make_fp", ii, jj)
                make_fp (ii, jdata)
            elif jj == 7 :
                log_iter ("run_fp", ii, jj)
                run_fp (ii, jdata, fp_absmachine)
            elif jj == 8 :
                log_iter ("post_fp", ii, jj)
                post_fp (ii, jdata)
            else :
                raise RuntimeError ("unknow task %d, something wrong" % jj)
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
    logging.getLogger("paramiko").setLevel(logging.WARNING)
    paramiko.util.log_to_file("<log_file_path>", level = "WARN")

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
    run_iter('param.json', MachineLocalGPU)
    # run_iter('param.json', MachineLocal)

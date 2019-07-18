#!/usr/bin/env python3

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
import dpdata
import numpy as np
import subprocess as sp
from dpgen import dlog
from dpgen.generator.lib.utils import make_iter_name
from dpgen.generator.lib.utils import create_path
from dpgen.generator.lib.utils import copy_file_list
from dpgen.generator.lib.utils import replace
from dpgen.generator.lib.utils import cmd_append_log
from dpgen.generator.lib.utils import log_iter
from dpgen.generator.lib.utils import record_iter
from dpgen.generator.lib.utils import log_task
from dpgen.generator.lib.lammps import make_lammps_input
from dpgen.generator.lib.vasp import write_incar_dict
from dpgen.generator.lib.vasp import make_vasp_incar_user_dict
from dpgen.generator.lib.pwscf import make_pwscf_input
from dpgen.generator.lib.pwscf import cvt_1frame
from dpgen.generator.lib.gaussian import make_gaussian_input, take_cluster
from dpgen.remote.RemoteJob import SSHSession, JobStatus, SlurmJob, PBSJob, LSFJob, CloudMachineJob
from dpgen.remote.group_jobs import ucloud_submit_jobs
from dpgen.remote.group_jobs import group_slurm_jobs
from dpgen.remote.group_jobs import group_local_jobs
from dpgen.util import sepline
from dpgen import ROOT_PATH
from pymatgen.io.vasp import Incar,Kpoints,Potcar
from dpgen.auto_test.lib.vasp import make_kspacing_kpoints

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
cvasp_file=os.path.join(ROOT_PATH,'generator/lib/cvasp.py')

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
    
def _check_empty_iter(iter_index, max_v = 0) :
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

def _check_skip_train(job) :
    try :
        skip = _get_param_alias(job, ['s_t', 'sk_tr', 'skip_train', 'skip_training']) 
    except ValueError :
        skip = False
    return skip


def poscar_to_conf(poscar, conf):
    sys = dpdata.System(poscar, fmt = 'vasp/poscar')
    sys.to_lammps_lmp(conf)


def dump_to_poscar(dump, poscar, type_map) :
    sys = dpdata.System(dump, fmt = 'lammps/dump', type_map = type_map)
    sys.to_vasp_poscar(poscar)

        
def make_train (iter_index, 
                jdata, 
                mdata) :    
    # load json param
    train_param = jdata['train_param']
    numb_models = jdata['numb_models']
    init_data_prefix = jdata['init_data_prefix']    
    init_data_prefix = os.path.abspath(init_data_prefix)
    init_data_sys_ = jdata['init_data_sys']    
    fp_task_min = jdata['fp_task_min']    
    model_devi_jobs = jdata['model_devi_jobs']
    
    if iter_index > 0 and _check_empty_iter(iter_index-1, fp_task_min) :
        log_task('prev data is empty, copy prev model')
        copy_model(numb_models, iter_index-1, iter_index)
        return
    elif iter_index > 0 and _check_skip_train(model_devi_jobs[iter_index-1]):
        log_task('skip training at step %d ' % (iter_index-1))
        copy_model(numb_models, iter_index-1, iter_index)
        return
    else :
        iter_name = make_iter_name(iter_index)
        work_path = os.path.join(iter_name, train_name)
        copy_flag = os.path.join(work_path, 'copied')
        if os.path.isfile(copy_flag) :
            os.remove(copy_flag)

    # establish work path
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    create_path(work_path)
    # link init data
    cwd = os.getcwd()
    os.chdir(work_path)
    os.symlink(os.path.abspath(init_data_prefix), 'data.init')
    # link iter data
    os.mkdir('data.iters')
    os.chdir('data.iters')
    for ii in range(iter_index) :
        os.symlink(os.path.relpath(os.path.join(cwd, make_iter_name(ii))), make_iter_name(ii))
    os.chdir(cwd)

    init_data_sys = []
    init_batch_size = []
    init_batch_size_ = list(jdata['init_batch_size'])
    sys_batch_size = jdata['sys_batch_size']
    for ii, ss in zip(init_data_sys_, init_batch_size_) :
        if 'init_multi_systems' in jdata and jdata['init_multi_systems']:
            for single_sys in os.listdir(os.path.join(work_path, 'data.init', ii)):
                init_data_sys.append(os.path.join('..', 'data.init', ii, single_sys))
                init_batch_size.append(ss)
        else:
            init_data_sys.append(os.path.join('..', 'data.init', ii))
            init_batch_size.append(ss)
    if iter_index > 0 :
        for ii in range(iter_index) :
            fp_path = os.path.join(make_iter_name(ii), fp_name)
            fp_data_sys = glob.glob(os.path.join(fp_path, "data.*"))            
            for jj in fp_data_sys :
                sys_idx = int(jj.split('.')[-1])
                if 'use_clusters' in jdata and jdata['use_clusters']:
                    nframes = 0
                    for sys_single in os.listdir(jj):
                        tmp_box = np.loadtxt(os.path.join(jj, sys_single, 'box.raw'))
                        tmp_box = np.reshape(tmp_box, [-1,9])
                        nframes += tmp_box.shape[0]
                    if nframes < fp_task_min :
                        log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                        continue
                    for sys_single in os.listdir(os.path.join(jj)):
                        init_data_sys.append(os.path.join('..', 'data.iters', jj, sys_single))
                        init_batch_size.append(sys_batch_size[sys_idx])
                else:
                    tmp_box = np.loadtxt(os.path.join(jj, 'box.raw'))
                    tmp_box = np.reshape(tmp_box, [-1,9])
                    nframes = tmp_box.shape[0]
                    if nframes < fp_task_min :
                        log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                        continue
                    init_data_sys.append(os.path.join('..', 'data.iters', jj))
                    init_batch_size.append(sys_batch_size[sys_idx])
    # establish tasks
    jinput = jdata['default_training_param']
    jinput['systems'] = init_data_sys    
    jinput['batch_size'] = init_batch_size
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(task_path)
        os.chdir(task_path)
        for ii in init_data_sys :
            if not os.path.isdir(ii) :
                raise RuntimeError ("data sys %s does not exists, cwd is %s" % (ii, os.getcwd()))
        os.chdir(cwd)
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
               mdata, 
               ssh_sess) :
    # load json param
    numb_models = jdata['numb_models']
    train_param = jdata['train_param']
    deepmd_path = mdata['deepmd_path']
    train_resources = mdata['train_resources']
    machine_type = mdata['train_machine']['machine_type']

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
    init_data_sys_ = jdata['init_data_sys']
    init_data_sys = []
    for ii in init_data_sys_ :
        init_data_sys.append(os.path.join('data.init', ii))
    fp_data_ = glob.glob(os.path.join('iter.*', '02.fp', 'data.*'))
    fp_data = []
    for ii in fp_data_:
        fp_data.append(os.path.join('data.iters', ii))
    trans_comm_data = []
    cwd = os.getcwd()
    os.chdir(work_path)
    for ii in init_data_sys :
        if 'init_multi_systems' in jdata and jdata['init_multi_systems']:
            for single_sys in os.listdir(os.path.join(ii)):
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'type.raw'))
        else:
            trans_comm_data += glob.glob(os.path.join(ii, 'set.*'))
            trans_comm_data += glob.glob(os.path.join(ii, 'type.raw'))
    for ii in fp_data :
        if 'use_clusters' in jdata and jdata['use_clusters']:
            for single_sys in os.listdir(os.path.join(ii)):
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'type.raw'))
        else:
            trans_comm_data += glob.glob(os.path.join(ii, 'set.*'))
            trans_comm_data += glob.glob(os.path.join(ii, 'type.raw'))
    os.chdir(cwd)

    if ssh_sess == None and machine_type == 'ucloud':
        ucloud_submit_jobs(mdata['train_machine'],
                           mdata['train_resources'],
                           command, 
                           work_path,
                           run_tasks,
                           1,
                           trans_comm_data,
                           forward_files,
                           backward_files)
    elif machine_type == 'slurm' :        
        group_slurm_jobs(ssh_sess,
                         train_resources,
                         command,
                         work_path,
                         run_tasks,
                         1,
                         trans_comm_data,
                         forward_files,
                         backward_files)
    elif machine_type == 'pbs' :        
        group_slurm_jobs(ssh_sess,
                         train_resources,
                         command,
                         work_path,
                         run_tasks,
                         1,
                         trans_comm_data,
                         forward_files,
                         backward_files,
                         remote_job = PBSJob)
    elif machine_type == 'lsf':        
        group_slurm_jobs(ssh_sess,
                           train_resources,
                           command,
                           work_path,
                           run_tasks,
                           1,
                           trans_comm_data,
                           forward_files,
                           backward_files,
                          remote_job = LSFJob)
    elif machine_type == 'local' :
        group_local_jobs(ssh_sess,
                         train_resources,
                         command,
                         work_path,
                         run_tasks,
                         1,
                         trans_comm_data,
                         forward_files,
                         backward_files)
    else :
        raise RuntimeError("unknow machine type")

    # exec_batch(command,
    #            1,
    #            train_nppn,
    #            train_ngpu,
    #            all_task,
    #            time_limit = train_tlimit,
    #            modules = train_modules,
    #            sources = train_sources)

def post_train (iter_index,
                jdata,
                mdata) :
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = mdata['deepmd_path']
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

def _get_param_alias(jdata, 
                     names) :
    for ii in names :
        if ii in jdata :
            return jdata[ii]
    raise ValueError("one of the keys %s should be in jdata %s" % (str(names), (json.dumps(jdata, indent=4))))

def parse_cur_job(cur_job) :
    ensemble = _get_param_alias(cur_job, ['ens', 'ensemble'])
    temps = [-1]
    press = [-1]
    if 'npt' in ensemble :
        temps = _get_param_alias(cur_job, ['Ts','temps'])
        press = _get_param_alias(cur_job, ['Ps','press'])
    elif 'nvt' == ensemble :
        temps = _get_param_alias(cur_job, ['Ts','temps'])
    nsteps = _get_param_alias(cur_job, ['nsteps'])
    trj_freq = _get_param_alias(cur_job, ['t_freq', 'trj_freq','traj_freq'])
    if 'pka_e' in cur_job :
        pka_e = _get_param_alias(cur_job, ['pka_e'])
    else :
        pka_e = None
    if 'dt' in cur_job :
        dt = _get_param_alias(cur_job, ['dt'])
    else :
        dt = None
    return ensemble, nsteps, trj_freq, temps, press, pka_e, dt

def make_model_devi (iter_index, 
                     jdata, 
                     mdata) :
    model_devi_dt = jdata['model_devi_dt']
    model_devi_neidelay = None
    if 'model_devi_neidelay' in jdata :
        model_devi_neidelay = jdata['model_devi_neidelay']
    model_devi_taut = 0.1
    if 'model_devi_taut' in jdata :
        model_devi_taut = jdata['model_devi_taut']
    model_devi_taup = 0.5
    if 'model_devi_taup' in jdata :
        model_devi_taup = jdata['model_devi_taup']
    model_devi_jobs = jdata['model_devi_jobs']
    if (iter_index >= len(model_devi_jobs)) :
        return False
    cur_job = model_devi_jobs[iter_index]
    # ensemble = model_devi_jobs['ensemble']
    # nsteps = model_devi_jobs['nsteps']
    # trj_freq = model_devi_jobs['trj_freq']
    # job_names = get_job_names (model_devi_jobs)
    # assert (iter_index < len(job_names)) 
    # cur_job_name = job_names[iter_index]    
    # cur_job = model_devi_jobs[cur_job_name]
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt = parse_cur_job(cur_job)
    if dt is not None :
        model_devi_dt = dt
    if "sys_configs_prefix" in jdata:
        sys_configs = []
        for sys_list in jdata["sys_configs"]:
            #assert (isinstance(sys_list, list) ), "Currently only support type list for sys in 'sys_conifgs' "
            temp_sys_list = [os.path.join(jdata["sys_configs_prefix"], sys) for sys in sys_list]
            sys_configs.append(temp_sys_list)
    else:
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
        cur_systems = [os.path.abspath(ii) for ii in cur_systems]
        conf_systems.append (cur_systems)
    mass_map = jdata['mass_map']

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
            if 'sys_format' in jdata:
                fmt = jdata['sys_format']
            else:
                fmt = 'vasp/poscar'
            system = dpdata.System(os.path.join(conf_path, poscar_name), fmt = fmt)
            system.to_lammps_lmp(os.path.join(conf_path, lmp_name))
            conf_counter += 1
        sys_counter += 1

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
                    cwd_ = os.getcwd()
                    os.chdir(task_path)
                    file_c = make_lammps_input(ensemble,
                                               loc_conf_name,
                                               task_model_list,
                                               nsteps,
                                               model_devi_dt,
                                               model_devi_neidelay,                                               
                                               trj_freq,
                                               mass_map,
                                               tt,
                                               tau_t = model_devi_taut,
                                               pres = pp, 
                                               tau_p = model_devi_taup, 
                                               pka_e = pka_e)
                    job = {}
                    job["ensemble"] = ensemble
                    job["press"] = pp
                    job["temps"] = tt
                    job["model_devi_dt"] =  model_devi_dt
                    with open('job.json', 'w') as _outfile:
                        json.dump(job, _outfile, indent = 4)

                    os.chdir(cwd_)
                    with open(os.path.join(task_path, 'input.lammps'), 'w') as fp :
                        fp.write(file_c)
                    task_counter += 1
            conf_counter += 1
        sys_counter += 1

    return True

def run_model_devi (iter_index, 
                    jdata,
                    mdata,
                    ssh_sess) :
    #rmprint("This module has been run !")
    lmp_exec = mdata['lmp_command']
    model_devi_group_size = mdata['model_devi_group_size']
    model_devi_resources = mdata['model_devi_resources']
    machine_type = mdata['model_devi_machine']['machine_type']

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))

    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = lmp_exec + " -i input.lammps"
    command = cmd_append_log(command, "model_devi.log")

    fp = open (os.path.join(work_path, 'cur_job.json'), 'r')
    cur_job = json.load (fp)
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt = parse_cur_job(cur_job)
    nframes = nsteps // trj_freq + 1
    
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
    #print("all_task is ", all_task)
    #print("run_tasks in run_model_deviation",run_tasks_)
    all_models = glob.glob(os.path.join(work_path, 'graph*pb'))
    model_names = [os.path.basename(ii) for ii in all_models]
    forward_files = ['conf.lmp', 'input.lammps', 'traj']
    backward_files = ['model_devi.out', 'model_devi.log', 'traj']

    dlog.info("group_size %d"%model_devi_group_size)
    if ssh_sess == None and machine_type == 'ucloud':
        dlog.info("The first situation!")
        ucloud_submit_jobs(mdata['model_devi_machine'],
                            mdata['model_devi_resources'],
                            command,
                            work_path,
                            run_tasks,
                            model_devi_group_size,
                            model_names,
                            forward_files,
                            backward_files)
    elif machine_type == 'slurm' :        
        dlog.info("The second situation!")
        group_slurm_jobs(ssh_sess,
                           model_devi_resources,
                           command,
                           work_path,
                           run_tasks,
                           model_devi_group_size,
                           model_names,
                           forward_files,
                           backward_files)
    elif machine_type == 'pbs' :        
        group_slurm_jobs(ssh_sess,
                           model_devi_resources,
                           command,
                           work_path,
                           run_tasks,
                           model_devi_group_size,
                           model_names,
                           forward_files,
                           backward_files,
                          remote_job = PBSJob)
    elif machine_type == 'lsf' :        
        group_slurm_jobs(ssh_sess,
                           model_devi_resources,
                           command,
                           work_path,
                           run_tasks,
                           model_devi_group_size,
                           model_names,
                           forward_files,
                           backward_files,
                          remote_job = LSFJob)
    elif machine_type == 'local' :        
        group_local_jobs(ssh_sess,
                           model_devi_resources,
                           command,
                           work_path,
                           run_tasks,
                           model_devi_group_size,
                           model_names,
                           forward_files,
                           backward_files)
    else :
        raise RuntimeError("unknow machine type")

    # exec_hosts_batch(exec_machine, command, model_devi_np, run_tasks, None, verbose = True, gpu = True)
    # exec_batch_group(command,
    #                  model_devi_nn, model_devi_nppn, model_devi_ngpu,
    #                  run_tasks,
    #                  group_size = model_devi_group_size,
    #                  time_limit = model_devi_tlimit,
    #                  modules = model_devi_modules,
    #                  sources = model_devi_sources)

def post_model_devi (iter_index, 
                     jdata, 
                     mdata) :
    pass

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
                         type_map,
                         cluster_cutoff = None):
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
        fp_rest_accurate = []
        fp_rest_failed = []
        modd_system_glob = os.path.join(modd_path, 'task.' + ss + '.*')
        modd_system_task = glob.glob(modd_system_glob)
        modd_system_task.sort()
        cc = 0
        for tt in modd_system_task :
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                all_conf = np.loadtxt(os.path.join(tt, 'model_devi.out'))
                sel_conf = []
                res_failed_conf = []
                res_accurate_conf = []
                for ii in range(all_conf.shape[0]) :
                    cc = int(all_conf[ii][0])
                    if cluster_cutoff is None:
                        if (all_conf[ii][1] < e_trust_hi and all_conf[ii][1] > e_trust_lo) or \
                        (all_conf[ii][4] < f_trust_hi and all_conf[ii][4] > f_trust_lo) and \
                        ii >= model_devi_skip :
                            fp_candidate.append([tt, cc])
                        elif (all_conf[ii][1] > e_trust_hi ) or (all_conf[ii][4] > f_trust_hi ):
                            fp_rest_failed.append([tt, cc])
                        elif (all_conf[ii][1] < e_trust_lo and all_conf[ii][4] < f_trust_lo ):
                            fp_rest_accurate.append([tt, cc])
                    else:
                        idx_candidate = np.where(np.logical_and(all_conf[ii][7:] < f_trust_hi, all_conf[ii][7:] > f_trust_lo))[0]
                        idx_rest_failed = np.where(all_conf[ii][7:] > f_trust_lo)[0]
                        idx_rest_accurate = np.where(all_conf[ii][7:] < f_trust_hi)[0]
                        for jj in idx_candidate:
                            fp_candidate.append([tt, cc, jj])
                        for jj in idx_rest_accurate:
                            fp_rest_accurate.append([tt, cc, jj])
                        for jj in idx_rest_failed:
                            fp_rest_failed.append([tt, cc, jj])
        random.shuffle(fp_candidate)
        random.shuffle(fp_rest_failed)
        random.shuffle(fp_rest_accurate)
        with open(os.path.join(work_path,'candidate.shuffled.%s.out'%ss), 'w') as fp:
            for ii in fp_candidate:
                fp.write(" ".join([str(nn) for nn in ii]) + "\n")
        with open(os.path.join(work_path,'rest_accurate.shuffled.%s.out'%ss), 'w') as fp:
            for ii in fp_rest_accurate:
                fp.write(" ".join([str(nn) for nn in ii]) + "\n")
        with open(os.path.join(work_path,'rest_failed.shuffled.%s.out'%ss), 'w') as fp:
            for ii in fp_rest_failed:
                fp.write(" ".join([str(nn) for nn in ii]) + "\n")
        numb_task = min(fp_task_max, len(fp_candidate))
        for cc in range(numb_task) :
            tt = fp_candidate[cc][0]
            ii = fp_candidate[cc][1]
            ss = os.path.basename(tt).split('.')[1]
            conf_name = os.path.join(tt, "traj")
            conf_name = os.path.join(conf_name, str(ii) + '.lammpstrj')
            conf_name = os.path.abspath(conf_name)

            # link job.json
            job_name = os.path.join(tt, "job.json")
            job_name = os.path.abspath(job_name)

            if cluster_cutoff is not None:
                # take clusters
                jj = fp_candidate[cc][2]
                poscar_name = '{}.cluster.{}.POSCAR'.format(conf_name, jj)
                new_system = take_cluster(conf_name, type_map, jj, cluster_cutoff)
                new_system.to_vasp_poscar(poscar_name)
            fp_task_name = make_fp_task_name(int(ss), cc)
            fp_task_path = os.path.join(work_path, fp_task_name)
            create_path(fp_task_path)
            fp_tasks.append(fp_task_path)
            cwd = os.getcwd()
            os.chdir(fp_task_path)
            if cluster_cutoff is None:
                os.symlink(os.path.relpath(conf_name), 'conf.dump')
                os.symlink(os.path.relpath(job_name), 'job.json')
            else:
                os.symlink(os.path.relpath(poscar_name), 'POSCAR')
            for pair in fp_link_files :
                os.symlink(pair[0], pair[1])
            os.chdir(cwd)
    if cluster_cutoff is None:
        cwd = os.getcwd()
        for ii in fp_tasks:
            os.chdir(ii)
            dump_to_poscar('conf.dump', 'POSCAR', type_map)
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

def _make_fp_vasp_kp (iter_index,jdata, incar):
    dincar=Incar.from_string(incar)
    standard_incar={}
    for key,val in dincar.items():
        standard_incar[key.upper()]=val
      
    try:
       kspacing = standard_incar['KSPACING'] 
    except:
       raise RuntimeError ("KSPACING must be given in INCAR")
    try:
       gamma = standard_incar['KGAMMA'] 
       if isinstance(gamma,bool):
          pass
       else:
          if gamma[0].upper()=="T":
             gamma=True
          else:
             gamma=False
    except:
       raise RuntimeError ("KGAMMA must be given in INCAR")

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        assert(os.path.exists('POSCAR'))
        ret=make_kspacing_kpoints('POSCAR', kspacing, gamma)
        kp=Kpoints.from_string(ret)
        kp.write_file("KPOINTS")
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
    model_devi_skip = jdata['model_devi_skip']
    e_trust_lo = jdata['model_devi_e_trust_lo']
    e_trust_hi = jdata['model_devi_e_trust_hi']
    f_trust_lo = jdata['model_devi_f_trust_lo']
    f_trust_hi = jdata['model_devi_f_trust_hi']
    type_map = jdata['type_map']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)

    #copy cvasp.py 
    shutil.copyfile(cvasp_file, os.path.join(work_path,'cvasp.py')) 

    modd_path = os.path.join(iter_name, model_devi_name)
    task_min = -1
    if os.path.isfile(os.path.join(modd_path, 'cur_job.json')) :
        cur_job = json.load(open(os.path.join(modd_path, 'cur_job.json'), 'r'))
        if 'task_min' in cur_job :
            task_min = cur_job['task_min']
    # make configs
    cluster_cutoff = jdata['cluster_cutoff'] if 'use_clusters' in jdata and jdata['use_clusters'] else None
    fp_tasks = _make_fp_vasp_inner(modd_path, work_path,
                                   model_devi_skip,
                                   e_trust_lo, e_trust_hi,
                                   f_trust_lo, f_trust_hi,
                                   task_min, fp_task_max,
                                   [],
                                   type_map,
                                   cluster_cutoff)
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
    # create incar
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    
   #fp_params=jdata["fp_params"]

    if 'fp_incar' in jdata.keys() :
        fp_incar_path = jdata['fp_incar']
        assert(os.path.exists(fp_incar_path))
        fp_incar_path = os.path.abspath(fp_incar_path)
        fr = open(fp_incar_path)
        incar = fr.read()
        fr.close()
        #incar= open(fp_incar_path).read()
    elif 'user_fp_params' in jdata.keys() :
        incar = write_incar_dict(jdata['user_fp_params'])
    else:
        incar = make_vasp_incar_user_dict(jdata['fp_params'])
    incar_file = os.path.join(work_path, 'INCAR')
    incar_file = os.path.abspath(incar_file)

    with open(incar_file, 'w') as fp:
        fp.write(incar)
    fp.close()
    _link_fp_vasp_incar(iter_index, jdata)
    # create potcar
    _link_fp_vasp_pp(iter_index, jdata)
    # create kpoints
    _make_fp_vasp_kp(iter_index, jdata, incar)
    # clean traj
    clean_traj = True
    if 'model_devi_clean_traj' in jdata :
        clean_traj = jdata['model_devi_clean_traj']
    if clean_traj:
        modd_path = os.path.join(iter_name, model_devi_name)
        md_trajs = glob.glob(os.path.join(modd_path, 'task*/traj'))
        for ii in md_trajs :
            shutil.rmtree(ii)

            
def make_fp_pwscf(iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return        
    # make pwscf input
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_pp_files = jdata['fp_pp_files']
    if 'user_fp_params' in jdata.keys() :        
        fp_params = jdata['user_fp_params']
        user_input = True
    else:
        fp_params = jdata['fp_params']
        user_input = False        
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        sys_data = dpdata.System('POSCAR').data
        sys_data['atom_masses'] = jdata['mass_map']
        ret = make_pwscf_input(sys_data, fp_pp_files, fp_params, user_input = user_input)
        with open('input', 'w') as fp:
            fp.write(ret)
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


def make_fp_gaussian(iter_index,
                     jdata):
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return        
    # make gaussian gjf file
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    if 'user_fp_params' in jdata.keys() :        
        fp_params = jdata['user_fp_params']
    else:
        fp_params = jdata['fp_params']
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        sys_data = dpdata.System('POSCAR').data
        ret = make_gaussian_input(sys_data, fp_params)
        with open('input', 'w') as fp:
            fp.write(ret)
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
             jdata, 
             mdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        make_fp_vasp(iter_index, jdata) 
    elif fp_style == "pwscf" :
        make_fp_pwscf(iter_index, jdata) 
    elif fp_style == "gaussian" :
        make_fp_gaussian(iter_index, jdata)
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

def _gaussian_check_fin(ii):
    if os.path.isfile(os.path.join(ii, 'output')) :
        with open(os.path.join(ii, 'output'), 'r') as fp :
            content = fp.read()
            count = content.count('Normal termination of Gaussian')
            if count != 1 :
                return False
    else :
        return False
    return True
    
        
def run_fp_inner (iter_index,
                  jdata,
                  mdata,
                  ssh_sess,
                  forward_files,
                  backward_files,
                  check_fin,
                  log_file = "log",
                  forward_common_files=[]) :
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']
    machine_type = mdata['fp_machine']['machine_type']
    # fp_command = ("OMP_NUM_THREADS=1 mpirun -n %d " % fp_np) + fp_command
    # cpu task in parallel
    # if (('numb_gpu' not in fp_resources) or (fp_resources['numb_gpu'] == 0)) and (machine_type == 'slurm'):
    #     fp_command = 'srun ' + fp_command
    # if (('numb_gpu' not in fp_resources) or (fp_resources['numb_gpu'] == 0)) and (machine_type == 'pbs'):
    #     fp_command = 'mpirun  ' + fp_command

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

    if ssh_sess == None and machine_type == 'ucloud':
        ucloud_submit_jobs(mdata['fp_machine'],
                            mdata['fp_resources'],
                            fp_command,
                            work_path,
                            run_tasks,
                            fp_group_size,
                            [],
                            forward_files,
                            backward_files)
    elif machine_type == 'slurm' :        
        group_slurm_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           forward_common_files,
                           forward_files,
                           backward_files)
    elif machine_type == 'pbs' :        
        group_slurm_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files,
                          remote_job = PBSJob)
    elif machine_type == 'lsf' :        
        group_slurm_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files,
                          remote_job = LSFJob)
    elif machine_type == 'local' :        
        group_local_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files)
    else :
        raise RuntimeError("unknow machine type")

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
            mdata,
            ssh_sess) :
    fp_style = jdata['fp_style']
    fp_pp_files = jdata['fp_pp_files']

    if fp_style == "vasp" :
        forward_files = ['POSCAR', 'INCAR', 'KPOINTS'] + fp_pp_files 
        backward_files = ['OUTCAR','vasprun.xml']
        if mdata["fp_resources"]['cvasp']:
            forward_common_files=['cvasp.py']
        else:
            forward_common_files=[]
        run_fp_inner(iter_index, jdata, mdata, ssh_sess, forward_files, backward_files, _vasp_check_fin,
                     forward_common_files=forward_common_files) 
    elif fp_style == "pwscf" :
        forward_files = ['input'] + fp_pp_files
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, ssh_sess, forward_files, backward_files, _qe_check_fin, log_file = 'output') 
    elif fp_style == "gaussian":
        forward_files = ['input']
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, ssh_sess, forward_files, backward_files, _gaussian_check_fin, log_file = 'output') 
    else :
        raise RuntimeError ("unsupported fp style") 


def post_fp_vasp (iter_index,
                  jdata):
    model_devi_jobs = jdata['model_devi_jobs']
    assert (iter_index < len(model_devi_jobs)) 

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
    
    cwd = os.getcwd()
    for ss in system_index :
        sys_outcars = glob.glob(os.path.join(work_path, "task.%s.*/OUTCAR"%ss))
        sys_outcars.sort()                

        flag=True
        for oo in sys_outcars :
            if flag:
                try:
                    _sys = dpdata.LabeledSystem(oo)
                except:
                    try:
                       _sys = dpdata.LabeledSystem(oo.replace('OUTCAR','vasprun.xml'))
                    except:
                       _sys = dpdata.LabeledSystem()
                if len(_sys)>1:
                    dlog.info(oo + "has more than one systems in OUTCAR")
                    all_sys = _sys.sub_system([0])
                    flag = False
                elif len(_sys) == 1:
                    all_sys = _sys
                    flag = False
                else:
                    pass
            else:
                try:
                    _sys = dpdata.LabeledSystem(oo)
                except:
                    try:
                       _sys = dpdata.LabeledSystem(oo.replace('OUTCAR','vasprun.xml'))
                    except:
                       _sys = dpdata.LabeledSystem()
                if len(_sys)>1:
                    dlog.info(oo + "has more than one systems in OUTCAR")
                    all_sys.append(_sys.sub_system([0])) 
                elif len(_sys) == 1:
                    all_sys.append(_sys)

        #print("len(all_sys)",len(all_sys))
        sys_data_path = os.path.join(work_path, 'data.%s'%ss)
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_outcars))


def post_fp_pwscf (iter_index,
                   jdata):
    model_devi_jobs = jdata['model_devi_jobs']
    assert (iter_index < len(model_devi_jobs)) 

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
    
    cwd = os.getcwd()
    for ss in system_index :
        sys_output = glob.glob(os.path.join(work_path, "task.%s.*/output"%ss))
        sys_input = glob.glob(os.path.join(work_path, "task.%s.*/input"%ss))
        sys_output.sort()
        sys_input.sort()

        flag=True
        for ii,oo in zip(sys_input,sys_output) :
            if flag:
                _sys = dpdata.LabeledSystem()
                _sys.data=cvt_1frame(ii,oo)
                if len(_sys)>0:
                   all_sys=_sys
                   flag=False
                else:
                   pass
            else:
                _sys = dpdata.LabeledSystem()
                _sys.data = cvt_1frame(ii,oo)
                if len(_sys)>0:
                   all_sys.append(_sys)

        sys_data_path = os.path.join(work_path, 'data.%s'%ss)
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))


def post_fp_gaussian (iter_index,
                      jdata):
    model_devi_jobs = jdata['model_devi_jobs']
    assert (iter_index < len(model_devi_jobs)) 

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
    
    cwd = os.getcwd()
    for ss in system_index :
        sys_output = glob.glob(os.path.join(work_path, "task.%s.*/output"%ss))
        sys_output.sort()
        for idx,oo in enumerate(sys_output) :
            sys = dpdata.LabeledSystem(oo, fmt = 'gaussian/log') 
            if idx == 0:
                if 'use_clusters' in jdata and jdata['use_clusters']:
                    all_sys = dpdata.MultiSystems(sys)
                else:
                    all_sys = sys
            else:
                all_sys.append(sys)
        sys_data_path = os.path.join(work_path, 'data.%s'%ss)
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))


def post_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        post_fp_vasp(iter_index, jdata) 
    elif fp_style == "pwscf" :
        post_fp_pwscf(iter_index, jdata) 
    elif fp_style == 'gaussian' :
        post_fp_gaussian(iter_index, jdata)
    else :
        raise RuntimeError ("unsupported fp style")            
    
def run_iter (json_file, machine_file) :
    with open (json_file, 'r') as fp :
        jdata = json.load (fp)
    with open (machine_file, 'r') as fp:
        mdata = json.load (fp)

    max_tasks = 10000
    numb_task = 9
    record = "record.dpgen"

    train_machine = mdata['train_machine']    
    if ('machine_type' in train_machine) and  \
       (train_machine['machine_type'] == 'ucloud'):
        train_ssh_sess = None
    else :
        train_ssh_sess = SSHSession(train_machine)
    model_devi_machine = mdata['model_devi_machine']    
    if ('machine_type' in model_devi_machine) and  \
       (model_devi_machine['machine_type'] == 'ucloud'):
        model_devi_ssh_sess = None
    else :
        model_devi_ssh_sess = SSHSession(model_devi_machine)

    fp_machine = mdata['fp_machine']    
    if ('machine_type' in fp_machine) and  \
       (fp_machine['machine_type'] == 'ucloud'):
        fp_ssh_sess = None
    else :
        fp_ssh_sess = SSHSession(fp_machine)

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        dlog.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    cont = True
    ii = -1
    while cont:
        ii += 1
        iter_name=make_iter_name(ii)
        sepline(iter_name,'=')
        for jj in range (numb_task) :
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] : 
                continue
            task_name="task %02d"%jj
            sepline(task_name,'-')
            if   jj == 0 :
                log_iter ("make_train", ii, jj)
                make_train (ii, jdata, mdata) 
            elif jj == 1 :
                log_iter ("run_train", ii, jj)
                run_train  (ii, jdata, mdata, train_ssh_sess)
            elif jj == 2 :
                log_iter ("post_train", ii, jj)
                post_train (ii, jdata, mdata)
            elif jj == 3 :
                log_iter ("make_model_devi", ii, jj)
                cont = make_model_devi (ii, jdata, mdata)
                if not cont :
                    break
            elif jj == 4 :
                log_iter ("run_model_devi", ii, jj)
                run_model_devi (ii, jdata, mdata, model_devi_ssh_sess)
            elif jj == 5 :
                log_iter ("post_model_devi", ii, jj)
                post_model_devi (ii, jdata, mdata)
            elif jj == 6 :
                log_iter ("make_fp", ii, jj)
                make_fp (ii, jdata, mdata)
            elif jj == 7 :
                log_iter ("run_fp", ii, jj)
                run_fp (ii, jdata, mdata, fp_ssh_sess)
            elif jj == 8 :
                log_iter ("post_fp", ii, jj)
                post_fp (ii, jdata)
            else :
                raise RuntimeError ("unknow task %d, something wrong" % jj)
            record_iter (record, ii, jj)
def gen_run(args) :
    if args.PARAM and args.MACHINE:
       dlog.info ("start running")
       run_iter (args.PARAM, args.MACHINE)
       dlog.info ("finished")    

def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("PARAM", type=str, 
                        help="The parameters of the generator")
    parser.add_argument("MACHINE", type=str, 
                        help="The settings of the machine running the generator")
    args = parser.parse_args()

    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')
    # logging.basicConfig (filename="compute_string.log", filemode="a", level=logging.INFO, format='%(asctime)s %(message)s')
    logging.getLogger("paramiko").setLevel(logging.WARNING)

    logging.info ("start running")
    run_iter (args.PARAM, args.MACHINE)
    logging.info ("finished!")

if __name__ == '__main__':
    _main()

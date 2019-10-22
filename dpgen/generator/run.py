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
import scipy.constants as pc
from collections import Counter
from distutils.version import LooseVersion
from dpgen import dlog
from dpgen import SHORT_CMD
from dpgen.generator.lib.utils import make_iter_name
from dpgen.generator.lib.utils import create_path
from dpgen.generator.lib.utils import copy_file_list
from dpgen.generator.lib.utils import replace
from dpgen.generator.lib.utils import log_iter
from dpgen.generator.lib.utils import record_iter
from dpgen.generator.lib.utils import log_task
from dpgen.generator.lib.lammps import make_lammps_input
from dpgen.generator.lib.vasp import write_incar_dict
from dpgen.generator.lib.vasp import make_vasp_incar_user_dict
from dpgen.generator.lib.vasp import incar_upper
from dpgen.generator.lib.pwscf import make_pwscf_input
#from dpgen.generator.lib.pwscf import cvt_1frame
from dpgen.generator.lib.siesta import make_siesta_input
from dpgen.generator.lib.gaussian import make_gaussian_input, take_cluster
from dpgen.generator.lib.cp2k import make_cp2k_input, make_cp2k_xyz
from dpgen.generator.lib.ele_temp import NBandsEsti
from dpgen.remote.RemoteJob import SSHSession, JobStatus, SlurmJob, PBSJob, LSFJob, CloudMachineJob, awsMachineJob
from dpgen.remote.group_jobs import ucloud_submit_jobs, aws_submit_jobs
from dpgen.remote.group_jobs import group_slurm_jobs
from dpgen.remote.group_jobs import group_local_jobs
from dpgen.remote.decide_machine import decide_train_machine, decide_fp_machine, decide_model_devi_machine
from dpgen.dispatcher.Dispatcher import Dispatcher, make_dispatcher
from dpgen.util import sepline
from dpgen import ROOT_PATH
from pymatgen.io.vasp import Incar,Kpoints,Potcar
from dpgen.auto_test.lib.vasp import make_kspacing_kpoints

template_name = 'template'
train_name = '00.train'
train_task_fmt = '%03d'
train_tmpl_path = os.path.join(template_name, train_name)
default_train_input_file = 'input.json'
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
    # train_param = jdata['train_param']
    train_input_file = default_train_input_file
    numb_models = jdata['numb_models']
    init_data_prefix = jdata['init_data_prefix']
    init_data_prefix = os.path.abspath(init_data_prefix)
    init_data_sys_ = jdata['init_data_sys']
    fp_task_min = jdata['fp_task_min']
    model_devi_jobs = jdata['model_devi_jobs']
    use_ele_temp = jdata.get('use_ele_temp', 0)    

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
    if 'init_batch_size' in jdata:
        init_batch_size_ = list(jdata['init_batch_size'])
    else:
        init_batch_size_ = ["auto" for aa in range(len(jdata['init_data_sys']))]
    if 'sys_batch_size' in jdata:
        sys_batch_size = jdata['sys_batch_size']
    else:
        sys_batch_size = ["auto" for aa in range(len(jdata['sys_configs']))]
    for ii, ss in zip(init_data_sys_, init_batch_size_) :
        if jdata.get('init_multi_systems', False):
            for single_sys in os.listdir(os.path.join(work_path, 'data.init', ii)):
                init_data_sys.append(os.path.join('..', 'data.init', ii, single_sys))
                init_batch_size.append(detect_batch_size(ss, os.path.join(work_path, 'data.init', ii, single_sys)))
        else:
            init_data_sys.append(os.path.join('..', 'data.init', ii))
            init_batch_size.append(detect_batch_size(ss, os.path.join(work_path, 'data.init', ii)))
    if iter_index > 0 :
        for ii in range(iter_index) :
            fp_path = os.path.join(make_iter_name(ii), fp_name)
            fp_data_sys = glob.glob(os.path.join(fp_path, "data.*"))
            for jj in fp_data_sys :
                sys_idx = int(jj.split('.')[-1])
                if jdata.get('use_clusters', False):
                    nframes = 0
                    for sys_single in os.listdir(jj):
                        tmp_box = np.loadtxt(os.path.join(jj, sys_single, 'box.raw'))
                        tmp_box = np.reshape(tmp_box, [-1,9])
                        nframes += tmp_box.shape[0]
                    if nframes < fp_task_min :
                        log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                        continue
                    for sys_single in os.listdir(jj):
                        init_data_sys.append(os.path.join('..', 'data.iters', jj, sys_single))
                        init_batch_size.append(detect_batch_size(sys_batch_size[sys_idx], os.path.join(jj, sys_single)))
                else:
                    tmp_box = np.loadtxt(os.path.join(jj, 'box.raw'))
                    tmp_box = np.reshape(tmp_box, [-1,9])
                    nframes = tmp_box.shape[0]
                    if nframes < fp_task_min :
                        log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                        continue
                    init_data_sys.append(os.path.join('..', 'data.iters', jj))
                    init_batch_size.append(detect_batch_size(sys_batch_size[sys_idx], jj))
    # establish tasks
    jinput = jdata['default_training_param']
    try:
        mdata["deepmd_version"]
    except:
        mdata = set_version(mdata)
    if LooseVersion(mdata["deepmd_version"]) < LooseVersion('1'):
        # 0.x
        jinput['systems'] = init_data_sys
        jinput['batch_size'] = init_batch_size
        if use_ele_temp:
            raise RuntimeError('the electron temperature is only supported by deepmd-kit >= 1.0.0, please upgrade your deepmd-kit')
    else:
        # 1.x
        jinput['training']['systems'] = init_data_sys
        jinput['training']['batch_size'] = init_batch_size
        # electron temperature
        if use_ele_temp == 0:
            pass
        elif use_ele_temp == 1:
            jinput['model']['fitting_net']['numb_fparam'] = 1
            jinput['model']['fitting_net'].pop('numb_aparam', None)
        elif use_ele_temp == 2:
            jinput['model']['fitting_net']['numb_aparam'] = 1
            jinput['model']['fitting_net'].pop('numb_fparam', None)
        else:
            raise RuntimeError('invalid setting for use_ele_temp ' + use_ele_temp)
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(task_path)
        os.chdir(task_path)
        for ii in init_data_sys :
            if not os.path.isdir(ii) :
                raise RuntimeError ("data sys %s does not exists, cwd is %s" % (ii, os.getcwd()))
        os.chdir(cwd)
        if LooseVersion(mdata["deepmd_version"]) < LooseVersion('1'):
            # 0.x
            jinput['seed'] = random.randrange(sys.maxsize) % (2**32)
        else:
            # 1.x
            jinput['model']['descriptor']['seed'] = random.randrange(sys.maxsize) % (2**32)
            jinput['model']['fitting_net']['seed'] = random.randrange(sys.maxsize) % (2**32)
            jinput['training']['seed'] = random.randrange(sys.maxsize) % (2**32)
        with open(os.path.join(task_path, train_input_file), 'w') as outfile:
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

def detect_batch_size(batch_size, system=None):
    if type(batch_size) == int:
        return batch_size
    elif batch_size == "auto":
        # automaticcaly set batch size, batch_size = 32 // atom_numb (>=1, <=fram_numb)
        s = dpdata.LabeledSystem(system, fmt='deepmd/npy')
        return int(min( np.ceil(32.0 / float(s["coords"].shape[1]) ), s["coords"].shape[0]))
    else:
        raise RuntimeError("Unsupported batch size")

def run_train (iter_index,
               jdata,
               mdata,
               dispatcher) :
    # load json param
    numb_models = jdata['numb_models']
    # train_param = jdata['train_param']
    train_input_file = default_train_input_file
    try:
        mdata["deepmd_version"]
    except:
        mdata = set_version(mdata)
    if LooseVersion(mdata["deepmd_version"]) < LooseVersion('1'):
        # 0.x
        deepmd_path = mdata['deepmd_path']
    else:
        # 1.x
        python_path = mdata['python_path']
    train_resources = mdata['train_resources']

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
    commands = []
    if LooseVersion(mdata["deepmd_version"]) < LooseVersion('1'):
        # 0.x
        command = os.path.join(deepmd_path, 'bin/dp_train %s' % train_input_file)
        commands.append(command)
        command = os.path.join(deepmd_path, 'bin/dp_frz')
        commands.append(command)        
    else:
        # 1.x
        command =  '%s -m deepmd train %s' % (python_path, train_input_file)
        commands.append(command)
        command = '%s -m deepmd freeze' % python_path
        commands.append(command)

    #_tasks = [os.path.basename(ii) for ii in all_task]
    # run_tasks = []
    # for ii in all_task:
    #     check_pb = os.path.join(ii, "frozen_model.pb")
    #     check_lcurve = os.path.join(ii, "lcurve.out")
    #     if os.path.isfile(check_pb) and os.path.isfile(check_lcurve):
    #         pass
    #     else:
    #         run_tasks.append(ii)
    run_tasks = [os.path.basename(ii) for ii in all_task]

    forward_files = [train_input_file]
    backward_files = ['frozen_model.pb', 'lcurve.out', 'train.log']
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
        if jdata.get('init_multi_systems', False):
            for single_sys in os.listdir(os.path.join(ii)):
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'type.raw'))
        else:
            trans_comm_data += glob.glob(os.path.join(ii, 'set.*'))
            trans_comm_data += glob.glob(os.path.join(ii, 'type.raw'))
    for ii in fp_data :
        if jdata.get('use_clusters', False):
            for single_sys in os.listdir(os.path.join(ii)):
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'type.raw'))
        else:
            trans_comm_data += glob.glob(os.path.join(ii, 'set.*'))
            trans_comm_data += glob.glob(os.path.join(ii, 'type.raw'))
    os.chdir(cwd)

    try:
        train_group_size = mdata['train_group_size']
    except:
        train_group_size = 1

    dispatcher.run_jobs(mdata['train_resources'],
                        commands,
                        work_path,
                        run_tasks,
                        train_group_size,
                        trans_comm_data,
                        forward_files,
                        backward_files,
                        outlog = 'train.log',
                        errlog = 'train.log')


def post_train (iter_index,
                jdata,
                mdata) :
    # load json param
    numb_models = jdata['numb_models']
    # paths
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    # check if is copied
    copy_flag = os.path.join(work_path, 'copied')
    if os.path.isfile(copy_flag) :
        log_task('copied model, do not post train')
        return
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
    use_ele_temp = jdata.get('use_ele_temp', 0)
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
            system = dpdata.System(os.path.join(conf_path, poscar_name), fmt = fmt, type_map = jdata['type_map'])
            system.to_lammps_lmp(os.path.join(conf_path, lmp_name))
            conf_counter += 1
        sys_counter += 1

    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        task_counter = 0
        for cc in ss :
            for tt_ in temps:
                if use_ele_temp:
                    if type(tt_) == list:
                        tt = tt_[0]
                        if use_ele_temp == 1:
                            te_f = tt_[1]
                            te_a = None
                        else:
                            te_f = None
                            te_a = tt_[1]
                    else:
                        assert(type(tt_) == float or type(tt_) == int)
                        tt = float(tt_)
                        if use_ele_temp == 1:
                            te_f = tt
                            te_a = None
                        else:
                            te_f = None
                            te_a = tt
                else :
                    tt = tt_
                    te_f = None
                    te_a = None
                for pp in press:
                    task_name = make_model_devi_task_name(sys_idx[sys_counter], task_counter)
                    conf_name = make_model_devi_conf_name(sys_idx[sys_counter], conf_counter) + '.lmp'
                    task_path = os.path.join(work_path, task_name)
                    # dlog.info(task_path)
                    create_path(task_path)
                    create_path(os.path.join(task_path, 'traj'))
                    loc_conf_name = 'conf.lmp'
                    os.symlink(os.path.join(os.path.join('..','confs'), conf_name),
                               os.path.join(task_path, loc_conf_name) )
                    cwd_ = os.getcwd()
                    os.chdir(task_path)
                    try:
                        mdata["deepmd_version"]
                    except:
                        mdata = set_version(mdata)
                    deepmd_version = mdata['deepmd_version']
                    file_c = make_lammps_input(ensemble,
                                               loc_conf_name,
                                               task_model_list,
                                               nsteps,
                                               model_devi_dt,
                                               model_devi_neidelay,
                                               trj_freq,
                                               mass_map,
                                               tt,
                                               jdata = jdata,
                                               tau_t = model_devi_taut,
                                               pres = pp,
                                               tau_p = model_devi_taup,
                                               pka_e = pka_e,
                                               ele_temp_f = te_f,
                                               ele_temp_a = te_a,
                                               deepmd_version = deepmd_version)
                    job = {}
                    job["ensemble"] = ensemble
                    job["press"] = pp
                    job["temps"] = tt
                    if te_f is not None:
                        job["ele_temp"] = te_f
                    if te_a is not None:
                        job["ele_temp"] = te_a
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
                    dispatcher) :
    #rmdlog.info("This module has been run !")
    lmp_exec = mdata['lmp_command']
    model_devi_group_size = mdata['model_devi_group_size']
    model_devi_resources = mdata['model_devi_resources']

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))

    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = lmp_exec + " -i input.lammps"
    commands = [command]

    fp = open (os.path.join(work_path, 'cur_job.json'), 'r')
    cur_job = json.load (fp)
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt = parse_cur_job(cur_job)
    nframes = nsteps // trj_freq + 1
    
    run_tasks_ = all_task
    # for ii in all_task:
    #     fres = os.path.join(ii, 'model_devi.out')
    #     if os.path.isfile(fres) :
    #         nlines = np.loadtxt(fres).shape[0]
    #         if nframes != nlines :
    #             run_tasks_.append(ii)
    #     else :
    #         run_tasks_.append(ii)

    run_tasks = [os.path.basename(ii) for ii in run_tasks_]
    #dlog.info("all_task is ", all_task)
    #dlog.info("run_tasks in run_model_deviation",run_tasks_)
    all_models = glob.glob(os.path.join(work_path, 'graph*pb'))
    model_names = [os.path.basename(ii) for ii in all_models]
    forward_files = ['conf.lmp', 'input.lammps', 'traj']
    backward_files = ['model_devi.out', 'model_devi.log', 'traj']

    dispatcher.run_jobs(mdata['model_devi_resources'],
                        commands,
                        work_path,
                        run_tasks,
                        model_devi_group_size,
                        model_names,
                        forward_files,
                        backward_files,
                        outlog = 'model_devi.log',
                        errlog = 'model_devi.log')


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
                         jdata):
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
    cluster_cutoff = jdata['cluster_cutoff'] if jdata.get('use_clusters', False) else None
    # skip save *.out if detailed_report_make_fp is False, default is True
    detailed_report_make_fp = jdata.get("detailed_report_make_fp", True)
    for ss in system_index :
        fp_candidate = []
        if detailed_report_make_fp:
            fp_rest_accurate = []
            fp_rest_failed = []
        modd_system_glob = os.path.join(modd_path, 'task.' + ss + '.*')
        modd_system_task = glob.glob(modd_system_glob)
        modd_system_task.sort()
        cc = 0
        counter = Counter()
        for tt in modd_system_task :
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                all_conf = np.loadtxt(os.path.join(tt, 'model_devi.out'))
                for ii in range(all_conf.shape[0]) :
                    if all_conf[ii][0] < model_devi_skip :
                        continue
                    cc = int(all_conf[ii][0])
                    if cluster_cutoff is None:
                        if (all_conf[ii][1] < e_trust_hi and all_conf[ii][1] >= e_trust_lo) or \
                           (all_conf[ii][4] < f_trust_hi and all_conf[ii][4] >= f_trust_lo) :
                            fp_candidate.append([tt, cc])
                            counter['candidate'] += 1
                        elif (all_conf[ii][1] >= e_trust_hi ) or (all_conf[ii][4] >= f_trust_hi ):
                            if detailed_report_make_fp:
                                fp_rest_failed.append([tt, cc])
                            counter['failed'] += 1
                        elif (all_conf[ii][1] < e_trust_lo and all_conf[ii][4] < f_trust_lo ):
                            if detailed_report_make_fp:
                                fp_rest_accurate.append([tt, cc])
                            counter['accurate'] += 1
                        else :
                            raise RuntimeError('md traj %s frame %d with f devi %f does not belong to either accurate, candidiate and failed, it should not happen' % (tt, ii, all_conf[ii][4]))
                    else:
                        idx_candidate = np.where(np.logical_and(all_conf[ii][7:] < f_trust_hi, all_conf[ii][7:] >= f_trust_lo))[0]
                        for jj in idx_candidate:
                            fp_candidate.append([tt, cc, jj])
                        counter['candidate'] += len(idx_candidate)
                        idx_rest_accurate = np.where(all_conf[ii][7:] < f_trust_lo)[0]
                        if detailed_report_make_fp:
                            for jj in idx_rest_accurate:
                                fp_rest_accurate.append([tt, cc, jj])
                        counter['accurate'] += len(idx_rest_accurate)
                        idx_rest_failed = np.where(all_conf[ii][7:] >= f_trust_hi)[0]
                        if detailed_report_make_fp:
                            for jj in idx_rest_failed:
                                fp_rest_failed.append([tt, cc, jj])
                        counter['failed'] += len(idx_rest_failed)
        # print a report
        fp_sum = sum(counter.values())
        for cc_key, cc_value in counter.items():
            dlog.info("system {0:s} {1:9s} : {2:6d} in {3:6d} {4:6.2f} %".format(ss, cc_key, cc_value, fp_sum, cc_value/fp_sum*100))
        random.shuffle(fp_candidate)
        if detailed_report_make_fp:
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
                new_system = take_cluster(conf_name, type_map, jj, jdata)
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
                np.save("atom_pref", new_system.data["atom_pref"])
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
    
def make_vasp_incar(jdata, filename):
    if 'fp_incar' in jdata.keys() :
        fp_incar_path = jdata['fp_incar']
        assert(os.path.exists(fp_incar_path))
        fp_incar_path = os.path.abspath(fp_incar_path)
        fr = open(fp_incar_path)
        incar = fr.read()
        fr.close()
    elif 'user_fp_params' in jdata.keys() :
        incar = write_incar_dict(jdata['user_fp_params'])
    else:
        incar = make_vasp_incar_user_dict(jdata['fp_params'])
    with open(filename, 'w') as fp:
        fp.write(incar)
    return incar    

def make_vasp_incar_ele_temp(jdata, filename, ele_temp, nbands_esti = None):
    with open(filename) as fp:
        incar = fp.read()
    incar = incar_upper(Incar.from_string(incar))
    incar['ISMEAR'] = -1
    incar['SIGMA'] = ele_temp * pc.Boltzmann / pc.electron_volt
    incar.write_file('INCAR')
    if nbands_esti is not None:
        nbands = nbands_esti.predict('.')
        with open(filename) as fp:
            incar = Incar.from_string(fp.read())
        incar['NBANDS'] = nbands
        incar.write_file('INCAR')

def _make_fp_vasp_incar (iter_index,
                         jdata,
                         nbands_esti = None) :
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        make_vasp_incar(jdata, 'INCAR')
        if os.path.exists('job.json'):
            with open('job.json') as fp:
                job_data = json.load(fp)
            if 'ele_temp' in job_data:
                make_vasp_incar_ele_temp(jdata, 'INCAR', 
                                         job_data['ele_temp'],
                                         nbands_esti = nbands_esti)
        os.chdir(cwd)

def _make_fp_vasp_kp (iter_index,jdata):
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        # get kspacing and kgamma from incar
        assert(os.path.exists('INCAR'))
        with open('INCAR') as fp:
            incar = fp.read()
        standard_incar = incar_upper(Incar.from_string(incar))
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
        # check poscar
        assert(os.path.exists('POSCAR'))
        # make kpoints
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

def sys_link_fp_vasp_pp (iter_index,
                         jdata) :
    fp_pp_path = jdata['fp_pp_path']
    fp_pp_files = jdata['fp_pp_files']
    fp_pp_path = os.path.abspath(fp_pp_path)
    type_map = jdata['type_map']
    assert(os.path.exists(fp_pp_path))
    assert(len(fp_pp_files) == len(type_map)), 'size of fp_pp_files should be the same as the size of type_map'

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    system_idx_str = [os.path.basename(ii).split('.')[1] for ii in fp_tasks]
    system_idx_str = list(set(system_idx_str))
    system_idx_str.sort()
    for ii in system_idx_str:
        potcars = []
        sys_tasks = glob.glob(os.path.join(work_path, 'task.%s.*' % ii))
        assert (len(sys_tasks) != 0)
        sys_poscar = os.path.join(sys_tasks[0], 'POSCAR')
        sys = dpdata.System(sys_poscar, fmt = 'vasp/poscar')
        for ele_name in sys['atom_names']:
            ele_idx = jdata['type_map'].index(ele_name)
            potcars.append(fp_pp_files[ele_idx])                
        with open(os.path.join(work_path,'POTCAR.%s' % ii), 'w') as fp_pot:
            for jj in potcars:
                with open(os.path.join(fp_pp_path, jj)) as fp:
                    fp_pot.write(fp.read())
        sys_tasks = glob.glob(os.path.join(work_path, 'task.%s.*' % ii))
        cwd = os.getcwd()
        for jj in sys_tasks:
            os.chdir(jj)
            os.symlink(os.path.join('..', 'POTCAR.%s' % ii), 'POTCAR')
            os.chdir(cwd)

def _make_fp_vasp_configs(iter_index,
                          jdata):
    fp_task_max = jdata['fp_task_max']
    model_devi_skip = jdata['model_devi_skip']
    e_trust_lo = 1e+10
    e_trust_hi = 1e+10
    f_trust_lo = jdata['model_devi_f_trust_lo']
    f_trust_hi = jdata['model_devi_f_trust_hi']
    type_map = jdata['type_map']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)

    #copy cvasp.py
    # Move cvasp interface to jdata
    if ('cvasp' in jdata) and (jdata['cvasp'] == True):
        shutil.copyfile(cvasp_file, os.path.join(work_path,'cvasp.py'))

    modd_path = os.path.join(iter_name, model_devi_name)
    task_min = -1
    if os.path.isfile(os.path.join(modd_path, 'cur_job.json')) :
        cur_job = json.load(open(os.path.join(modd_path, 'cur_job.json'), 'r'))
        if 'task_min' in cur_job :
            task_min = cur_job['task_min']
    # make configs
    fp_tasks = _make_fp_vasp_inner(modd_path, work_path,
                                   model_devi_skip,
                                   e_trust_lo, e_trust_hi,
                                   f_trust_lo, f_trust_hi,
                                   task_min, fp_task_max,
                                   [],
                                   type_map,
                                   jdata)
    return fp_tasks

def make_fp_vasp (iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return
    # abs path for fp_incar if it exists
    if 'fp_incar' in jdata:
        jdata['fp_incar'] = os.path.abspath(jdata['fp_incar'])
    # get nbands esti if it exists
    if 'fp_nbands_esti_data' in jdata:
        nbe = NBandsEsti(jdata['fp_nbands_esti_data'])
    else:
        nbe = None
    # order is critical!
    # 1, create potcar
    sys_link_fp_vasp_pp(iter_index, jdata)
    # 2, create incar
    _make_fp_vasp_incar(iter_index, jdata, nbands_esti = nbe)
    # 3, create kpoints
    _make_fp_vasp_kp(iter_index, jdata)


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


def make_fp_siesta(iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return
    # make siesta input
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
        ret = make_siesta_input(sys_data, fp_pp_files, fp_params)
        with open('input', 'w') as fp:
            fp.write(ret)
        os.chdir(cwd)
    # link pp files
    _link_fp_vasp_pp(iter_index, jdata)
        
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

def make_fp_cp2k (iter_index,
                  jdata):
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return
    # make cp2k input
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
        # make input for every task
        cp2k_input = make_cp2k_input(sys_data, fp_params)
        with open('input.inp', 'w') as fp:
            fp.write(cp2k_input)
            fp.close()
        # make coord.xyz used by cp2k for every task
        cp2k_coord = make_cp2k_xyz(sys_data)
        with open('coord.xyz', 'w') as fp:
            fp.write(cp2k_coord)
            fp.close()
        os.chdir(cwd)

    # link pp files
    _link_fp_vasp_pp(iter_index, jdata)

def make_fp (iter_index,
             jdata,
             mdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        make_fp_vasp(iter_index, jdata)
    elif fp_style == "pwscf" :
        make_fp_pwscf(iter_index, jdata)
    elif fp_style == "siesta" :
        make_fp_siesta(iter_index, jdata)
    elif fp_style == "gaussian" :
        make_fp_gaussian(iter_index, jdata)
    elif fp_style == "cp2k" :
        make_fp_cp2k(iter_index, jdata)
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


def _siesta_check_fin(ii) :
    if os.path.isfile(os.path.join(ii, 'output')) :
        with open(os.path.join(ii, 'output'), 'r') as fp :
            content = fp.read()
            count = content.count('End of run')
            if count != 1 :
                return False
    else :
        return False
    return True

def _gaussian_check_fin(ii):
    if os.path.isfile(os.path.join(ii, 'output')) :
        with open(os.path.join(ii, 'output'), 'r') as fp :
            content = fp.read()
            count = content.count('termination')
            if count == 0 :
                return False
    else :
        return False
    return True

def _cp2k_check_fin(ii):
    if os.path.isfile(os.path.join(ii, 'output')) :
        with open(os.path.join(ii, 'output'), 'r') as fp :
            content = fp.read()
            count = content.count('SCF run converged')
            if count == 0 :
                return False
    else :
        return False
    return True

def run_fp_inner (iter_index,
                  jdata,
                  mdata,
                  dispatcher,
                  forward_files,
                  backward_files,
                  check_fin,
                  log_file = "log",
                  forward_common_files=[]) :
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    fp_run_tasks = fp_tasks
    # for ii in fp_tasks :
    #     if not check_fin(ii) :
    #         fp_run_tasks.append(ii)
    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]

    dispatcher.run_jobs(mdata['fp_resources'],
                        [fp_command],
                        work_path,
                        run_tasks,
                        fp_group_size,
                        forward_common_files,
                        forward_files,
                        backward_files,
                        outlog = log_file,
                        errlog = log_file)


def run_fp (iter_index,
            jdata,
            mdata,
            dispatcher) :
    fp_style = jdata['fp_style']
    fp_pp_files = jdata['fp_pp_files']

    if fp_style == "vasp" :
        forward_files = ['POSCAR', 'INCAR', 'POTCAR']
        backward_files = ['OUTCAR','vasprun.xml']
        # Move cvasp interface to jdata
        if ('cvasp' in jdata) and (jdata['cvasp'] == True):
            mdata['fp_resources']['cvasp'] = True
        if ('cvasp' in  mdata["fp_resources"] ) and (mdata["fp_resources"]["cvasp"]==True):
            #dlog.info("cvasp is on !")
            forward_common_files=['cvasp.py']
            forward_files.append('KPOINTS')
        else:
            forward_common_files=[]
        run_fp_inner(iter_index, jdata, mdata, dispatcher, forward_files, backward_files, _vasp_check_fin,
                     forward_common_files=forward_common_files)
    elif fp_style == "pwscf" :
        forward_files = ['input'] + fp_pp_files
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, dispatcher, forward_files, backward_files, _qe_check_fin, log_file = 'output')
    elif fp_style == "siesta":
        forward_files = ['input'] + fp_pp_files
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, dispatcher, forward_files, backward_files, _siesta_check_fin, log_file='output')
    elif fp_style == "gaussian":
        forward_files = ['input']
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, dispatcher, forward_files, backward_files, _gaussian_check_fin, log_file = 'output')
    elif fp_style == "cp2k":
        forward_files = ['input.inp', 'coord.xyz']
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, dispatcher, forward_files, backward_files, _cp2k_check_fin, log_file = 'output')
    else :
        raise RuntimeError ("unsupported fp style")


def post_fp_vasp (iter_index,
                  jdata,
                  rfailed=None):

    ratio_failed =  rfailed if rfailed else jdata.get('ratio_failed',0.05)
    model_devi_jobs = jdata['model_devi_jobs']
    assert (iter_index < len(model_devi_jobs))
    use_ele_temp = jdata.get('use_ele_temp', 0)

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

    tcount=0
    icount=0
    for ss in system_index :
        sys_outcars = glob.glob(os.path.join(work_path, "task.%s.*/OUTCAR"%ss))
        sys_outcars.sort()
        tcount += len(sys_outcars)
        all_sys = None
        all_te = []
        for oo in sys_outcars :
            try:
                _sys = dpdata.LabeledSystem(oo, type_map = jdata['type_map'])
            except:
                dlog.info('Try to parse from vasprun.xml')
                try:
                   _sys = dpdata.LabeledSystem(oo.replace('OUTCAR','vasprun.xml'), type_map = jdata['type_map'])
                except:
                   _sys = dpdata.LabeledSystem()
                   dlog.info('Failed fp path: %s'%oo.replace('OUTCAR',''))
            if len(_sys) == 1:
                if all_sys is None:
                    all_sys = _sys
                else:
                    all_sys.append(_sys)
                # save ele_temp, if any
                with open(oo.replace('OUTCAR', 'job.json')) as fp:
                    job_data = json.load(fp)
                if 'ele_temp' in job_data:
                    assert(use_ele_temp)
                    ele_temp = job_data['ele_temp']
                    all_te.append(ele_temp)
            else:
                icount+=1
        all_te = np.array(all_te)
        if all_sys is not None:
           sys_data_path = os.path.join(work_path, 'data.%s'%ss)
           all_sys.to_deepmd_raw(sys_data_path)
           all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_outcars))
           if all_te.size > 0:
               assert(len(all_sys) == all_sys.get_nframes())
               assert(len(all_sys) == all_te.size)
               all_te = np.reshape(all_te, [-1,1])
               if use_ele_temp == 0:
                   raise RuntimeError('should not get ele temp at setting: use_ele_temp == 0')
               elif use_ele_temp == 1:
                   np.savetxt(os.path.join(sys_data_path, 'fparam.raw'), all_te)
                   np.save(os.path.join(sys_data_path, 'set.000', 'fparam.npy'), all_te)
               elif use_ele_temp == 2:
                   tile_te = np.tile(all_te, [1, all_sys.get_natoms()])
                   np.savetxt(os.path.join(sys_data_path, 'aparam.raw'), tile_te)
                   np.save(os.path.join(sys_data_path, 'set.000', 'aparam.npy'), tile_te)
               else:
                   raise RuntimeError('invalid setting of use_ele_temp ' + use_ele_temp)

    dlog.info("failed frame number: %s "%icount)
    dlog.info("total frame number: %s "%tcount)
    reff=icount/tcount
    dlog.info('ratio of failed frame:  {:.2%}'.format(reff))

    if reff>ratio_failed:
       raise RuntimeError("find too many unsuccessfully terminated jobs")


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
                _sys = dpdata.LabeledSystem(oo, fmt = 'qe/pw/scf', type_map = jdata['type_map'])
                if len(_sys)>0:
                   all_sys=_sys
                   flag=False
                else:
                   pass
            else:
                _sys = dpdata.LabeledSystem(oo, fmt = 'qe/pw/scf', type_map = jdata['type_map'])
                if len(_sys)>0:
                   all_sys.append(_sys)

        sys_data_path = os.path.join(work_path, 'data.%s'%ss)
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))

def post_fp_siesta (iter_index,
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
        for idx, oo in enumerate(sys_output):
            _sys = dpdata.LabeledSystem()
            _sys.data['atom_names'], \
            _sys.data['atom_numbs'], \
            _sys.data['atom_types'], \
            _sys.data['cells'], \
            _sys.data['coords'], \
            _sys.data['energies'], \
            _sys.data['forces'], \
            _sys.data['virials'] \
            = dpdata.siesta.output.obtain_frame(oo)
            if idx == 0:
                all_sys = _sys
            else:
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
            if len(sys) > 0:
                sys.check_type_map(type_map = jdata['type_map'])
            if jdata.get('use_atom_pref', False):
                sys.data['atom_pref'] = np.load(os.path.join(os.path.dirname(oo), "atom_pref.npy"))
            if idx == 0:
                if jdata.get('use_clusters', False):
                    all_sys = dpdata.MultiSystems(sys, type_map = jdata['type_map'])
                else:
                    all_sys = sys
            else:
                all_sys.append(sys)
        sys_data_path = os.path.join(work_path, 'data.%s'%ss)
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))


def post_fp_cp2k (iter_index,
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
            sys = dpdata.LabeledSystem(oo, fmt = 'cp2k/output')
            if len(sys) > 0:
                sys.check_type_map(type_map = jdata['type_map'])
            if idx == 0:
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
    elif fp_style == "siesta":
        post_fp_siesta(iter_index, jdata)
    elif fp_style == 'gaussian' :
        post_fp_gaussian(iter_index, jdata)
    elif fp_style == 'cp2k' :
        post_fp_cp2k(iter_index, jdata)
    else :
        raise RuntimeError ("unsupported fp style")
    # clean traj
    iter_name = make_iter_name(iter_index)
    clean_traj = True
    if 'model_devi_clean_traj' in jdata :
        clean_traj = jdata['model_devi_clean_traj']
    if clean_traj:
        modd_path = os.path.join(iter_name, model_devi_name)
        md_trajs = glob.glob(os.path.join(modd_path, 'task*/traj'))
        for ii in md_trajs :
            shutil.rmtree(ii)

def set_version(mdata):
    if 'deepmd_path' in mdata:
        deepmd_version = '0.1'
    elif 'python_path' in mdata:
        deepmd_version = '1'
    elif 'train' in mdata:
        if 'deepmd_path' in mdata['train'][0]:
            deepmd_version = '0.1'
        elif 'python_path' in mdata['train'][0]:
            deepmd_version = '1'
        else:
            deepmd_version = '0.1'
    else:
        deepmd_version = '0.1'
    # set
    mdata['deepmd_version'] = deepmd_version
    return mdata



    

def run_iter (param_file, machine_file) :
    try:
       import ruamel
       from monty.serialization import loadfn,dumpfn
       warnings.simplefilter('ignore', ruamel.yaml.error.MantissaNoDotYAML1_1Warning)
       jdata=loadfn(param_file)
       mdata=loadfn(machine_file)
    except:
       with open (param_file, 'r') as fp :
           jdata = json.load (fp)
       with open (machine_file, 'r') as fp:
           mdata = json.load (fp)

    if jdata.get('pretty_print',False):
       #assert(jdata["pretty_format"] in ['json','yaml'])
       fparam=SHORT_CMD+'_'+param_file.split('.')[0]+'.'+jdata.get('pretty_format','json')
       dumpfn(jdata,fparam,indent=4)
       fmachine=SHORT_CMD+'_'+machine_file.split('.')[0]+'.'+jdata.get('pretty_format','json')
       dumpfn(mdata,fmachine,indent=4)

    max_tasks = 10000
    numb_task = 9
    record = "record.dpgen"
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
                mdata  = decide_train_machine(mdata)
                disp = make_dispatcher(mdata['train_machine'])
                run_train  (ii, jdata, mdata, disp)
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
                mdata = decide_model_devi_machine(mdata)
                disp = make_dispatcher(mdata['model_devi_machine'])
                run_model_devi (ii, jdata, mdata, disp)
            elif jj == 5 :
                log_iter ("post_model_devi", ii, jj)
                post_model_devi (ii, jdata, mdata)
            elif jj == 6 :
                log_iter ("make_fp", ii, jj)
                make_fp (ii, jdata, mdata)
            elif jj == 7 :
                log_iter ("run_fp", ii, jj)
                mdata = decide_fp_machine(mdata)
                disp = make_dispatcher(mdata['fp_machine'])
                run_fp (ii, jdata, mdata, disp)
            elif jj == 8 :
                log_iter ("post_fp", ii, jj)
                post_fp (ii, jdata)
            else :
                raise RuntimeError ("unknown task %d, something wrong" % jj)
            record_iter (record, ii, jj)


def gen_run(args) :
    if args.PARAM and args.MACHINE:
        if args.debug:
            dlog.setLevel(logging.DEBUG)
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

#!/usr/bin/env python3

"""
init: data
iter:
        00.train
        01.model_devi
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
import logging.handlers
import queue
import warnings
import shutil
import time
import copy
import dpdata
import numpy as np
import subprocess as sp
import scipy.constants as pc
from collections import Counter
from distutils.version import LooseVersion
from numpy.linalg  import norm
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
from dpgen.generator.lib.abacus_pw_scf import make_abacus_pw_scf_stru, make_abacus_pw_scf_input, make_abacus_pw_scf_kpt
#from dpgen.generator.lib.pwscf import cvt_1frame
from dpgen.generator.lib.pwmat import make_pwmat_input_dict
from dpgen.generator.lib.pwmat import write_input_dict
from dpgen.generator.lib.pwmat import make_pwmat_input_user_dict
from dpgen.generator.lib.pwmat import input_upper
from dpgen.generator.lib.siesta import make_siesta_input
from dpgen.generator.lib.gaussian import make_gaussian_input, take_cluster
from dpgen.generator.lib.cp2k import make_cp2k_input, make_cp2k_input_from_external, make_cp2k_xyz
from dpgen.generator.lib.ele_temp import NBandsEsti
from dpgen.remote.RemoteJob import SSHSession, JobStatus, SlurmJob, PBSJob, LSFJob, CloudMachineJob, awsMachineJob
from dpgen.remote.group_jobs import ucloud_submit_jobs, aws_submit_jobs
from dpgen.remote.group_jobs import group_slurm_jobs
from dpgen.remote.group_jobs import group_local_jobs
from dpgen.remote.decide_machine import decide_train_machine, decide_fp_machine, decide_model_devi_machine
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, make_dispatcher, make_submission
from dpgen.util import sepline
from dpgen import ROOT_PATH
from pymatgen.io.vasp import Incar,Kpoints,Potcar
from dpgen.auto_test.lib.vasp import make_kspacing_kpoints
try:
    from gromacs.fileformats.mdp import MDP
except ImportError:
    dlog.info("GromacsWrapper>=0.8.0 is needed for DP-GEN + Gromacs.")
    pass

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
    for ii in range(numb_model):
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


def dump_to_poscar(dump, poscar, type_map, fmt = "lammps/dump") :
    sys = dpdata.System(dump, fmt = fmt, type_map = type_map)
    sys.to_vasp_poscar(poscar)

def dump_to_deepmd_raw(dump, deepmd_raw, type_map, fmt='gromacs/gro', charge=None):
    system = dpdata.System(dump, fmt = fmt, type_map = type_map)
    system.to_deepmd_raw(deepmd_raw)
    if charge is not None:
        with open(os.path.join(deepmd_raw, "charge"), 'w+') as f:
            f.write(str(charge))


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
    training_iter0_model = jdata.get('training_iter0_model_path', [])
    training_init_model = jdata.get('training_init_model', False)
    training_reuse_iter = jdata.get('training_reuse_iter')
    training_reuse_old_ratio = jdata.get('training_reuse_old_ratio', None)

    if 'training_reuse_stop_batch' in jdata.keys():
        training_reuse_stop_batch = jdata['training_reuse_stop_batch']
    elif 'training_reuse_numb_steps' in jdata.keys():
        training_reuse_stop_batch = jdata['training_reuse_numb_steps']
    else:
        training_reuse_stop_batch = 40000
        
    training_reuse_start_lr = jdata.get('training_reuse_start_lr', 1e-4)
    training_reuse_start_pref_e = jdata.get('training_reuse_start_pref_e', 0.1)
    training_reuse_start_pref_f = jdata.get('training_reuse_start_pref_f', 100)
    model_devi_activation_func = jdata.get('model_devi_activation_func', None)

    if training_reuse_iter is not None and training_reuse_old_ratio is None:
        raise RuntimeError("training_reuse_old_ratio not found but is mandatory when using init-model (training_reuse_iter is detected in param).\n" \
        "It defines the ratio of the old-data picking probability to the all-data(old-data plus new-data) picking probability in training after training_reuse_iter.\n" \
        "Denoting the index of the current iter as N (N >= training_reuse_iter ), old-data refers to those existed before the N-1 iter, and new-data refers to that obtained by the N-1 iter.\n" \
        "A recommended strategy is making the new-to-old ratio close to 10 times of the default value, to reasonably increase the sensitivity of the model to the new-data.\n" \
        "By default, the picking probability of data from one system or one iter is proportional to the number of batches (the number of frames divided by batch_size) of that systems or iter.\n" \
        "Detailed discussion about init-model (in Chinese) please see https://mp.weixin.qq.com/s/qsKMZ0j270YhQKvwXUiFvQ")

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
        if len(init_data_sys_) > len(init_batch_size_):
            warnings.warn("The batch sizes are not enough. Assume auto for those not spefified.")
            init_batch_size.extend(["auto" for aa in range(len(init_data_sys_)-len(init_batch_size))])
    else:
        init_batch_size_ = ["auto" for aa in range(len(jdata['init_data_sys']))]
    if 'sys_batch_size' in jdata:
        sys_batch_size = jdata['sys_batch_size']
    else:
        sys_batch_size = ["auto" for aa in range(len(jdata['sys_configs']))]

    # make sure all init_data_sys has the batch size -- for the following `zip`
    assert (len(init_data_sys_) <= len(init_batch_size_))
    for ii, ss in zip(init_data_sys_, init_batch_size_) :
        if jdata.get('init_multi_systems', False):
            for single_sys in os.listdir(os.path.join(work_path, 'data.init', ii)):
                init_data_sys.append(os.path.join('..', 'data.init', ii, single_sys))
                init_batch_size.append(detect_batch_size(ss, os.path.join(work_path, 'data.init', ii, single_sys)))
        else:
            init_data_sys.append(os.path.join('..', 'data.init', ii))
            init_batch_size.append(detect_batch_size(ss, os.path.join(work_path, 'data.init', ii)))
    old_range = None
    if iter_index > 0 :
        for ii in range(iter_index) :
            if ii == iter_index - 1:
                old_range = len(init_data_sys)
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
                    nframes = dpdata.System(jj, 'deepmd/npy').get_nframes()
                    if nframes < fp_task_min :
                        log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                        continue
                    init_data_sys.append(os.path.join('..', 'data.iters', jj))
                    init_batch_size.append(detect_batch_size(sys_batch_size[sys_idx], jj))
    # establish tasks
    jinput = jdata['default_training_param']
    try:
        mdata["deepmd_version"]
    except KeyError:
        mdata = set_version(mdata)
    # setup data systems
    if LooseVersion(mdata["deepmd_version"]) >= LooseVersion('1') and LooseVersion(mdata["deepmd_version"]) < LooseVersion('2'):
        # 1.x
        jinput['training']['systems'] = init_data_sys
        jinput['training']['batch_size'] = init_batch_size
        jinput['model']['type_map'] = jdata['type_map']
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
            raise RuntimeError('invalid setting for use_ele_temp ' + str(use_ele_temp))
    elif LooseVersion(mdata["deepmd_version"]) >= LooseVersion('2') and LooseVersion(mdata["deepmd_version"]) < LooseVersion('3'):
        # 2.x
        jinput['training']['training_data'] = {}
        jinput['training']['training_data']['systems'] = init_data_sys
        jinput['training']['training_data']['batch_size'] = init_batch_size
        jinput['model']['type_map'] = jdata['type_map']
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
            raise RuntimeError('invalid setting for use_ele_temp ' + str(use_ele_temp))
    else:
        raise RuntimeError("DP-GEN currently only supports for DeePMD-kit 1.x version!" )
    # set training reuse model
    if training_reuse_iter is not None and iter_index >= training_reuse_iter:
        if LooseVersion('1') <= LooseVersion(mdata["deepmd_version"]) < LooseVersion('2'):
            jinput['training']['stop_batch'] = training_reuse_stop_batch
            jinput['training']['auto_prob_style'] \
                ="prob_sys_size; 0:%d:%f; %d:%d:%f" \
                %(old_range, training_reuse_old_ratio, old_range, len(init_data_sys), 1.-training_reuse_old_ratio)
        elif LooseVersion('2') <= LooseVersion(mdata["deepmd_version"]) < LooseVersion('3'):
            jinput['training']['numb_steps'] = training_reuse_stop_batch
            jinput['training']['training_data']['auto_prob'] \
                ="prob_sys_size; 0:%d:%f; %d:%d:%f" \
                %(old_range, training_reuse_old_ratio, old_range, len(init_data_sys), 1.-training_reuse_old_ratio)
        else:
            raise RuntimeError("Unsupported DeePMD-kit version: %s" % mdata["deepmd_version"])
        if jinput['loss'].get('start_pref_e') is not None:
            jinput['loss']['start_pref_e'] = training_reuse_start_pref_e
        if jinput['loss'].get('start_pref_f') is not None:
            jinput['loss']['start_pref_f'] = training_reuse_start_pref_f
        jinput['learning_rate']['start_lr'] = training_reuse_start_lr
        

    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(task_path)
        os.chdir(task_path)
        for jj in init_data_sys :
            if not os.path.isdir(jj) :
                raise RuntimeError ("data sys %s does not exists, cwd is %s" % (jj, os.getcwd()))
        os.chdir(cwd)
        # set random seed for each model
        if LooseVersion(mdata["deepmd_version"]) >= LooseVersion('1') and LooseVersion(mdata["deepmd_version"]) < LooseVersion('3'):
            # 1.x
            if jinput['model']['descriptor']['type'] == 'hybrid':
                for desc in jinput['model']['descriptor']['list']:
                    desc['seed'] = random.randrange(sys.maxsize) % (2**32)
            else:
                jinput['model']['descriptor']['seed'] = random.randrange(sys.maxsize) % (2**32)
            jinput['model']['fitting_net']['seed'] = random.randrange(sys.maxsize) % (2**32)
            jinput['training']['seed'] = random.randrange(sys.maxsize) % (2**32)
        else:
            raise RuntimeError("DP-GEN currently only supports for DeePMD-kit 1.x version!" )
        # set model activation function
        if model_devi_activation_func is not None:
            if LooseVersion(mdata["deepmd_version"]) < LooseVersion('1'):
                raise RuntimeError('model_devi_activation_func does not suppport deepmd version', mdata['deepmd_version'])
            assert(type(model_devi_activation_func) is list and len(model_devi_activation_func) == numb_models)
            if len(np.array(model_devi_activation_func).shape) == 2 :                                    # 2-dim list for emd/fitting net-resolved assignment of actF
                jinput['model']['descriptor']['activation_function'] = model_devi_activation_func[ii][0]
                jinput['model']['fitting_net']['activation_function'] = model_devi_activation_func[ii][1]
            if len(np.array(model_devi_activation_func).shape) == 1 :                                    # for backward compatibility, 1-dim list, not net-resolved
                jinput['model']['descriptor']['activation_function'] = model_devi_activation_func[ii]
                jinput['model']['fitting_net']['activation_function'] = model_devi_activation_func[ii]
        # dump the input.json
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
            _link_old_models(work_path, old_model_files, ii)
    else:
        if type(training_iter0_model) == str:
            training_iter0_model = [training_iter0_model]
        iter0_models = []
        for ii in training_iter0_model:
            model_is = glob.glob(ii)
            model_is.sort()
            iter0_models += [os.path.abspath(ii) for ii in model_is]
        if training_init_model:
            assert(numb_models == len(iter0_models)), "training_iter0_model should be provided, and the number of models should be equal to %d" % numb_models
        for ii in range(len(iter0_models)):
            old_model_files = glob.glob(os.path.join(iter0_models[ii], 'model.ckpt*'))
            _link_old_models(work_path, old_model_files, ii)


def _link_old_models(work_path, old_model_files, ii):
    """
    link the `ii`th old model given by `old_model_files` to
    the `ii`th training task in `work_path`
    """
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
               mdata) :
    # print("debug:run_train:mdata", mdata)
    # load json param
    numb_models = jdata['numb_models']
    # train_param = jdata['train_param']
    train_input_file = default_train_input_file
    training_reuse_iter = jdata.get('training_reuse_iter')
    training_init_model = jdata.get('training_init_model', False)
    if training_reuse_iter is not None and iter_index >= training_reuse_iter:
        training_init_model = True
    try:
        mdata["deepmd_version"]
    except KeyError:
        mdata = set_version(mdata)

    
    train_command = mdata.get('train_command', 'dp')
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
    if LooseVersion(mdata["deepmd_version"]) >= LooseVersion('1') and LooseVersion(mdata["deepmd_version"]) < LooseVersion('3'):
        
        # 1.x
        ## Commands are like `dp train` and `dp freeze`
        ## train_command should not be None
        assert(train_command)
        command =  '%s train %s' % (train_command, train_input_file)
        if training_init_model:
            command = "{ if [ ! -f model.ckpt.index ]; then %s --init-model old/model.ckpt; else %s --restart model.ckpt; fi }" % (command, command)
        else:
            command = "{ if [ ! -f model.ckpt.index ]; then %s; else %s --restart model.ckpt; fi }" % (command, command)
        command = "/bin/sh -c '%s'" % command
        commands.append(command)
        command = '%s freeze' % train_command
        commands.append(command)
    else:
        raise RuntimeError("DP-GEN currently only supports for DeePMD-kit 1.x version!" )

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
    if training_init_model:
        forward_files += [os.path.join('old', 'model.ckpt.meta'),
                          os.path.join('old', 'model.ckpt.index'),
                          os.path.join('old', 'model.ckpt.data-00000-of-00001')
        ]
    backward_files = ['frozen_model.pb', 'lcurve.out', 'train.log']
    backward_files+= ['model.ckpt.meta', 'model.ckpt.index', 'model.ckpt.data-00000-of-00001', 'checkpoint']
    init_data_sys_ = jdata['init_data_sys']
    init_data_sys = []
    for ii in init_data_sys_ :
        init_data_sys.append(os.path.join('data.init', ii))
    trans_comm_data = []
    cwd = os.getcwd()
    os.chdir(work_path)
    fp_data = glob.glob(os.path.join('data.iters', 'iter.*', '02.fp', 'data.*'))
    for ii in init_data_sys :
        if jdata.get('init_multi_systems', False):
            for single_sys in os.listdir(os.path.join(ii)):
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'type*.raw'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'nopbc'))
        else:
            trans_comm_data += glob.glob(os.path.join(ii, 'set.*'))
            trans_comm_data += glob.glob(os.path.join(ii, 'type*.raw'))
            trans_comm_data += glob.glob(os.path.join(ii, 'nopbc'))
    for ii in fp_data :
        if jdata.get('use_clusters', False):
            for single_sys in os.listdir(os.path.join(ii)):
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'type*.raw'))
                trans_comm_data += glob.glob(os.path.join(ii, single_sys, 'nopbc'))
        else:
            trans_comm_data += glob.glob(os.path.join(ii, 'set.*'))
            trans_comm_data += glob.glob(os.path.join(ii, 'type*.raw'))
            trans_comm_data += glob.glob(os.path.join(ii, 'nopbc'))
    os.chdir(cwd)

    try:
        train_group_size = mdata['train_group_size']
    except:
        train_group_size = 1

    api_version = mdata.get('api_version', '0.9')
    # print('debug:commands', commands)

    if LooseVersion(api_version) < LooseVersion('1.0'):
        warnings.warn(f"the dpdispatcher will be updated to new version."
            f"And the interface may be changed. Please check the documents for more details")
        dispatcher = make_dispatcher(mdata['train_machine'], mdata['train_resources'], work_path, run_tasks, train_group_size)
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

    elif LooseVersion(api_version) >= LooseVersion('1.0'):
        submission = make_submission(
            mdata['train_machine'],
            mdata['train_resources'],
            commands=commands,
            work_path=work_path,
            run_tasks=run_tasks,
            group_size=train_group_size,
            forward_common_files=trans_comm_data,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog = 'train.log',
            errlog = 'train.log')
        submission.run_submission()

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
    elif 'nvt' == ensemble or 'nve' == ensemble:
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

def expand_matrix_values(target_list, cur_idx = 0):
    nvar = len(target_list)
    if cur_idx == nvar :
        return [[]]
    else :
        res = []
        prev = expand_matrix_values(target_list, cur_idx+1)
        for ii in target_list[cur_idx]:
            tmp = copy.deepcopy(prev)
            for jj in tmp:
               jj.insert(0, ii)
               res.append(jj)
        return res

def parse_cur_job_revmat(cur_job, use_plm = False):
    templates = [cur_job['template']['lmp']]
    if use_plm :
        templates.append(cur_job['template']['plm'])
    revise_keys = []
    revise_values = []
    if 'rev_mat' not in cur_job.keys():
        cur_job['rev_mat'] = {}
    if 'lmp' not in cur_job['rev_mat'].keys():
        cur_job['rev_mat']['lmp'] = {}
    for ii in cur_job['rev_mat']['lmp'].keys():
        revise_keys.append(ii)
        revise_values.append(cur_job['rev_mat']['lmp'][ii])
    n_lmp_keys = len(revise_keys)
    if use_plm:
        if 'plm' not in cur_job['rev_mat'].keys():
            cur_job['rev_mat']['plm'] = {}
        for ii in cur_job['rev_mat']['plm'].keys():
            revise_keys.append(ii)
            revise_values.append(cur_job['rev_mat']['plm'][ii])
    revise_matrix = expand_matrix_values(revise_values)
    return revise_keys, revise_matrix, n_lmp_keys


def parse_cur_job_sys_revmat(cur_job, sys_idx, use_plm=False):
    templates = [cur_job['template']['lmp']]
    if use_plm:
        templates.append(cur_job['template']['plm'])
    sys_revise_keys = []
    sys_revise_values = []
    if 'sys_rev_mat' not in cur_job.keys():
        cur_job['sys_rev_mat'] = {}
    local_rev = cur_job['sys_rev_mat'].get(str(sys_idx), {})
    if 'lmp' not in local_rev.keys():
        local_rev['lmp'] = {}
    for ii in local_rev['lmp'].keys():
        sys_revise_keys.append(ii)
        sys_revise_values.append(local_rev['lmp'][ii])
    n_sys_lmp_keys = len(sys_revise_keys)
    if use_plm:
        if 'plm' not in local_rev.keys():
            local_rev['plm'] = {}
        for ii in local_rev['plm'].keys():
            sys_revise_keys.append(ii)
            sys_revise_values.append(local_rev['plm'][ii])
    sys_revise_matrix = expand_matrix_values(sys_revise_values)
    return sys_revise_keys, sys_revise_matrix, n_sys_lmp_keys

def find_only_one_key(lmp_lines, key):
    found = []
    for idx in range(len(lmp_lines)):
        words = lmp_lines[idx].split()
        nkey = len(key)
        if len(words) >= nkey and words[:nkey] == key :
            found.append(idx)
    if len(found) > 1:
        raise RuntimeError('found %d keywords %s' % (len(found), key))
    if len(found) == 0:
        raise RuntimeError('failed to find keyword %s' % (key))
    return found[0]


def revise_lmp_input_model(lmp_lines, task_model_list, trj_freq, deepmd_version = '1'):
    idx = find_only_one_key(lmp_lines, ['pair_style', 'deepmd'])
    graph_list = ' '.join(task_model_list)
    if LooseVersion(deepmd_version) < LooseVersion('1'):
        lmp_lines[idx] = "pair_style      deepmd %s %d model_devi.out\n" % (graph_list, trj_freq)
    else:
        lmp_lines[idx] = "pair_style      deepmd %s out_freq %d out_file model_devi.out\n" % (graph_list, trj_freq)
    return lmp_lines


def revise_lmp_input_dump(lmp_lines, trj_freq):
    idx = find_only_one_key(lmp_lines, ['dump', 'dpgen_dump'])
    lmp_lines[idx] = "dump            dpgen_dump all custom %d traj/*.lammpstrj id type x y z\n" % trj_freq
    return lmp_lines


def revise_lmp_input_plm(lmp_lines, in_plm, out_plm = 'output.plumed'):
    idx = find_only_one_key(lmp_lines, ['fix', 'dpgen_plm'])
    lmp_lines[idx] = "fix            dpgen_plm all plumed plumedfile %s outfile %s\n" % (in_plm, out_plm)
    return lmp_lines


def revise_by_keys(lmp_lines, keys, values):
    for kk,vv in zip(keys, values):
        for ii in range(len(lmp_lines)):
            lmp_lines[ii] = lmp_lines[ii].replace(kk, str(vv))
    return lmp_lines


def make_model_devi (iter_index,
                     jdata,
                     mdata) :
    # The MD engine to perform model deviation
    # Default is lammps
    model_devi_engine = jdata.get('model_devi_engine', "lammps")
    model_devi_jobs = jdata['model_devi_jobs']
    if (iter_index >= len(model_devi_jobs)) :
        return False
    cur_job = model_devi_jobs[iter_index]
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

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = glob.glob(os.path.join(train_path, "graph*pb"))
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
            if model_devi_engine == "lammps":
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
                if jdata.get('model_devi_nopbc', False):
                    system.remove_pbc()
                system.to_lammps_lmp(os.path.join(conf_path, lmp_name))
            elif model_devi_engine == "gromacs":
                pass
            conf_counter += 1
        sys_counter += 1

    input_mode = "native"
    if "template" in cur_job:
        input_mode = "revise_template"
    use_plm = jdata.get('model_devi_plumed', False)
    use_plm_path = jdata.get('model_devi_plumed_path', False)
    if input_mode == "native":
        if model_devi_engine == "lammps":
            _make_model_devi_native(iter_index, jdata, mdata, conf_systems)
        elif model_devi_engine == "gromacs":
            _make_model_devi_native_gromacs(iter_index, jdata, mdata, conf_systems)
        else:
            raise RuntimeError("unknown model_devi engine", model_devi_engine)
    elif input_mode == "revise_template":
        _make_model_devi_revmat(iter_index, jdata, mdata, conf_systems)
    else:
        raise RuntimeError('unknown model_devi input mode', input_mode)

    return True


def _make_model_devi_revmat(iter_index, jdata, mdata, conf_systems):
    model_devi_jobs = jdata['model_devi_jobs']
    if (iter_index >= len(model_devi_jobs)) :
        return False
    cur_job = model_devi_jobs[iter_index]
    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")
    mass_map = jdata['mass_map']
    use_plm = jdata.get('model_devi_plumed', False)
    use_plm_path = jdata.get('model_devi_plumed_path', False)
    trj_freq = _get_param_alias(cur_job, ['t_freq', 'trj_freq','traj_freq'])

    rev_keys, rev_mat, num_lmp = parse_cur_job_revmat(cur_job, use_plm = use_plm)
    lmp_templ = cur_job['template']['lmp']
    lmp_templ = os.path.abspath(lmp_templ)
    if use_plm:
        plm_templ = cur_job['template']['plm']
        plm_templ = os.path.abspath(plm_templ)
        if use_plm_path:
            plm_path_templ = cur_job['template']['plm_path']
            plm_path_templ = os.path.abspath(plm_path_templ)

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = glob.glob(os.path.join(train_path, "graph*pb"))
    task_model_list = []
    for ii in models:
        task_model_list.append(os.path.join('..', os.path.basename(ii)))
    work_path = os.path.join(iter_name, model_devi_name)
    try:
        mdata["deepmd_version"]
    except KeyError:
        mdata = set_version(mdata)
    deepmd_version = mdata['deepmd_version']

    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        task_counter = 0
        for cc in ss :
            sys_rev = cur_job.get('sys_rev_mat', None)
            total_rev_keys = rev_keys
            total_rev_mat = rev_mat
            total_num_lmp = num_lmp
            if sys_rev is not None:
                total_rev_mat = []
                sys_rev_keys, sys_rev_mat, sys_num_lmp = parse_cur_job_sys_revmat(cur_job,
                                                                                  sys_idx=sys_idx[sys_counter],
                                                                                  use_plm=use_plm)
                _lmp_keys = rev_keys[:num_lmp] + sys_rev_keys[:sys_num_lmp]
                if use_plm:
                    _plm_keys = rev_keys[num_lmp:] + sys_rev_keys[sys_num_lmp:]
                    _lmp_keys += _plm_keys
                total_rev_keys = _lmp_keys
                total_num_lmp = num_lmp + sys_num_lmp
                for pub in rev_mat:
                    for pri in sys_rev_mat:
                        _lmp_mat = pub[:num_lmp] + pri[:sys_num_lmp]
                        if use_plm:
                            _plm_mat = pub[num_lmp:] + pri[sys_num_lmp:]
                            _lmp_mat += _plm_mat
                        total_rev_mat.append(_lmp_mat)
            for ii in range(len(total_rev_mat)):
                total_rev_item = total_rev_mat[ii]
                task_name = make_model_devi_task_name(sys_idx[sys_counter], task_counter)
                conf_name = make_model_devi_conf_name(sys_idx[sys_counter], conf_counter) + '.lmp'
                task_path = os.path.join(work_path, task_name)
                # create task path
                create_path(task_path)
                create_path(os.path.join(task_path, 'traj'))
                # link conf
                loc_conf_name = 'conf.lmp'
                os.symlink(os.path.join(os.path.join('..','confs'), conf_name),
                           os.path.join(task_path, loc_conf_name) )
                cwd_ = os.getcwd()
                # chdir to task path
                os.chdir(task_path)
                shutil.copyfile(lmp_templ, 'input.lammps')
                # revise input of lammps
                with open('input.lammps') as fp:
                    lmp_lines = fp.readlines()
                lmp_lines = revise_lmp_input_model(lmp_lines, task_model_list, trj_freq, deepmd_version = deepmd_version)
                lmp_lines = revise_lmp_input_dump(lmp_lines, trj_freq)
                lmp_lines = revise_by_keys(
                    lmp_lines, total_rev_keys[:total_num_lmp], total_rev_item[:total_num_lmp]
                )
                # revise input of plumed
                if use_plm:
                    lmp_lines = revise_lmp_input_plm(lmp_lines, 'input.plumed')
                    shutil.copyfile(plm_templ, 'input.plumed')
                    with open('input.plumed') as fp:
                        plm_lines = fp.readlines()
                    # allow using the same list as lmp
                    # user should not use the same key name for plm
                    plm_lines = revise_by_keys(
                        plm_lines, total_rev_keys, total_rev_item
                    )
                    with open('input.plumed', 'w') as fp:
                        fp.write(''.join(plm_lines))
                    if use_plm_path:
                       shutil.copyfile(plm_path_templ, 'plmpath.pdb')
                # dump input of lammps
                with open('input.lammps', 'w') as fp:
                    fp.write(''.join(lmp_lines))
                with open('job.json', 'w') as fp:
                    job = {}
                    for ii,jj in zip(total_rev_keys, total_rev_item) : job[ii] = jj
                    json.dump(job, fp, indent = 4)
                os.chdir(cwd_)
                task_counter += 1
            conf_counter += 1
        sys_counter += 1


def _make_model_devi_native(iter_index, jdata, mdata, conf_systems):
    model_devi_jobs = jdata['model_devi_jobs']
    if (iter_index >= len(model_devi_jobs)) :
        return False
    cur_job = model_devi_jobs[iter_index]
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt = parse_cur_job(cur_job)
    if dt is not None :
        model_devi_dt = dt
    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")

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
    mass_map = jdata['mass_map']
    nopbc = jdata.get('model_devi_nopbc', False)

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = glob.glob(os.path.join(train_path, "graph*pb"))
    task_model_list = []
    for ii in models:
        task_model_list.append(os.path.join('..', os.path.basename(ii)))
    work_path = os.path.join(iter_name, model_devi_name)

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
                    except KeyError:
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
                                               nopbc = nopbc,
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

def _make_model_devi_native_gromacs(iter_index, jdata, mdata, conf_systems):
    # only support for deepmd v2.0
    if LooseVersion(mdata['deepmd_version']) < LooseVersion('2.0'):
        raise RuntimeError("Only support deepmd-kit 2.x for model_devi_engine='gromacs'")
    model_devi_jobs = jdata['model_devi_jobs']
    if (iter_index >= len(model_devi_jobs)) :
        return False
    cur_job = model_devi_jobs[iter_index]
    dt = cur_job.get("dt", None)
    if dt is not None:
        model_devi_dt = dt
    else:
        model_devi_dt = jdata['model_devi_dt']
    nsteps = cur_job.get("nsteps", None)
    lambdas = cur_job.get("lambdas", [])
    temps = cur_job.get("temps", [])
    if not lambdas:
        lambdas = [1.0]
    else:
        for ll in lambdas:
            if ll > 1:
                raise RuntimeError("lambda is larger than 1.0")
    if not temps:
        temps = [298.0]
    if nsteps is None:
        raise RuntimeError("nsteps is None, you should set nsteps in model_devi_jobs!")
    # Currently Gromacs engine is not supported for different temperatures!
    # If you want to change temperatures, you should change it in mdp files.
  
    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")

    mass_map = jdata['mass_map']

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = glob.glob(os.path.join(train_path, "graph*pb"))
    task_model_list = []
    for ii in models:
        task_model_list.append(os.path.join('..', os.path.basename(ii)))
    work_path = os.path.join(iter_name, model_devi_name)

    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        task_counter = 0
        for cc in ss :
            for ll in lambdas:
                for tt in temps:
                    task_name = make_model_devi_task_name(sys_idx[sys_counter], task_counter)
                    task_path = os.path.join(work_path, task_name)
                    create_path(task_path)
                    gromacs_settings = jdata.get("gromacs_settings" , "")
                    for key,file in gromacs_settings.items():
                        if key != "traj_filename" and key != "mdp_filename" and key != "group_name" and key != "maxwarn":
                            os.symlink(os.path.join(cc,file), os.path.join(task_path, file))
                    # input.json for DP-Gromacs
                    with open(os.path.join(cc, "input.json")) as f:
                        input_json = json.load(f)
                    input_json["graph_file"] = models[0]
                    input_json["lambda"] = ll
                    with open(os.path.join(task_path,'input.json'), 'w') as _outfile:
                        json.dump(input_json, _outfile, indent = 4)

                    # trj_freq
                    trj_freq = cur_job.get("trj_freq", 10)
                    mdp = MDP()
                    mdp.read(os.path.join(cc, gromacs_settings['mdp_filename']))
                    mdp['nstcomm'] = trj_freq
                    mdp['nstxout'] = trj_freq
                    mdp['nstlog'] = trj_freq
                    mdp['nstenergy'] = trj_freq
                    # dt
                    mdp['dt'] = model_devi_dt
                    # nsteps
                    mdp['nsteps'] = nsteps
                    # temps
                    if "ref_t" in list(mdp.keys()):
                        mdp["ref_t"] = tt
                    else:
                        mdp["ref-t"] = tt
                    mdp.write(os.path.join(task_path, gromacs_settings['mdp_filename']))

                    cwd_ = os.getcwd()
                    os.chdir(task_path)
                    job = {}
                    job["trj_freq"] = cur_job["trj_freq"]
                    job["model_devi_dt"] =  model_devi_dt
                    job["nsteps"] = nsteps
                    with open('job.json', 'w') as _outfile:
                        json.dump(job, _outfile, indent = 4)
                    os.chdir(cwd_) 
                    task_counter += 1
            conf_counter += 1
        sys_counter += 1



def run_model_devi (iter_index,
                    jdata,
                    mdata) :
    #rmdlog.info("This module has been run !")
    lmp_exec = mdata['lmp_command']
    # Anguse: lmp_exec name should be changed to model_devi_exec.
    # We should also change make_dispatcher
    # For now, I will use this name for gromacs command

    model_devi_group_size = mdata['model_devi_group_size']
    model_devi_resources = mdata['model_devi_resources']
    use_plm = jdata.get('model_devi_plumed', False)
    use_plm_path = jdata.get('model_devi_plumed_path', False)

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))

    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    fp = open (os.path.join(work_path, 'cur_job.json'), 'r')
    cur_job = json.load (fp)

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

    model_devi_engine = jdata.get("model_devi_engine", "lammps")
    if model_devi_engine == "lammps":
        command = "{ if [ ! -f dpgen.restart.10000 ]; then %s -i input.lammps -v restart 0; else %s -i input.lammps -v restart 1; fi }" % (lmp_exec, lmp_exec)
        command = "/bin/sh -c '%s'" % command
        commands = [command]
        forward_files = ['conf.lmp', 'input.lammps', 'traj']
        backward_files = ['model_devi.out', 'model_devi.log', 'traj']
        if use_plm:
            forward_files += ['input.plumed']
           # backward_files += ['output.plumed']
            backward_files += ['output.plumed','COLVAR','dump.0.xyz']
            if use_plm_path:
                forward_files += ['plmpath.pdb']
    elif model_devi_engine == "gromacs":
        
        gromacs_settings = jdata.get("gromacs_settings", {})
        mdp_filename = gromacs_settings.get("mdp_filename", "md.mdp")
        topol_filename = gromacs_settings.get("topol_filename", "processed.top")
        conf_filename = gromacs_settings.get("conf_filename", "conf.gro")
        index_filename = gromacs_settings.get("index_filename", "index.raw")
        type_filename = gromacs_settings.get("type_filename", "type.raw")
        ndx_filename = gromacs_settings.get("ndx_filename", "")
        # Initial reference to process pbc condition.
        # Default is em.tpr
        ref_filename = gromacs_settings.get("ref_filename", "em.tpr")
        deffnm = gromacs_settings.get("deffnm", "deepmd")
        maxwarn = gromacs_settings.get("maxwarn", 1)
        traj_filename = gromacs_settings.get("traj_filename", "deepmd_traj.gro")
        grp_name = gromacs_settings.get("group_name", "Other")
        trj_freq = cur_job.get("trj_freq", 10)

        command = "%s grompp -f %s -p %s -c %s -o %s -maxwarn %d" % (lmp_exec, mdp_filename, topol_filename, conf_filename, deffnm, maxwarn)
        command += "&& %s mdrun -deffnm %s -cpi" %(lmp_exec, deffnm)
        if ndx_filename:
            command += f"&& echo -e \"{grp_name}\\n{grp_name}\\n\" | {lmp_exec} trjconv -s {ref_filename} -f {deffnm}.trr -n {ndx_filename} -o {traj_filename} -pbc mol -ur compact -center"
        else:
            command += "&& echo -e \"%s\\n%s\\n\" | %s trjconv -s %s -f %s.trr -o %s -pbc mol -ur compact -center" % (grp_name, grp_name, lmp_exec, ref_filename, deffnm, traj_filename)
        command += "&& if [ ! -d traj ]; then \n mkdir traj; fi\n"
        command += f"python -c \"import dpdata;system = dpdata.System('{traj_filename}', fmt='gromacs/gro'); [system.to_gromacs_gro('traj/%d.gromacstrj' % (i * {trj_freq}), frame_idx=i) for i in range(system.get_nframes())]; system.to_deepmd_npy('traj_deepmd')\""
        command += f"&& dp model-devi -m ../graph.000.pb ../graph.001.pb ../graph.002.pb ../graph.003.pb -s traj_deepmd -o model_devi.out -f {trj_freq}"
        commands = [command]

        forward_files = [mdp_filename, topol_filename, conf_filename, index_filename, ref_filename, type_filename, "input.json", "job.json" ]
        if ndx_filename: forward_files.append(ndx_filename)
        backward_files = ["%s.tpr" % deffnm, "%s.log" %deffnm , traj_filename, 'model_devi.out', "traj", "traj_deepmd" ]


    cwd = os.getcwd()

    api_version = mdata.get('api_version', '0.9')
    if LooseVersion(api_version) < LooseVersion('1.0'):
        warnings.warn(f"the dpdispatcher will be updated to new version."
            f"And the interface may be changed. Please check the documents for more details")
        dispatcher = make_dispatcher(mdata['model_devi_machine'], mdata['model_devi_resources'], work_path, run_tasks, model_devi_group_size)
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

    elif LooseVersion(api_version) >= LooseVersion('1.0'):
        submission = make_submission(
            mdata['model_devi_machine'],
            mdata['model_devi_resources'],
            commands=commands,
            work_path=work_path,
            run_tasks=run_tasks,
            group_size=model_devi_group_size,
            forward_common_files=model_names,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog = 'model_devi.log',
            errlog = 'model_devi.log')
        submission.run_submission()

def post_model_devi (iter_index,
                     jdata,
                     mdata) :
    pass


def _to_face_dist(box_):
    box = np.reshape(box_, [3,3])
    vol = np.abs(np.linalg.det(box))
    dists = []
    for [ii,jj] in [[0, 1], [1, 2], [2, 0]]:
        vv = np.cross(box[ii], box[jj])
        dists.append(vol / np.linalg.norm(vv))
    return np.array(dists)

def check_cluster(conf_name,
                  fp_cluster_vacuum,
                  fmt='lammps/dump'):
    sys = dpdata.System(conf_name, fmt)
    assert(sys.get_nframes() == 1)
    cell=sys.data['cells'][0]
    coord=sys.data['coords'][0]
    xlim=max(coord[:,0])-min(coord[:,0])
    ylim=max(coord[:,1])-min(coord[:,1])
    zlim=max(coord[:,2])-min(coord[:,2])
    a,b,c=map(norm,[cell[0,:],cell[1,:],cell[2,:]])
    min_vac=min([a-xlim,b-ylim,c-zlim])
    #print([a-xlim,b-ylim,c-zlim])
    #_,r3d=miniball.get_bounding_ball(coord)

    if min_vac < fp_cluster_vacuum:
       is_bad = True
    else:
       is_bad = False
    return is_bad

def check_bad_box(conf_name,
                   criteria,
                   fmt = 'lammps/dump'):
    all_c = criteria.split(';')
    sys = dpdata.System(conf_name, fmt)
    assert(sys.get_nframes() == 1)
    is_bad = False
    for ii in all_c:
        [key, value] = ii.split(':')
        if key == 'length_ratio':
            lengths = np.linalg.norm(sys['cells'][0], axis = 1)
            ratio = np.max(lengths) / np.min(lengths)
            if ratio > float(value):
                is_bad = True
        elif key == 'height_ratio':
            lengths = np.linalg.norm(sys['cells'][0], axis = 1)
            dists = _to_face_dist(sys['cells'][0])
            ratio = np.max(lengths) / np.min(dists)
            if ratio > float(value):
                is_bad = True
        else:
            raise RuntimeError('unknow key', key)
    return is_bad

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

    charges_recorder = []  # record charges for each fp_task
    charges_map = jdata.get("sys_charges", [])

    cluster_cutoff = jdata['cluster_cutoff'] if jdata.get('use_clusters', False) else None
    # skip save *.out if detailed_report_make_fp is False, default is True
    detailed_report_make_fp = jdata.get("detailed_report_make_fp", True)
    # skip bad box criteria
    skip_bad_box = jdata.get('fp_skip_bad_box')
    # skip discrete structure in cluster
    fp_cluster_vacuum = jdata.get('fp_cluster_vacuum',None)
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
        counter['candidate'] = 0
        counter['failed'] = 0
        counter['accurate'] = 0
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

        # set number of tasks
        accurate_ratio = float(counter['accurate']) / float(fp_sum)
        fp_accurate_threshold = jdata.get('fp_accurate_threshold', 1)
        fp_accurate_soft_threshold = jdata.get('fp_accurate_soft_threshold', fp_accurate_threshold)
        if accurate_ratio < fp_accurate_soft_threshold :
            this_fp_task_max = fp_task_max
        elif accurate_ratio >= fp_accurate_soft_threshold and accurate_ratio < fp_accurate_threshold:
            this_fp_task_max = int(fp_task_max * (accurate_ratio - fp_accurate_threshold) / (fp_accurate_soft_threshold - fp_accurate_threshold))
        else:
            this_fp_task_max = 0
        numb_task = min(this_fp_task_max, len(fp_candidate))
        if (numb_task < fp_task_min):
            numb_task = 0
        dlog.info("system {0:s} accurate_ratio: {1:8.4f}    thresholds: {2:6.4f} and {3:6.4f}   eff. task min and max {4:4d} {5:4d}   number of fp tasks: {6:6d}".format(ss, accurate_ratio, fp_accurate_soft_threshold, fp_accurate_threshold, fp_task_min, this_fp_task_max, numb_task))
        # make fp tasks
        model_devi_engine = jdata.get("model_devi_engine", "lammps")
        count_bad_box = 0
        count_bad_cluster = 0
        for cc in range(numb_task) :
            tt = fp_candidate[cc][0]
            ii = fp_candidate[cc][1]
            ss = os.path.basename(tt).split('.')[1]
            conf_name = os.path.join(tt, "traj")
            if model_devi_engine == "lammps":
                conf_name = os.path.join(conf_name, str(ii) + '.lammpstrj')
            elif model_devi_engine == "gromacs":
                conf_name = os.path.join(conf_name, str(ii) + '.gromacstrj')
            else:
                raise RuntimeError("unknown model_devi engine", model_devi_engine)
            conf_name = os.path.abspath(conf_name)
            if skip_bad_box is not None:
                skip = check_bad_box(conf_name, skip_bad_box)
                if skip:
                    count_bad_box += 1
                    continue

            if fp_cluster_vacuum is not None:
                assert fp_cluster_vacuum >0
                skip_cluster = check_cluster(conf_name, fp_cluster_vacuum)
                if skip_cluster:
                    count_bad_cluster +=1
                    continue

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
            if charges_map:
                charges_recorder.append(charges_map[int(ss)])
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
        if count_bad_box > 0:
            dlog.info("system {0:s} skipped {1:6d} confs with bad box, {2:6d} remains".format(ss, count_bad_box, numb_task - count_bad_box))
        if count_bad_cluster > 0:
            dlog.info("system {0:s} skipped {1:6d} confs with bad cluster, {2:6d} remains".format(ss, count_bad_cluster, numb_task - count_bad_cluster))
    if cluster_cutoff is None:
        cwd = os.getcwd()
        for idx, task in enumerate(fp_tasks):
            os.chdir(task)
            if model_devi_engine == "lammps":
                dump_to_poscar('conf.dump', 'POSCAR', type_map, fmt = "lammps/dump")
                if charges_map:
                    warnings.warn('"sys_charges" keyword only support for gromacs engine now.')
            elif model_devi_engine == "gromacs":
                # dump_to_poscar('conf.dump', 'POSCAR', type_map, fmt = "gromacs/gro")
                if charges_map:
                    dump_to_deepmd_raw('conf.dump', 'deepmd.raw', type_map, fmt='gromacs/gro', charge=charges_recorder[idx])
                else:
                    dump_to_deepmd_raw('conf.dump', 'deepmd.raw', type_map, fmt='gromacs/gro', charge=None)
            else:
                raise RuntimeError("unknown model_devi engine", model_devi_engine)
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

def make_pwmat_input(jdata, filename):
    if 'fp_incar' in jdata.keys() :
        fp_incar_path = jdata['fp_incar']
        assert(os.path.exists(fp_incar_path))
        fp_incar_path = os.path.abspath(fp_incar_path)
        fr = open(fp_incar_path)
        input = fr.read()
        fr.close()
    elif 'user_fp_params' in jdata.keys() :
        fp_params = jdata['user_fp_params']
        node1 = fp_params['node1']
        node2 = fp_params['node2']
        atom_config = fp_params['in.atom']
        ecut = fp_params['ecut']
        e_error = fp_params['e_error']
        rho_error = fp_params['rho_error']
        kspacing = fp_params['kspacing']
        flag_symm = fp_params['flag_symm']
        os.system("command -v poscar2config.x | wc -l > 1.txt")
        fc = open('1.txt')
        flag_command = fc.read()
        fc.close()
        if int(flag_command) == 1 :
            os.system('poscar2config.x < POSCAR > tmp.config')
        else:
            os.system('cp ../../../out_data_post_fp_pwmat/02.fp/task.000.000000/poscar2config.x ./')
            os.system('./poscar2config.x < POSCAR > tmp.config')
        os.system('rm -rf tmp.config')
        input_dict = make_pwmat_input_dict(node1, node2, atom_config, ecut, e_error,
                                           rho_error, icmix = None, smearing = None,
                                           sigma = None, kspacing = kspacing,
                                           flag_symm = flag_symm
        )

        input = write_input_dict(input_dict)
    else:
        input = make_pwmat_input_user_dict(jdata['fp_params'])
    if 'IN.PSP' in input or 'in.psp' in input:
        with open(filename, 'w') as fp:
            fp.write(input)
            fp.write('job=scf\n')
            if 'OUT.MLMD' in input or 'out.mlmd' in input:
                return input
            else:
                fp.write('OUT.MLMD = T')
                return input
    else:
        with open(filename, 'w') as fp:
            fp.write(input)
            fp.write('job=scf\n')
            fp_pp_files = jdata['fp_pp_files']
            for idx, ii in enumerate(fp_pp_files) :
                fp.write('IN.PSP%d = %s\n' %(idx+1, ii))
            if 'OUT.MLMD' in input or 'out.mlmd' in input:
                return input
            else:
                fp.write('OUT.MLMD = T')
                return input

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

def make_fp_vasp_incar (iter_index,
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

def _make_fp_pwmat_input (iter_index,
                         jdata) :
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        make_pwmat_input(jdata, 'etot.input')
        os.system("sed -i '1,2c 4 1' etot.input")
        os.chdir(cwd)

def make_fp_vasp_cp_cvasp(iter_index,jdata):
    # Move cvasp interface to jdata
    if ('cvasp' in jdata) and (jdata['cvasp'] == True):
       pass
    else:
       return
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        #copy cvasp.py
        shutil.copyfile(cvasp_file, 'cvasp.py')
        os.chdir(cwd)

def make_fp_vasp_kp (iter_index,jdata):
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_aniso_kspacing = jdata.get('fp_aniso_kspacing')

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
        if fp_aniso_kspacing is None:
            try:
                kspacing = standard_incar['KSPACING']
            except KeyError:
                raise RuntimeError ("KSPACING must be given in INCAR")
        else :
            kspacing = fp_aniso_kspacing
        try:
            gamma = standard_incar['KGAMMA']
            if isinstance(gamma,bool):
                pass
            else:
                if gamma[0].upper()=="T":
                    gamma=True
                else:
                    gamma=False
        except KeyError:
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
    make_fp_vasp_incar(iter_index, jdata, nbands_esti = nbe)
    # 3, create kpoints
    make_fp_vasp_kp(iter_index, jdata)
    # 4, copy cvasp
    make_fp_vasp_cp_cvasp(iter_index,jdata)


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

def make_fp_abacus_pw_scf(iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return
    # make abacus/pw/scf input
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_pp_files = jdata['fp_pp_files']
    if 'user_fp_params' in jdata.keys() :
        fp_params = jdata['user_fp_params']
        #user_input = True
    else:
        raise RuntimeError("Key 'user_fp_params' and its value have to be specified in parameter json file.")
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        sys_data = dpdata.System('POSCAR').data
        if 'mass_map' in jdata:
            sys_data['atom_masses'] = jdata['mass_map']
        ret_input = make_abacus_pw_scf_input(fp_params)
        with open('INPUT', 'w') as fp:
            fp.write(ret_input)
        ret_kpt = make_abacus_pw_scf_kpt(fp_params)
        with open("KPT", "w") as fp:
            fp.write(ret_kpt)
        ret_stru = make_abacus_pw_scf_stru(sys_data, fp_pp_files)
        with open("STRU", "w") as fp:
            fp.write(ret_stru)

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

    model_devi_engine = jdata.get('model_devi_engine', 'lammps')
    for ii in fp_tasks:
        os.chdir(ii)
        if model_devi_engine == "lammps":
            sys_data = dpdata.System('POSCAR').data
        elif model_devi_engine == "gromacs":
            sys_data = dpdata.System("deepmd.raw", fmt='deepmd/raw').data
            if os.path.isfile('deepmd.raw/charge'):
                sys_data['charge'] = int(np.loadtxt('deepmd.raw/charge', dtype=int))
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
    # some users might use own inputs
    # specify the input path string
    elif 'external_input_path' in jdata.keys() :
        fp_params = None
        exinput_path = os.path.abspath(jdata['external_input_path'])
    else:
        fp_params = jdata['fp_params']
    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        sys_data = dpdata.System('POSCAR').data
        # make input for every task
        # if fp_params exits, make keys
        if fp_params:
            cp2k_input = make_cp2k_input(sys_data, fp_params)
        else:
        # else read from user input
            cp2k_input = make_cp2k_input_from_external(sys_data, exinput_path)
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

def make_fp_pwmat (iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return
    # abs path for fp_incar if it exists
    if 'fp_incar' in jdata:
        jdata['fp_incar'] = os.path.abspath(jdata['fp_incar'])
    # order is critical!
    # 1, link pp files
    _link_fp_vasp_pp(iter_index, jdata)
    # 2, create pwmat input
    _make_fp_pwmat_input(iter_index, jdata)

def make_fp (iter_index,
             jdata,
             mdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        make_fp_vasp(iter_index, jdata)
    elif fp_style == "pwscf" :
        make_fp_pwscf(iter_index, jdata)
    elif fp_style == "abacus/scf" :
        make_fp_abacus_pw_scf(iter_index, jdata)
    elif fp_style == "siesta" :
        make_fp_siesta(iter_index, jdata)
    elif fp_style == "gaussian" :
        make_fp_gaussian(iter_index, jdata)
    elif fp_style == "cp2k" :
        make_fp_cp2k(iter_index, jdata)
    elif fp_style == "pwmat" :
        make_fp_pwmat(iter_index, jdata)
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

def _abacus_pw_scf_check_fin(ii) :
    if os.path.isfile(os.path.join(ii, 'OUT.ABACUS/running_scf.log')) :
        with open(os.path.join(ii, 'OUT.ABACUS/running_scf.log'), 'r') as fp :
            content = fp.read()
            count = content.count('!FINAL_ETOT_IS')
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

def _pwmat_check_fin (ii) :
    if os.path.isfile(os.path.join(ii, 'REPORT')) :
        with open(os.path.join(ii, 'REPORT'), 'r') as fp :
            content = fp.read()
            count = content.count('time')
            if count != 1 :
                return False
    else :
        return False
    return True

def run_fp_inner (iter_index,
                  jdata,
                  mdata,
                  forward_files,
                  backward_files,
                  check_fin,
                  log_file = "fp.log",
                  forward_common_files=[]) :
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']
    mark_failure = fp_resources.get('mark_failure', False)

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

    api_version = mdata.get('api_version', '0.9')
    if LooseVersion(api_version) < LooseVersion('1.0'):
        warnings.warn(f"the dpdispatcher will be updated to new version."
            f"And the interface may be changed. Please check the documents for more details")
        dispatcher = make_dispatcher(mdata['fp_machine'], mdata['fp_resources'], work_path, run_tasks, fp_group_size)
        dispatcher.run_jobs(mdata['fp_resources'],
                        [fp_command],
                        work_path,
                        run_tasks,
                        fp_group_size,
                        forward_common_files,
                        forward_files,
                        backward_files,
                        mark_failure = mark_failure,
                        outlog = log_file,
                        errlog = log_file)

    elif LooseVersion(api_version) >= LooseVersion('1.0'):
        submission = make_submission(
            mdata['fp_machine'],
            mdata['fp_resources'],
            commands=[fp_command],
            work_path=work_path,
            run_tasks=run_tasks,
            group_size=fp_group_size,
            forward_common_files=forward_common_files,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog = log_file,
            errlog = log_file)
        submission.run_submission()


def run_fp (iter_index,
            jdata,
            mdata) :
    fp_style = jdata['fp_style']
    fp_pp_files = jdata['fp_pp_files']

    if fp_style == "vasp" :
        forward_files = ['POSCAR', 'INCAR', 'POTCAR','KPOINTS']
        backward_files = ['OUTCAR','vasprun.xml']
        # Move cvasp interface to jdata
        if ('cvasp' in jdata) and (jdata['cvasp'] == True):
            mdata['fp_resources']['cvasp'] = True
        if ('cvasp' in  mdata["fp_resources"] ) and (mdata["fp_resources"]["cvasp"]==True):
            dlog.info("cvasp is on !")
            forward_files.append('cvasp.py')
            forward_common_files=[]
        else:
            forward_common_files=[]
        run_fp_inner(iter_index, jdata, mdata,  forward_files, backward_files, _vasp_check_fin,
                     forward_common_files=forward_common_files)
    elif fp_style == "pwscf" :
        forward_files = ['input'] + fp_pp_files
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata,  forward_files, backward_files, _qe_check_fin, log_file = 'output')
    elif fp_style == "abacus/scf":
        forward_files = ["INPUT", "STRU", "KPT"] + fp_pp_files
        backward_files = ["output", "OUT.ABACUS"]
        run_fp_inner(iter_index, jdata, mdata,  forward_files, backward_files, _abacus_pw_scf_check_fin, log_file = 'output')
    elif fp_style == "siesta":
        forward_files = ['input'] + fp_pp_files
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata,  forward_files, backward_files, _siesta_check_fin, log_file='output')
    elif fp_style == "gaussian":
        forward_files = ['input']
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, forward_files, backward_files, _gaussian_check_fin, log_file = 'output')
    elif fp_style == "cp2k":
        forward_files = ['input.inp', 'coord.xyz']
        backward_files = ['output']
        run_fp_inner(iter_index, jdata, mdata, forward_files, backward_files, _cp2k_check_fin, log_file = 'output')
    elif fp_style == "pwmat" :
        forward_files = ['atom.config', 'etot.input'] + fp_pp_files
        backward_files = ['REPORT', 'OUT.MLMD', 'output']
        run_fp_inner(iter_index, jdata, mdata, forward_files, backward_files, _pwmat_check_fin, log_file = 'output')
    else :
        raise RuntimeError ("unsupported fp style")


def post_fp_check_fail(iter_index,
                       jdata,
                       rfailed = None) :
    ratio_failed =  rfailed if rfailed else jdata.get('ratio_failed',0.05)
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return
    # check fail according to tag_failure
    fp_failed_tags = glob.glob(os.path.join(work_path, 'task.*', 'tag_failure*'))
    fp_failed_tasks = [os.path.dirname(ii) for ii in fp_failed_tags]
    fp_failed_tasks = list(set(fp_failed_tasks))

    ntask = len(fp_tasks)
    nfail = len(fp_failed_tasks)
    rfail = float(nfail) / float(ntask)
    dlog.info("failed tasks: %6d in %6d  %6.2f %% " % (nfail, ntask, rfail * 100.))
    if rfail > ratio_failed:
       raise RuntimeError("find too many unsuccessfully terminated jobs")


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
                   raise RuntimeError('invalid setting of use_ele_temp ' + str(use_ele_temp))

    rfail=float(icount)/float(tcount)
    dlog.info("failed frame: %6d in %6d  %6.2f %% " % (icount, tcount, rfail * 100.))

    if rfail>ratio_failed:
       raise RuntimeError("find too many unsuccessfully terminated jobs. Too many FP tasks are not converged. Please check your input parameters (e.g. INCAR) or configuration (e.g. POSCAR) in directories \'iter.*.*/02.fp/task.*.*/.\'")


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

def post_fp_abacus_pw_scf (iter_index,
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
        sys_output = glob.glob(os.path.join(work_path, "task.%s.*"%ss))
        sys_input = glob.glob(os.path.join(work_path, "task.%s.*/INPUT"%ss))
        sys_output.sort()
        sys_input.sort()

        flag=True
        for ii,oo in zip(sys_input,sys_output) :
            if flag:
                _sys = dpdata.LabeledSystem(oo, fmt = 'abacus/scf', type_map = jdata['type_map'])
                if len(_sys)>0:
                   all_sys=_sys
                   flag=False
                else:
                   pass
            else:
                _sys = dpdata.LabeledSystem(oo, fmt = 'abacus/scf', type_map = jdata['type_map'])
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
            # TODO : UnboundLocalError sometimes occurs and I cannot figure it out.
            try:
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
            except UnboundLocalError as e:
                pass
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


def post_fp_pwmat (iter_index,
                  jdata,
                  rfailed=None):

    ratio_failed =  rfailed if rfailed else jdata.get('ratio_failed',0.05)
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

    tcount=0
    icount=0
    for ss in system_index :
        sys_output = glob.glob(os.path.join(work_path, "task.%s.*/OUT.MLMD"%ss))
        sys_output.sort()
        tcount += len(sys_output)
        all_sys = None
        for oo in sys_output :
            _sys = dpdata.LabeledSystem(oo, type_map = jdata['type_map'])
            if len(_sys) == 1:
                if all_sys is None:
                    all_sys = _sys
                else:
                    all_sys.append(_sys)
            else:
                icount+=1
        if all_sys is not None:
           sys_data_path = os.path.join(work_path, 'data.%s'%ss)
           all_sys.to_deepmd_raw(sys_data_path)
           all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))
    dlog.info("failed frame number: %s "%icount)
    dlog.info("total frame number: %s "%tcount)
    reff=icount/tcount
    dlog.info('ratio of failed frame:  {:.2%}'.format(reff))

    if reff>ratio_failed:
       raise RuntimeError("find too many unsuccessfully terminated jobs")


def post_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']
    post_fp_check_fail(iter_index, jdata)
    if fp_style == "vasp" :
        post_fp_vasp(iter_index, jdata)
    elif fp_style == "pwscf" :
        post_fp_pwscf(iter_index, jdata)
    elif fp_style == "abacus/scf":
        post_fp_abacus_pw_scf(iter_index, jdata)
    elif fp_style == "siesta":
        post_fp_siesta(iter_index, jdata)
    elif fp_style == 'gaussian' :
        post_fp_gaussian(iter_index, jdata)
    elif fp_style == 'cp2k' :
        post_fp_cp2k(iter_index, jdata)
    elif fp_style == 'pwmat' :
        post_fp_pwmat(iter_index, jdata)
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
    
    deepmd_version = '1'
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

    if mdata.get('handlers', None):
        if mdata['handlers'].get('smtp', None):
            que = queue.Queue(-1)
            queue_handler = logging.handlers.QueueHandler(que)
            smtp_handler = logging.handlers.SMTPHandler(**mdata['handlers']['smtp'])
            listener = logging.handlers.QueueListener(que, smtp_handler)
            dlog.addHandler(queue_handler)
            listener.start()

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
            sepline("{} {}".format(iter_name, task_name),'-')
            if   jj == 0 :
                log_iter ("make_train", ii, jj)
                make_train (ii, jdata, mdata)
            elif jj == 1 :
                log_iter ("run_train", ii, jj)
                mdata  = decide_train_machine(mdata)
                run_train  (ii, jdata, mdata)
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
                run_model_devi (ii, jdata, mdata)

            elif jj == 5 :
                log_iter ("post_model_devi", ii, jj)
                post_model_devi (ii, jdata, mdata)
            elif jj == 6 :
                log_iter ("make_fp", ii, jj)
                make_fp (ii, jdata, mdata)
            elif jj == 7 :
                log_iter ("run_fp", ii, jj)
                mdata = decide_fp_machine(mdata)
                run_fp (ii, jdata, mdata)
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

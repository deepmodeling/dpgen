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
import itertools
import copy
import dpdata
import numpy as np
import subprocess as sp
import scipy.constants as pc
from collections import Counter
from collections.abc import Iterable
from distutils.version import LooseVersion
from typing import List
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
from dpgen.generator.lib.utils import symlink_user_forward_files
from dpgen.generator.lib.lammps import make_lammps_input, get_dumped_forces, get_all_dumped_forces
from dpgen.generator.lib.make_calypso import _make_model_devi_native_calypso,_make_model_devi_buffet
from dpgen.generator.lib.run_calypso import gen_structures,analysis,run_calypso_model_devi
from dpgen.generator.lib.parse_calypso import _parse_calypso_input,_parse_calypso_dis_mtx
from dpgen.generator.lib.vasp import write_incar_dict
from dpgen.generator.lib.vasp import make_vasp_incar_user_dict
from dpgen.generator.lib.vasp import incar_upper
from dpgen.generator.lib.pwscf import make_pwscf_input
from dpgen.generator.lib.abacus_scf import make_abacus_scf_stru, make_abacus_scf_input, make_abacus_scf_kpt
from dpgen.generator.lib.abacus_scf import get_abacus_input_parameters
#from dpgen.generator.lib.pwscf import cvt_1frame
from dpgen.generator.lib.pwmat import make_pwmat_input_dict
from dpgen.generator.lib.pwmat import write_input_dict
from dpgen.generator.lib.pwmat import make_pwmat_input_user_dict
from dpgen.generator.lib.pwmat import input_upper
from dpgen.generator.lib.siesta import make_siesta_input
from dpgen.generator.lib.gaussian import make_gaussian_input, take_cluster
from dpgen.generator.lib.cp2k import make_cp2k_input, make_cp2k_input_from_external, make_cp2k_xyz
from dpgen.generator.lib.ele_temp import NBandsEsti
from dpgen.remote.decide_machine import convert_mdata
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, make_dispatcher, make_submission
from dpgen.util import sepline, expand_sys_str, normalize
from dpgen import ROOT_PATH
from pymatgen.io.vasp import Incar,Kpoints,Potcar
from dpgen.auto_test.lib.vasp import make_kspacing_kpoints
from .arginfo import run_jdata_arginfo


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
# for calypso 
calypso_run_opt_name = 'gen_stru_analy'
calypso_model_devi_name = 'model_devi_results'
calypso_run_model_devi_file = os.path.join(ROOT_PATH,'generator/lib/calypso_run_model_devi.py')
check_outcar_file = os.path.join(ROOT_PATH,'generator/lib/calypso_check_outcar.py')
run_opt_file = os.path.join(ROOT_PATH,'generator/lib/calypso_run_opt.py')

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
    # check the number of collected data
    sys_data = glob.glob(os.path.join(fp_path, "data.*"))
    empty_sys = []
    for ii in sys_data :
        nframe = 0
        sys_paths = expand_sys_str(ii)
        for single_sys in sys_paths:
            sys = dpdata.LabeledSystem(os.path.join(single_sys), fmt = 'deepmd/npy')
            nframe += len(sys)
        empty_sys.append(nframe < max_v)
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


# def dump_to_poscar(dump, poscar, type_map, fmt = "lammps/dump") :
#    sys = dpdata.System(dump, fmt = fmt, type_map = type_map)
#    sys.to_vasp_poscar(poscar)

def dump_to_deepmd_raw(dump, deepmd_raw, type_map, fmt='gromacs/gro', charge=None):
    system = dpdata.System(dump, fmt = fmt, type_map = type_map)
    system.to_deepmd_raw(deepmd_raw)
    if charge is not None:
        with open(os.path.join(deepmd_raw, "charge"), 'w') as f:
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
        training_reuse_stop_batch = 400000
        
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

    model_devi_engine = jdata.get('model_devi_engine', "lammps")
    if iter_index > 0 and _check_empty_iter(iter_index-1, fp_task_min) :
        log_task('prev data is empty, copy prev model')
        copy_model(numb_models, iter_index-1, iter_index)
        return
    elif model_devi_engine != 'calypso' and iter_index > 0 and _check_skip_train(model_devi_jobs[iter_index-1]):
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
        sys_paths = expand_sys_str(os.path.join(init_data_prefix, ii))
        for single_sys in sys_paths:
            init_data_sys.append(os.path.normpath(os.path.join('..', 'data.init', ii, os.path.relpath(single_sys, os.path.join(init_data_prefix, ii)))))
            init_batch_size.append(detect_batch_size(ss, single_sys))
    old_range = None
    if iter_index > 0 :
        for ii in range(iter_index) :
            if ii == iter_index - 1:
                old_range = len(init_data_sys)
            fp_path = os.path.join(make_iter_name(ii), fp_name)
            fp_data_sys = glob.glob(os.path.join(fp_path, "data.*"))
            if model_devi_engine == 'calypso':
                _modd_path = os.path.join(make_iter_name(ii), model_devi_name, calypso_model_devi_name)
                sys_list = glob.glob(os.path.join(_modd_path, "*.structures"))
                sys_batch_size = ["auto" for aa in range(len(sys_list))]
            for jj in fp_data_sys :
                sys_idx = int(jj.split('.')[-1])
                sys_paths = expand_sys_str(jj)
                nframes = 0
                for sys_single in sys_paths:
                    nframes += dpdata.LabeledSystem(sys_single, fmt="deepmd/npy").get_nframes()
                if nframes < fp_task_min :
                    log_task('nframes (%d) in data sys %s is too small, skip' % (nframes, jj))
                    continue
                for sys_single in sys_paths:
                    init_data_sys.append(os.path.normpath(os.path.join('..', 'data.iters', sys_single)))
                    init_batch_size.append(detect_batch_size(sys_batch_size[sys_idx], sys_single))
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
        raise RuntimeError("DP-GEN currently only supports for DeePMD-kit 1.x or 2.x version!" )
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
            # HDF5 path contains #
            if not (os.path.isdir(jj) if "#" not in jj else os.path.isfile(jj.split("#")[0])):
                raise RuntimeError ("data sys %s does not exists, cwd is %s" % (jj, os.getcwd()))
        os.chdir(cwd)
        # set random seed for each model
        if LooseVersion(mdata["deepmd_version"]) >= LooseVersion('1') and LooseVersion(mdata["deepmd_version"]) < LooseVersion('3'):
            # 1.x
            if jinput['model']['descriptor']['type'] == 'hybrid':
                for desc in jinput['model']['descriptor']['list']:
                    desc['seed'] = random.randrange(sys.maxsize) % (2**32)
            elif jinput['model']['descriptor']['type'] == 'loc_frame':
                pass
            else:
                jinput['model']['descriptor']['seed'] = random.randrange(sys.maxsize) % (2**32)
            jinput['model']['fitting_net']['seed'] = random.randrange(sys.maxsize) % (2**32)
            if 'type_embedding' in jinput['model']:
                jinput['model']['type_embedding']['seed'] = random.randrange(sys.maxsize) % (2**32)
            jinput['training']['seed'] = random.randrange(sys.maxsize) % (2**32)
        else:
            raise RuntimeError("DP-GEN currently only supports for DeePMD-kit 1.x or 2.x version!" )
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
    # Copy user defined forward files
    symlink_user_forward_files(mdata=mdata, task_type="train", work_path=work_path)
    


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
        # check if h5 file
        format = 'deepmd/npy' if "#" not in system else 'deepmd/hdf5'
        s = dpdata.LabeledSystem(system, fmt=format)
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
        if jdata.get("dp_compress", False):
            commands.append("%s compress" % train_command)
    else:
        raise RuntimeError("DP-GEN currently only supports for DeePMD-kit 1.x or 2.x version!" )

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
    if jdata.get("dp_compress", False):
        backward_files.append('frozen_model_compressed.pb')
    init_data_sys_ = jdata['init_data_sys']
    init_data_sys = []
    for ii in init_data_sys_ :
        init_data_sys.append(os.path.join('data.init', ii))
    trans_comm_data = []
    cwd = os.getcwd()
    os.chdir(work_path)
    fp_data = glob.glob(os.path.join('data.iters', 'iter.*', '02.fp', 'data.*'))
    for ii in itertools.chain(init_data_sys, fp_data) :
        sys_paths = expand_sys_str(ii)
        for single_sys in sys_paths:
            if "#" not in single_sys:
                trans_comm_data += glob.glob(os.path.join(single_sys, 'set.*'))
                trans_comm_data += glob.glob(os.path.join(single_sys, 'type*.raw'))
                trans_comm_data += glob.glob(os.path.join(single_sys, 'nopbc'))
            else:
                # H5 file
                trans_comm_data.append(single_sys.split("#")[0])
    # remove duplicated files
    trans_comm_data = list(set(trans_comm_data))
    os.chdir(cwd)

    try:
        train_group_size = mdata['train_group_size']
    except Exception:
        train_group_size = 1

    api_version = mdata.get('api_version', '0.9')
    
    user_forward_files = mdata.get("train" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files += mdata.get("train" + "_user_backward_files", [])
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
        if not jdata.get("dp_compress", False):
            model_name = 'frozen_model.pb'
        else:
            model_name = 'frozen_model_compressed.pb'
        task_file = os.path.join(train_task_fmt % ii, model_name)
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
    if model_devi_engine != 'calypso':
        if (iter_index >= len(model_devi_jobs)) :
            return False
    else:
        # mode 1: generate structures according to the user-provided input.dat file, so calypso_input_path and model_devi_max_iter are needed
        run_mode = 1
        if "calypso_input_path" in jdata:
            try:
                maxiter = jdata.get('model_devi_max_iter')
            except KeyError:
                raise KeyError('calypso_input_path key exists so you should provide model_devi_max_iter key to control the max iter number')
        # mode 2: control each iteration to generate structures in specific way by providing model_devi_jobs key
        else:
            try:
                maxiter = max(model_devi_jobs[-1].get('times'))
                run_mode = 2
            except KeyError:
                raise KeyError('did not find model_devi_jobs["times"] key')
        if (iter_index > maxiter) :
            dlog.info(f'iter_index is {iter_index} and maxiter is {maxiter}')
            return False

    if "sys_configs_prefix" in jdata:
        sys_configs = []
        for sys_list in jdata["sys_configs"]:
            #assert (isinstance(sys_list, list) ), "Currently only support type list for sys in 'sys_conifgs' "
            temp_sys_list = [os.path.join(jdata["sys_configs_prefix"], sys) for sys in sys_list]
            sys_configs.append(temp_sys_list)
    else:
        sys_configs = jdata['sys_configs']
    shuffle_poscar = jdata.get('shuffle_poscar', False)

    if model_devi_engine != 'calypso':
        cur_job = model_devi_jobs[iter_index]
        sys_idx = expand_idx(cur_job['sys_idx'])
    else:
        cur_job = {'model_devi_engine':'calypso','input.dat':'user_provided'}
        sys_idx = []

    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")
    conf_systems = []
    for idx in sys_idx :
        cur_systems = []
        ss = sys_configs[idx]
        for ii in ss :
            cur_systems += sorted(glob.glob(ii))
        # cur_systems should not be sorted, as we may add specific constrict to the similutions 
        #cur_systems.sort()
        cur_systems = [os.path.abspath(ii) for ii in cur_systems]
        conf_systems.append (cur_systems)

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = sorted(glob.glob(os.path.join(train_path, "graph*pb")))
    work_path = os.path.join(iter_name, model_devi_name)
    create_path(work_path)
    if model_devi_engine == 'calypso':
        _calypso_run_opt_path = os.path.join(work_path,calypso_run_opt_name)
        calypso_model_devi_path = os.path.join(work_path,calypso_model_devi_name)
        create_path(calypso_model_devi_path)
        # run model devi script
        calypso_run_model_devi_script = os.path.join(calypso_model_devi_path,'calypso_run_model_devi.py')
        shutil.copyfile(calypso_run_model_devi_file,calypso_run_model_devi_script)
        # Create work path list
        calypso_run_opt_path = []

        # mode 1: generate structures according to the user-provided input.dat file,
        # so calypso_input_path and model_devi_max_iter are needed
        if run_mode == 1:
            if jdata.get('vsc', False) and len(jdata.get('type_map')) > 1:
                # [input.dat.Li.250, input.dat.Li.300]
                one_ele_inputdat_list = glob.glob(
                        f"{jdata.get('calypso_input_path')}/input.dat.{jdata.get('type_map')[0]}.*"
                        )
                if len(one_ele_inputdat_list) == 0:
                    number_of_pressure = 1
                else: 
                    number_of_pressure = len(list(set(one_ele_inputdat_list)))

                # calypso_run_opt_path = ['gen_struc_analy.000','gen_struc_analy.001']
                for temp_idx in range(number_of_pressure):
                    calypso_run_opt_path.append('%s.%03d'%(_calypso_run_opt_path, temp_idx))
            elif not jdata.get('vsc', False):
                calypso_run_opt_path.append('%s.%03d'%(_calypso_run_opt_path, 0))
                        
        # mode 2: control each iteration to generate structures in specific way 
        # by providing model_devi_jobs key
        elif run_mode == 2:
            for iiidx, jobbs in enumerate(model_devi_jobs):
                if iter_index in jobbs.get('times'):
                    cur_job = model_devi_jobs[iiidx]
                    
            pressures_list = cur_job.get('PSTRESS', [0.0001])
            for temp_idx in range(len(pressures_list)):
                calypso_run_opt_path.append('%s.%03d'%(_calypso_run_opt_path, temp_idx))
        # to different directory
        # calypso_run_opt_path = ['gen_struc_analy.000','gen_struc_analy.001','gen_struc_analy.002',]
        for temp_calypso_run_opt_path in calypso_run_opt_path:
            create_path(temp_calypso_run_opt_path)
            # run confs opt script
            run_opt_script = os.path.join(temp_calypso_run_opt_path,'calypso_run_opt.py')
            shutil.copyfile(run_opt_file,run_opt_script)
            # check outcar script
            check_outcar_script = os.path.join(temp_calypso_run_opt_path,'check_outcar.py')
            shutil.copyfile(check_outcar_file,check_outcar_script)

    for mm in models :
        model_name = os.path.basename(mm)
        if model_devi_engine != 'calypso':
            os.symlink(mm, os.path.join(work_path, model_name))
        else:
            for temp_calypso_run_opt_path in calypso_run_opt_path:
                models_path = os.path.join(temp_calypso_run_opt_path, model_name)
                if not os.path.exists(models_path):
                    os.symlink(mm, models_path)

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
            elif model_devi_engine == "amber":
                # Jinzhe's specific Amber version
                conf_name = make_model_devi_conf_name(sys_idx[sys_counter], conf_counter)
                rst7_name = conf_name + '.rst7'
                # link restart file
                os.symlink(cc, os.path.join(conf_path, rst7_name))
            conf_counter += 1
        sys_counter += 1

    input_mode = "native"
    if "calypso_input_path" in jdata:
        input_mode = "buffet"
    if "template" in cur_job:
        input_mode = "revise_template"
    use_plm = jdata.get('model_devi_plumed', False)
    use_plm_path = jdata.get('model_devi_plumed_path', False)
    if input_mode == "native":
        if model_devi_engine == "lammps":
            _make_model_devi_native(iter_index, jdata, mdata, conf_systems)
        elif model_devi_engine == "gromacs":
            _make_model_devi_native_gromacs(iter_index, jdata, mdata, conf_systems)
        elif model_devi_engine == "amber":
            _make_model_devi_amber(iter_index, jdata, mdata, conf_systems)
        elif model_devi_engine == "calypso":
            _make_model_devi_native_calypso(iter_index,model_devi_jobs, calypso_run_opt_path)  # generate input.dat automatic in each iter
        else:
            raise RuntimeError("unknown model_devi engine", model_devi_engine)
    elif input_mode == "revise_template":
        _make_model_devi_revmat(iter_index, jdata, mdata, conf_systems)
    elif input_mode == "buffet":
        _make_model_devi_buffet(jdata,calypso_run_opt_path)  # generate confs according to the input.dat provided
    else:
        raise RuntimeError('unknown model_devi input mode', input_mode)
    #Copy user defined forward_files
    symlink_user_forward_files(mdata=mdata, task_type="model_devi", work_path=work_path)
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
    models = sorted(glob.glob(os.path.join(train_path, "graph*pb")))
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
                model_devi_merge_traj = jdata.get('model_devi_merge_traj', False)
                if not model_devi_merge_traj :
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
                # only revise the line "pair_style deepmd" if the user has not written the full line (checked by then length of the line)
                template_has_pair_deepmd=1
                for line_idx,line_context in enumerate(lmp_lines):
                    if (line_context[0] != "#") and ("pair_style" in line_context) and ("deepmd" in line_context):
                        template_has_pair_deepmd=0
                        template_pair_deepmd_idx=line_idx
                if template_has_pair_deepmd == 0:
                    if LooseVersion(deepmd_version) < LooseVersion('1'):
                        if len(lmp_lines[template_pair_deepmd_idx].split()) !=  (len(models) + len(["pair_style","deepmd","10", "model_devi.out"])):
                            lmp_lines = revise_lmp_input_model(lmp_lines, task_model_list, trj_freq, deepmd_version = deepmd_version)
                    else:
                        if len(lmp_lines[template_pair_deepmd_idx].split()) != (len(models) + len(["pair_style","deepmd","out_freq", "10", "out_file", "model_devi.out"])):
                            lmp_lines = revise_lmp_input_model(lmp_lines, task_model_list, trj_freq, deepmd_version = deepmd_version)
                #use revise_lmp_input_model to raise error message if "part_style" or "deepmd" not found
                else:
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
                    model_devi_merge_traj = jdata.get('model_devi_merge_traj', False)
                    if not model_devi_merge_traj :
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
    try:
        from gromacs.fileformats.mdp import MDP
    except ImportError as e:
        raise RuntimeError("GromacsWrapper>=0.8.0 is needed for DP-GEN + Gromacs.") from e
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
    lambdas = cur_job.get("lambdas", [1.0])
    temps = cur_job.get("temps", [298.0])

    for ll in lambdas:
        assert (ll >= 0.0 and ll <= 1.0), "Lambda should be in [0,1]"

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

def _make_model_devi_amber(iter_index: int, jdata: dict, mdata: dict, conf_systems: list):
    """Make amber's MD inputs.

    Parameters
    ----------
    iter_index : int
        iter index
    jdata : dict
        run parameters. The following parameters will be used in this method:
            model_devi_jobs : list[dict]
                The list including the dict for information of each cycle:
                    sys_idx : list[int]
                        list of systems to run
                    trj_freq : int
                        freq to dump trajectory
            low_level : str
                low level method
            cutoff : float
                cutoff radius of the DPRc model
            parm7_prefix : str
                The path prefix to AMBER PARM7 files
            parm7 : list[str]
                List of paths to AMBER PARM7 files. Each file maps to a system.
            mdin_prefix : str 
                The path prefix to AMBER mdin files
            mdin : list[str]
                List of paths to AMBER mdin files. Each files maps to a system.
                The following keywords will be replaced by the actual value:
                    @freq@ : freq to dump trajectory
                    @nstlim@ : total time step to run
                    @qm_region@ : AMBER mask of the QM region
                    @qm_theory@ : The QM theory, such as DFTB2
                    @qm_charge@ : The total charge of the QM theory, such as -2
                    @rcut@ : cutoff radius of the DPRc model
                    @GRAPH_FILE0@, @GRAPH_FILE1@, ... : graph files
            qm_region : list[str]
                AMBER mask of the QM region. Each mask maps to a system.
            qm_charge : list[int]
                Charge of the QM region. Each charge maps to a system.
            nsteps : list[int]
                The number of steps to run. Each number maps to a system.
            r : list[list[float]] or list[list[list[float]]]
                Constrict values for the enhanced sampling. The first dimension maps to systems.
                The second dimension maps to confs in each system. The third dimension is the
                constrict value. It can be a single float for 1D or list of floats for nD.
            disang_prefix : str
                The path prefix to disang prefix.
            disang : list[str]
                List of paths to AMBER disang files. Each file maps to a sytem.
                The keyword RVAL will be replaced by the constrict values, or RVAL1, RVAL2, ...
                for an nD system.
    mdata : dict
        machine parameters. Nothing will be used in this method.
    conf_systems : list
        conf systems

    References
    ----------
    .. [1] Development of Range-Corrected Deep Learning Potentials for Fast, Accurate Quantum
       Mechanical/Molecular Mechanical Simulations of Chemical Reactions in Solution, 
       Jinzhe Zeng, Timothy J. Giese, Şölen Ekesan, and Darrin M. York, Journal of Chemical
       Theory and Computation 2021 17 (11), 6993-7009 

    inputs: restart (coords), param, mdin, graph, disang (optional)

    """
    model_devi_jobs = jdata['model_devi_jobs']
    if (iter_index >= len(model_devi_jobs)) :
        return False
    cur_job = model_devi_jobs[iter_index]   
    sys_idx = expand_idx(cur_job['sys_idx'])
    if (len(sys_idx) != len(list(set(sys_idx)))) :
        raise RuntimeError("system index should be uniq")

    iter_name = make_iter_name(iter_index)
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    work_path = os.path.join(iter_name, model_devi_name)
    # parm7 - list
    parm7 = jdata['parm7']
    parm7_prefix = jdata.get("parm7_prefix", "")
    parm7 = [os.path.join(parm7_prefix, pp) for pp in parm7]

    # link parm file
    for ii, pp in enumerate(parm7):
        os.symlink(pp, os.path.join(work_path, 'qmmm%d.parm7'%ii))
    # TODO: consider writing input in json instead of a given file
    # mdin 
    mdin = jdata['mdin']
    mdin_prefix = jdata.get("mdin_prefix", "")
    mdin = [os.path.join(mdin_prefix, pp) for pp in mdin]

    qm_region = jdata['qm_region']
    qm_charge = jdata['qm_charge']
    nsteps = jdata['nsteps']

    for ii, pp in enumerate(mdin):
        with open(pp) as f, open(os.path.join(work_path, 'init%d.mdin'%ii), 'w') as fw:
            mdin_str = f.read()
            # freq, nstlim, qm_region, qm_theory, qm_charge, rcut, graph
            mdin_str = mdin_str.replace("@freq@", str(cur_job.get('trj_freq', 50))) \
                               .replace("@nstlim@", str(nsteps[ii])) \
                               .replace("@qm_region@", qm_region[ii]) \
                               .replace("@qm_charge@", str(qm_charge[ii])) \
                               .replace("@qm_theory@", jdata['low_level']) \
                               .replace("@rcut@", str(jdata['cutoff']))
            models = sorted(glob.glob(os.path.join(train_path, "graph.*.pb")))
            task_model_list = []
            for ii in models:
                task_model_list.append(os.path.join('..', os.path.basename(ii)))
            # graph
            for jj, mm in enumerate(task_model_list):
                # replace graph
                mdin_str = mdin_str.replace("@GRAPH_FILE%d@" % jj, mm)
            fw.write(mdin_str)
    # disang - list
    disang = jdata['disang']
    disang_prefix = jdata.get("disang_prefix", "")
    disang = [os.path.join(disang_prefix, pp) for pp in disang]

    for sys_counter, ss in enumerate(conf_systems):
        for idx_cc, cc in enumerate(ss) :
            task_counter = idx_cc
            conf_counter = idx_cc

            task_name = make_model_devi_task_name(sys_idx[sys_counter], task_counter)
            conf_name = make_model_devi_conf_name(sys_idx[sys_counter], conf_counter)
            task_path = os.path.join(work_path, task_name)
            # create task path
            create_path(task_path)
            # link restart file
            loc_conf_name = 'init.rst7'
            os.symlink(os.path.join(os.path.join('..','confs'), conf_name + ".rst7"),
                    os.path.join(task_path, loc_conf_name) )
            cwd_ = os.getcwd()
            # chdir to task path
            os.chdir(task_path)
            
            # reaction coordinates of umbrella sampling
            # TODO: maybe consider a better name instead of `r`?
            if 'r' in jdata:
                r=jdata['r'][sys_idx[sys_counter]][conf_counter]
                # r can either be a float or a list of float (for 2D coordinates)
                if not isinstance(r, Iterable) or isinstance(r, str):
                    r = [r]
                # disang file should include RVAL, RVAL2, ...
                with open(disang[sys_idx[sys_counter]]) as f, open('TEMPLATE.disang', 'w') as fw:
                    tl = f.read()
                    for ii, rr in enumerate(r):
                        if isinstance(rr, Iterable) and not isinstance(rr, str):
                            raise RuntimeError("rr should not be iterable! sys: %d rr: %s r: %s" % (sys_idx[sys_counter], str(rr), str(r)))
                        tl = tl.replace("RVAL"+str(ii+1), str(rr))
                    if len(r) == 1:
                        tl = tl.replace("RVAL", str(r[0]))
                    fw.write(tl)

            with open('job.json', 'w') as fp:
                json.dump(cur_job, fp, indent = 4)
            os.chdir(cwd_)

def run_md_model_devi (iter_index,
                       jdata,
                       mdata) :

    #rmdlog.info("This module has been run !")
    model_devi_exec = mdata['model_devi_command']

    model_devi_group_size = mdata['model_devi_group_size']
    model_devi_resources = mdata['model_devi_resources']
    use_plm = jdata.get('model_devi_plumed', False)
    use_plm_path = jdata.get('model_devi_plumed_path', False)
    model_devi_merge_traj = jdata.get('model_devi_merge_traj', False)

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
        command = "{ if [ ! -f dpgen.restart.10000 ]; then %s -i input.lammps -v restart 0; else %s -i input.lammps -v restart 1; fi }" % (model_devi_exec, model_devi_exec)
        command = "/bin/sh -c '%s'" % command
        commands = [command]
        
        forward_files = ['conf.lmp', 'input.lammps']
        backward_files = ['model_devi.out', 'model_devi.log']
        if model_devi_merge_traj :
            backward_files += ['all.lammpstrj']
        else :
            forward_files += ['traj']     
            backward_files += ['traj']          

        if use_plm:
            forward_files += ['input.plumed']
           # backward_files += ['output.plumed']
            backward_files += ['output.plumed','COLVAR']
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

        command = "%s grompp -f %s -p %s -c %s -o %s -maxwarn %d" % (model_devi_exec, mdp_filename, topol_filename, conf_filename, deffnm, maxwarn)
        command += "&& %s mdrun -deffnm %s -cpi" %(model_devi_exec, deffnm)
        if ndx_filename:
            command += f"&& echo -e \"{grp_name}\\n{grp_name}\\n\" | {model_devi_exec} trjconv -s {ref_filename} -f {deffnm}.trr -n {ndx_filename} -o {traj_filename} -pbc mol -ur compact -center"
        else:
            command += "&& echo -e \"%s\\n%s\\n\" | %s trjconv -s %s -f %s.trr -o %s -pbc mol -ur compact -center" % (grp_name, grp_name, model_devi_exec, ref_filename, deffnm, traj_filename)
        command += "&& if [ ! -d traj ]; then \n mkdir traj; fi\n"
        command += f"python -c \"import dpdata;system = dpdata.System('{traj_filename}', fmt='gromacs/gro'); [system.to_gromacs_gro('traj/%d.gromacstrj' % (i * {trj_freq}), frame_idx=i) for i in range(system.get_nframes())]; system.to_deepmd_npy('traj_deepmd')\""
        command += f"&& dp model-devi -m ../graph.000.pb ../graph.001.pb ../graph.002.pb ../graph.003.pb -s traj_deepmd -o model_devi.out -f {trj_freq}"
        commands = [command]

        forward_files = [mdp_filename, topol_filename, conf_filename, index_filename, ref_filename, type_filename, "input.json", "job.json" ]
        if ndx_filename: forward_files.append(ndx_filename)
        backward_files = ["%s.tpr" % deffnm, "%s.log" %deffnm , traj_filename, 'model_devi.out', "traj", "traj_deepmd" ]
    elif model_devi_engine == "amber":
        commands = [(
            "TASK=$(basename $(pwd)) && "
            "SYS1=${TASK:5:3} && "
            "SYS=$((10#$SYS1)) && "
        )+ model_devi_exec + (
            " -O -p ../qmmm$SYS.parm7 -c init.rst7 -i ../init$SYS.mdin -o rc.mdout -r rc.rst7 -x rc.nc -inf rc.mdinfo -ref init.rst7"
        )]
        forward_files = ['init.rst7', 'TEMPLATE.disang']
        backward_files = ['rc.mdout', 'rc.nc', 'rc.rst7', 'TEMPLATE.dumpave']
        model_names.extend(["qmmm*.parm7", "init*.mdin"])

    cwd = os.getcwd()

    user_forward_files = mdata.get("model_devi" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files += mdata.get("model_devi" + "_user_backward_files", [])
    api_version = mdata.get('api_version', '0.9')
    if(len(run_tasks) == 0): 
        raise RuntimeError("run_tasks for model_devi should not be empty! Please check your files.") 
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

def run_model_devi(iter_index,jdata,mdata):

    model_devi_engine = jdata.get("model_devi_engine", "lammps")
    if model_devi_engine != "calypso":
        run_md_model_devi(iter_index,jdata,mdata)
    else:
        run_calypso_model_devi(iter_index,jdata,mdata)
    
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
        #
        elif key == 'wrap_ratio':
            ratio=[sys['cells'][0][1][0]/sys['cells'][0][0][0],sys['cells'][0][2][1]/sys['cells'][0][1][1],sys['cells'][0][2][0]/sys['cells'][0][0][0]]
            if np.max(np.abs(ratio)) > float(value):
                is_bad = True
        elif key == 'tilt_ratio':
            ratio=[sys['cells'][0][1][0]/sys['cells'][0][1][1],sys['cells'][0][2][1]/sys['cells'][0][2][2],sys['cells'][0][2][0]/sys['cells'][0][2][2]]
            if np.max(np.abs(ratio)) > float(value):
                is_bad= True
        else:
            raise RuntimeError('unknow key', key)
    return is_bad

def _read_model_devi_file(
        task_path : str,
        model_devi_f_avg_relative : bool = False,
        model_devi_merge_traj : bool = False
):
    model_devi = np.loadtxt(os.path.join(task_path, 'model_devi.out'))
    if model_devi_f_avg_relative :
        if(model_devi_merge_traj is True) : 
            all_traj = os.path.join(task_path, 'all.lammpstrj')
            all_f = get_all_dumped_forces(all_traj)
        else :
            trajs = glob.glob(os.path.join(task_path, 'traj', '*.lammpstrj'))
            all_f = []
            for ii in trajs:
                all_f.append(get_dumped_forces(ii))     

        all_f = np.array(all_f)
        all_f = all_f.reshape([-1,3])
        avg_f = np.sqrt(np.average(np.sum(np.square(all_f), axis = 1)))
        model_devi[:,4:7] = model_devi[:,4:7] / avg_f
        np.savetxt(os.path.join(task_path, 'model_devi_avgf.out'), model_devi, fmt='%16.6e')
    return model_devi


def _select_by_model_devi_standard(
        modd_system_task: List[str],
        f_trust_lo : float,
        f_trust_hi : float,
        v_trust_lo : float,
        v_trust_hi : float,
        cluster_cutoff : float, 
        model_devi_engine : str,
        model_devi_skip : int = 0,
        model_devi_f_avg_relative : bool = False,
        model_devi_merge_traj : bool = False, 
        detailed_report_make_fp : bool = True,
):
    if model_devi_engine == 'calypso':
        iter_name = modd_system_task[0].split('/')[0]
        _work_path = os.path.join(iter_name, model_devi_name)
        # calypso_run_opt_path = os.path.join(_work_path,calypso_run_opt_name)
        calypso_run_opt_path = glob.glob('%s/%s.*'%(_work_path, calypso_run_opt_name))[0]
        numofspecies = _parse_calypso_input('NumberOfSpecies',calypso_run_opt_path)
        min_dis = _parse_calypso_dis_mtx(numofspecies,calypso_run_opt_path)
    fp_candidate = []
    if detailed_report_make_fp:
        fp_rest_accurate = []
        fp_rest_failed = []
    cc = 0
    counter = Counter()
    counter['candidate'] = 0
    counter['failed'] = 0
    counter['accurate'] = 0
    for tt in modd_system_task :
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            all_conf = _read_model_devi_file(tt, model_devi_f_avg_relative, model_devi_merge_traj)

            if all_conf.shape == (7,):
                all_conf = all_conf.reshape(1,all_conf.shape[0])
            elif model_devi_engine == 'calypso' and all_conf.shape == (8,):
                all_conf = all_conf.reshape(1,all_conf.shape[0])
            for ii in range(all_conf.shape[0]) :
                if all_conf[ii][0] < model_devi_skip :
                    continue
                cc = int(all_conf[ii][0])
                if cluster_cutoff is None:
                    if model_devi_engine == 'calypso':
                        if float(all_conf[ii][-1]) <= float(min_dis):
                            if detailed_report_make_fp:
                                fp_rest_failed.append([tt, cc])
                            counter['failed'] += 1
                            continue
                    if (all_conf[ii][1] < v_trust_hi and all_conf[ii][1] >= v_trust_lo) or \
                       (all_conf[ii][4] < f_trust_hi and all_conf[ii][4] >= f_trust_lo) :
                        fp_candidate.append([tt, cc])
                        counter['candidate'] += 1
                    elif (all_conf[ii][1] >= v_trust_hi ) or (all_conf[ii][4] >= f_trust_hi ):
                        if detailed_report_make_fp:
                            fp_rest_failed.append([tt, cc])
                        counter['failed'] += 1
                    elif (all_conf[ii][1] < v_trust_lo and all_conf[ii][4] < f_trust_lo ):
                        if detailed_report_make_fp:
                            fp_rest_accurate.append([tt, cc])
                        counter['accurate'] += 1
                    else :
                        if model_devi_engine == 'calypso':
                            dlog.info('ase opt traj %s frame %d with f devi %f does not belong to either accurate, candidiate and failed '% (tt, ii, all_conf[ii][4]))
                        else:
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
    
    return fp_rest_accurate, fp_candidate, fp_rest_failed, counter



def _select_by_model_devi_adaptive_trust_low(
        modd_system_task: List[str],
        f_trust_hi : float,
        numb_candi_f : int,
        perc_candi_f : float,
        v_trust_hi : float,
        numb_candi_v : int,
        perc_candi_v : float,
        model_devi_skip : int = 0,
        model_devi_f_avg_relative : bool = False,
        model_devi_merge_traj : bool = False, 
):
    """
    modd_system_task    model deviation tasks belonging to one system
    f_trust_hi
    numb_candi_f        number of candidate due to the f model deviation
    perc_candi_f        percentage of candidate due to the f model deviation
    v_trust_hi
    numb_candi_v        number of candidate due to the v model deviation
    perc_candi_v        percentage of candidate due to the v model deviation
    model_devi_skip
    
    returns
    accur               the accurate set
    candi               the candidate set
    failed              the failed set
    counter             counters, number of elements in the sets
    f_trust_lo          adapted trust level of f
    v_trust_lo          adapted trust level of v
    """
    idx_v = 1
    idx_f = 4
    accur = set()
    candi = set()
    failed = []
    coll_v = []
    coll_f = []
    for tt in modd_system_task:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model_devi = np.loadtxt(os.path.join(tt, 'model_devi.out'))
            model_devi = _read_model_devi_file(tt, model_devi_f_avg_relative, model_devi_merge_traj)
            for ii in range(model_devi.shape[0]) :
                if model_devi[ii][0] < model_devi_skip :
                    continue
                cc = int(model_devi[ii][0])
                # tt: name of task folder
                # cc: time step of the frame
                md_v = model_devi[ii][idx_v]
                md_f = model_devi[ii][idx_f]
                if md_f > f_trust_hi or md_v > v_trust_hi:
                    failed.append([tt, cc])
                else:
                    coll_v.append([model_devi[ii][idx_v], tt, cc])
                    coll_f.append([model_devi[ii][idx_f], tt, cc])
                    # now accur takes all non-failed frames,
                    # will be substracted by candidate lat er
                    accur.add((tt, cc))
    # sort
    coll_v.sort()
    coll_f.sort()
    assert(len(coll_v) == len(coll_f))
    # calcuate numbers
    numb_candi_v = max(numb_candi_v, int(perc_candi_v * 0.01 * len(coll_v)))
    numb_candi_f = max(numb_candi_f, int(perc_candi_f * 0.01 * len(coll_f)))
    # adjust number of candidate
    if len(coll_v) < numb_candi_v:
        numb_candi_v = len(coll_v)
    if len(coll_f) < numb_candi_f:
        numb_candi_f = len(coll_f)
    # compute trust lo
    if numb_candi_v == 0:
        v_trust_lo = v_trust_hi
    else:
        v_trust_lo = coll_v[-numb_candi_v][0]
    if numb_candi_f == 0:
        f_trust_lo = f_trust_hi
    else:
        f_trust_lo = coll_f[-numb_candi_f][0]        
    # add to candidate set
    for ii in range(len(coll_v) - numb_candi_v, len(coll_v)):
        candi.add(tuple(coll_v[ii][1:]))
    for ii in range(len(coll_f) - numb_candi_f, len(coll_f)):
        candi.add(tuple(coll_f[ii][1:]))
    # accurate set is substracted by the candidate set
    accur = accur - candi
    # convert to list
    candi = [list(ii) for ii in candi]
    accur = [list(ii) for ii in accur]
    # counters
    counter = Counter()
    counter['candidate'] = len(candi)
    counter['failed'] = len(failed)
    counter['accurate'] = len(accur)

    return accur, candi, failed, counter, f_trust_lo, v_trust_lo
    

def _make_fp_vasp_inner (iter_index, 
                         modd_path,
                         work_path,
                         model_devi_skip,
                         v_trust_lo,
                         v_trust_hi,
                         f_trust_lo,
                         f_trust_hi,
                         fp_task_min,
                         fp_task_max,
                         fp_link_files,
                         type_map,
                         jdata):
    """
    iter_index          int             iter index
    modd_path           string          path of model devi
    work_path           string          path of fp
    fp_task_max         int             max number of tasks
    fp_link_files       [string]        linked files for fp, POTCAR for example
    fp_params           map             parameters for fp
    """

    # --------------------------------------------------------------------------------------------------------------------------------------
    model_devi_engine = jdata.get('model_devi_engine', 'lammps')
    if model_devi_engine == 'calypso':
        iter_name = work_path.split('/')[0]
        _work_path = os.path.join(iter_name, model_devi_name)
        # calypso_run_opt_path = os.path.join(_work_path,calypso_run_opt_name)
        calypso_run_opt_path = glob.glob('%s/%s.*'%(_work_path, calypso_run_opt_name))[0]
        numofspecies = _parse_calypso_input('NumberOfSpecies',calypso_run_opt_path)
        min_dis = _parse_calypso_dis_mtx(numofspecies,calypso_run_opt_path)

        calypso_total_fp_num = 300
        modd_path = os.path.join(modd_path,calypso_model_devi_name)
        model_devi_skip = -1
        with open(os.path.join(modd_path,'Model_Devi.out'),'r') as summfile:
            summary = np.loadtxt(summfile)
        summaryfmax = summary[:,-4]
        dis  = summary[:,-1]
        acc  = np.where((summaryfmax <= f_trust_lo) & (dis > float(min_dis)))
        fail = np.where((summaryfmax >  f_trust_hi) | (dis <= float(min_dis)))
        nnan = np.where(np.isnan(summaryfmax))

        acc_num  = len(acc[0])
        fail_num = len(fail[0])
        nan_num  = len(nnan[0])
        tot = len(summaryfmax) - nan_num
        candi_num = tot - acc_num - fail_num
        dlog.info("summary  accurate_ratio: {0:8.4f}%  candidata_ratio: {1:8.4f}%  failed_ratio: {2:8.4f}%  in {3:d} structures".format(
                   acc_num*100/tot,candi_num*100/tot,fail_num*100/tot,tot ))
    # --------------------------------------------------------------------------------------------------------------------------------------

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

    cluster_cutoff = jdata.get('cluster_cutoff', None)
    model_devi_adapt_trust_lo = jdata.get('model_devi_adapt_trust_lo', False)
    model_devi_f_avg_relative = jdata.get('model_devi_f_avg_relative', False)
    model_devi_merge_traj = jdata.get('model_devi_merge_traj', False)
    # skip save *.out if detailed_report_make_fp is False, default is True
    detailed_report_make_fp = jdata.get("detailed_report_make_fp", True)
    # skip bad box criteria
    skip_bad_box = jdata.get('fp_skip_bad_box')
    # skip discrete structure in cluster
    fp_cluster_vacuum = jdata.get('fp_cluster_vacuum',None)

    def _trust_limitation_check(sys_idx, lim):
        if isinstance(lim, list):
            sys_lim = lim[sys_idx]
        elif isinstance(lim, dict):
            sys_lim = lim[str(sys_idx)]
        else:
            sys_lim = lim
        return sys_lim

    for ss in system_index:
        modd_system_glob = os.path.join(modd_path, 'task.' + ss + '.*')
        modd_system_task = glob.glob(modd_system_glob)
        modd_system_task.sort()
        if model_devi_engine in ('lammps', 'gromacs', 'calypso'):
            # convert global trust limitations to local ones
            f_trust_lo_sys = _trust_limitation_check(int(ss), f_trust_lo)
            f_trust_hi_sys = _trust_limitation_check(int(ss), f_trust_hi)
            v_trust_lo_sys = _trust_limitation_check(int(ss), v_trust_lo)
            v_trust_hi_sys = _trust_limitation_check(int(ss), v_trust_hi)

            # assumed e -> v
            if not model_devi_adapt_trust_lo:
                fp_rest_accurate, fp_candidate, fp_rest_failed, counter \
                    =  _select_by_model_devi_standard(
                        modd_system_task,
                        f_trust_lo_sys, f_trust_hi_sys,
                        v_trust_lo_sys, v_trust_hi_sys,
                        cluster_cutoff, 
                        model_devi_engine,
                        model_devi_skip,
                        model_devi_f_avg_relative = model_devi_f_avg_relative,
                        model_devi_merge_traj = model_devi_merge_traj, 
                        detailed_report_make_fp = detailed_report_make_fp,
                    )
            else:
                numb_candi_f = jdata.get('model_devi_numb_candi_f', 10)
                numb_candi_v = jdata.get('model_devi_numb_candi_v', 0)
                perc_candi_f = jdata.get('model_devi_perc_candi_f', 0.)
                perc_candi_v = jdata.get('model_devi_perc_candi_v', 0.)
                fp_rest_accurate, fp_candidate, fp_rest_failed, counter, f_trust_lo_ad, v_trust_lo_ad \
                    =  _select_by_model_devi_adaptive_trust_low(
                        modd_system_task,
                        f_trust_hi_sys, numb_candi_f, perc_candi_f,
                        v_trust_hi_sys, numb_candi_v, perc_candi_v,
                        model_devi_skip = model_devi_skip,
                        model_devi_f_avg_relative = model_devi_f_avg_relative,
                        model_devi_merge_traj = model_devi_merge_traj, 
                    )
                dlog.info("system {0:s} {1:9s} : f_trust_lo {2:6.3f}   v_trust_lo {3:6.3f}".format(ss, 'adapted', f_trust_lo_ad, v_trust_lo_ad))
        elif model_devi_engine == "amber":
            counter = Counter()
            counter['candidate'] = 0
            counter['failed'] = 0
            counter['accurate'] = 0
            fp_rest_accurate = []
            fp_candidate = []
            fp_rest_failed = []
            for tt in modd_system_task :
                cc = 0
                with open(os.path.join(tt, "rc.mdout")) as f:
                    skip_first = False
                    first_active = True
                    for line in f:
                        if line.startswith("     ntx     =       1"):
                            skip_first = True
                        if line.startswith("Active learning frame written with max. frc. std.:"):
                            if skip_first and first_active:
                                first_active = False
                                continue
                            model_devi = float(line.split()[-2]) * dpdata.unit.EnergyConversion("kcal_mol", "eV").value()
                            if model_devi < f_trust_lo:
                                # accurate
                                if detailed_report_make_fp:
                                    fp_rest_accurate.append([tt, cc])
                                counter['accurate'] += 1
                            elif model_devi > f_trust_hi:
                                # failed
                                if detailed_report_make_fp:
                                    fp_rest_failed.append([tt, cc])
                                counter['failed'] += 1
                            else:
                                # candidate
                                fp_candidate.append([tt, cc])
                                counter['candidate'] += 1
                            cc += 1

        else:
            raise RuntimeError('unknown model_devi_engine', model_devi_engine)

        # print a report
        fp_sum = sum(counter.values())

        if fp_sum == 0:
            dlog.info('system {0:s} has no fp task, maybe the model devi is nan %'.format(ss))
            continue
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
        # ----------------------------------------------------------------------------
        if model_devi_engine == 'calypso':
            calypso_intend_fp_num_temp = (len(fp_candidate)/candi_num)*calypso_total_fp_num
            if calypso_intend_fp_num_temp < 1:
                calypso_intend_fp_num = 1
            else:
                calypso_intend_fp_num = int(calypso_intend_fp_num_temp)
        # ----------------------------------------------------------------------------
        numb_task = min(this_fp_task_max, len(fp_candidate))
        if (numb_task < fp_task_min):
            numb_task = 0

        # ---------------------------------------------------------------------------- 
        if (model_devi_engine == 'calypso' and len(jdata.get('type_map')) == 1) or \
           (model_devi_engine == 'calypso' and len(jdata.get('type_map')) > 1 and candi_num <= calypso_total_fp_num):
            numb_task = min(this_fp_task_max, len(fp_candidate))
            if (numb_task < fp_task_min):
                numb_task = 0
        elif (model_devi_engine == 'calypso' and len(jdata.get('type_map')) > 1 and candi_num > calypso_total_fp_num):
            numb_task = calypso_intend_fp_num
            if (len(fp_candidate) < numb_task):
                numb_task = 0
        # ----------------------------------------------------------------------------
        dlog.info("system {0:s} accurate_ratio: {1:8.4f}    thresholds: {2:6.4f} and {3:6.4f}   eff. task min and max {4:4d} {5:4d}   number of fp tasks: {6:6d}".format(ss, accurate_ratio, fp_accurate_soft_threshold, fp_accurate_threshold, fp_task_min, this_fp_task_max, numb_task))
        # make fp tasks
        
        # read all.lammpstrj, save in all_sys for each system_index
        all_sys = []
        trj_freq = None
        if model_devi_merge_traj :
            for ii in modd_system_task :
                all_traj = os.path.join(ii, 'all.lammpstrj')
                all_sys_per_task = dpdata.System(all_traj, fmt = 'lammps/dump', type_map = type_map)
                all_sys.append(all_sys_per_task)
            model_devi_jobs = jdata['model_devi_jobs']
            cur_job = model_devi_jobs[iter_index]
            trj_freq = int(_get_param_alias(cur_job, ['t_freq', 'trj_freq', 'traj_freq']))
        
        count_bad_box = 0
        count_bad_cluster = 0
        fp_candidate = sorted(fp_candidate[:numb_task])

        for cc in range(numb_task) :
            tt = fp_candidate[cc][0]
            ii = fp_candidate[cc][1]
            ss = os.path.basename(tt).split('.')[1]
            conf_name = os.path.join(tt, "traj")
            conf_sys = None
            if model_devi_engine == "lammps":
                if model_devi_merge_traj :
                    conf_sys = all_sys[int(os.path.basename(tt).split('.')[-1])][int(int(ii) / trj_freq)]
                else :
                    conf_name = os.path.join(conf_name, str(ii) + '.lammpstrj')
                ffmt = 'lammps/dump'
            elif model_devi_engine == "gromacs":
                conf_name = os.path.join(conf_name, str(ii) + '.gromacstrj')
                ffmt = 'lammps/dump'
            elif model_devi_engine == "amber":
                conf_name = os.path.join(tt, "rc.nc")
                rst_name = os.path.abspath(os.path.join(tt, "init.rst7"))
            elif model_devi_engine == "calypso":
                conf_name = os.path.join(conf_name, str(ii) + '.poscar')
                ffmt = 'vasp/poscar'
            else:
                raise RuntimeError("unknown model_devi engine", model_devi_engine)
            conf_name = os.path.abspath(conf_name)
            if skip_bad_box is not None:
                skip = check_bad_box(conf_name, skip_bad_box, fmt=ffmt)
                if skip:
                    count_bad_box += 1
                    continue

            if fp_cluster_vacuum is not None:
                assert fp_cluster_vacuum >0
                skip_cluster = check_cluster(conf_name, fp_cluster_vacuum)
                if skip_cluster:
                    count_bad_cluster +=1
                    continue

            if model_devi_engine != 'calypso':
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
                if model_devi_engine == "lammps":                   
                    if model_devi_merge_traj:
                        conf_sys.to("lammps/lmp", "conf.dump")
                    else: 
                        os.symlink(os.path.relpath(conf_name), 'conf.dump')
                    os.symlink(os.path.relpath(job_name), 'job.json')
                elif model_devi_engine == "gromacs":
                    os.symlink(os.path.relpath(conf_name), 'conf.dump')
                    os.symlink(os.path.relpath(job_name), 'job.json')
                elif model_devi_engine == "amber":
                    # read and write with ase
                    from ase.io.netcdftrajectory import NetCDFTrajectory, write_netcdftrajectory
                    if cc > 0 and tt == fp_candidate[cc-1][0]:
                        # same MD task, use the same file
                        pass
                    else:
                        # not the same file
                        if cc > 0:
                            # close the old file
                            netcdftraj.close()
                        netcdftraj = NetCDFTrajectory(conf_name)
                    # write nc file
                    write_netcdftrajectory('rc.nc', netcdftraj[ii])
                    if cc >= numb_task - 1:
                        netcdftraj.close()
                    # link restart since it's necessary to start Amber
                    os.symlink(os.path.relpath(rst_name), 'init.rst7')
                    os.symlink(os.path.relpath(job_name), 'job.json')
                elif model_devi_engine == "calypso":
                    os.symlink(os.path.relpath(conf_name), 'POSCAR')
                    fjob = open('job.json','w+')
                    fjob.write('{"model_devi_engine":"calypso"}')
                    fjob.close()
                    #os.system('touch job.json')
                else:
                    raise RuntimeError('unknown model_devi_engine', model_devi_engine)
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
    if model_devi_engine == 'calypso':
        dlog.info("summary  accurate_ratio: {0:8.4f}%  candidata_ratio: {1:8.4f}%  failed_ratio: {2:8.4f}%  in {3:d} structures".format( acc_num*100/tot,candi_num*100/tot,fail_num*100/tot,tot ))
    if cluster_cutoff is None:
        cwd = os.getcwd()
        for idx, task in enumerate(fp_tasks):
            os.chdir(task)
            if model_devi_engine == "lammps":
                sys = None
                if model_devi_merge_traj:
                    sys = dpdata.System('conf.dump', fmt = "lammps/lmp", type_map = type_map)
                else :
                    sys = dpdata.System('conf.dump', fmt = "lammps/dump", type_map = type_map)
                sys.to_vasp_poscar('POSCAR')
                # dump to poscar 

                if charges_map:
                    warnings.warn('"sys_charges" keyword only support for gromacs engine now.')
            elif model_devi_engine == "gromacs":
                # dump_to_poscar('conf.dump', 'POSCAR', type_map, fmt = "gromacs/gro")
                if charges_map:
                    dump_to_deepmd_raw('conf.dump', 'deepmd.raw', type_map, fmt='gromacs/gro', charge=charges_recorder[idx])
                else:
                    dump_to_deepmd_raw('conf.dump', 'deepmd.raw', type_map, fmt='gromacs/gro', charge=None)
            elif model_devi_engine in ("amber", 'calypso'):
                pass
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

def _link_fp_abacus_orb_descript (iter_index,
                      jdata) :
    # assume lcao orbital files, numerical descrptors and model for dpks are all in fp_pp_path.
    fp_pp_path = jdata['fp_pp_path']

    fp_orb_files = jdata['fp_orb_files']
    assert(os.path.exists(fp_pp_path))
    fp_dpks_descriptor = None
    fp_dpks_model = None
    if "fp_dpks_descriptor" in jdata:
        fp_dpks_descriptor = jdata["fp_dpks_descriptor"]
    if "user_fp_params" in jdata:
        if "deepks_model" in jdata["user_fp_params"]:
            fp_dpks_model = jdata["user_fp_params"]["deepks_model"]

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
        for jj in fp_orb_files:
            orb_file = os.path.join(fp_pp_path, jj)
            os.symlink(orb_file, jj)
        if fp_dpks_descriptor is not None:
            descrptor = os.path.join(fp_pp_path, fp_dpks_descriptor)
            os.symlink(descrptor, fp_dpks_descriptor)
        if fp_dpks_model is not None:
            dpks_model = os.path.join(fp_pp_path, fp_dpks_model)
            os.symlink(dpks_model, fp_dpks_model)
        os.chdir(cwd)

def _make_fp_vasp_configs(iter_index,
                          jdata):
    fp_task_max = jdata['fp_task_max']
    model_devi_skip = jdata['model_devi_skip']
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
    else:
        cur_job = {}
    # support iteration dependent trust levels
    v_trust_lo = cur_job.get('model_devi_v_trust_lo', jdata.get('model_devi_v_trust_lo', 1e10))
    v_trust_hi = cur_job.get('model_devi_v_trust_hi', jdata.get('model_devi_v_trust_hi', 1e10))
    if cur_job.get('model_devi_f_trust_lo') is not None:
        f_trust_lo = cur_job.get('model_devi_f_trust_lo')
    else:
        f_trust_lo = jdata['model_devi_f_trust_lo']
    if cur_job.get('model_devi_f_trust_hi') is not None:
        f_trust_hi = cur_job.get('model_devi_f_trust_hi')
    else:
        f_trust_hi = jdata['model_devi_f_trust_hi']

    # make configs
    fp_tasks = _make_fp_vasp_inner(iter_index, 
                                   modd_path, work_path,
                                   model_devi_skip,
                                   v_trust_lo, v_trust_hi,
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
        sys_data['atom_masses'] = []
        pps = []
        for iii in sys_data['atom_names']:
            sys_data['atom_masses'].append(jdata['mass_map'][jdata['type_map'].index(iii)])
            pps.append(fp_pp_files[jdata['type_map'].index(iii)])
        ret = make_pwscf_input(sys_data, pps, fp_params, user_input = user_input)
        with open('input', 'w') as fp:
            fp.write(ret)
        os.chdir(cwd)
    # link pp files
    _link_fp_vasp_pp(iter_index, jdata)

def make_fp_abacus_scf(iter_index,
                  jdata) :
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    if len(fp_tasks) == 0 :
        return
    # make abacus/pw/scf input
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_pp_files = jdata['fp_pp_files']
    fp_orb_files = None
    fp_dpks_descriptor = None
    # get paramters for writting INPUT file
    fp_params = {}
    if 'user_fp_params' in jdata.keys() :
        fp_params = jdata['user_fp_params']
        # for lcao 
        if 'basis_type' in fp_params:
            if fp_params['basis_type'] == 'lcao':
                assert('fp_orb_files' in jdata and type(jdata['fp_orb_files']) == list and len(jdata['fp_orb_files']) == len(fp_pp_files))
                fp_orb_files = jdata['fp_orb_files']
        dpks_out_labels = fp_params.get('deepks_out_labels',0)
        dpks_scf = fp_params.get('deepks_scf',0)
        if dpks_out_labels or dpks_scf:
            assert('fp_dpks_descriptor' in jdata and type(jdata['fp_dpks_descriptor']) == str)
            fp_dpks_descriptor = jdata['fp_dpks_descriptor']
        #user_input = True
        ret_input = make_abacus_scf_input(fp_params)
    elif 'fp_incar' in jdata.keys():
        fp_input_path = jdata['fp_incar']
        assert(os.path.exists(fp_input_path))
        fp_input_path = os.path.abspath(fp_input_path)
        fp_params = get_abacus_input_parameters(fp_input_path)
        ret_input = make_abacus_scf_input(fp_params)
    else:
        raise RuntimeError("Set 'user_fp_params' or 'fp_incar' in json file to make INPUT of ABACUS")
    # get paramters for writting KPT file
    if 'kspacing' not in fp_params.keys():
        if 'gamma_only' in fp_params.keys():
            if fp_params["gamma_only"]==1:
                gamma_param = {"k_points":[1,1,1,0,0,0]}
                ret_kpt = make_abacus_scf_kpt(gamma_param)
            else:
                if 'k_points' in jdata.keys() :
                    ret_kpt = make_abacus_scf_kpt(jdata)
                elif 'fp_kpt_file' in jdata.keys():
                    fp_kpt_path = jdata['fp_kpt_file']
                    assert(os.path.exists(fp_kpt_path))
                    fp_kpt_path = os.path.abspath(fp_kpt_path)
                    fk = open(fp_kpt_path)
                    ret_kpt = fk.read()
                    fk.close()
                else:
                    raise RuntimeError("Cannot find any k-points information")
        else:
            if 'k_points' in jdata.keys() :
                ret_kpt = make_abacus_scf_kpt(jdata)
            elif 'fp_kpt_file' in jdata.keys():
                fp_kpt_path = jdata['fp_kpt_file']
                assert(os.path.exists(fp_kpt_path))
                fp_kpt_path = os.path.abspath(fp_kpt_path)
                fk = open(fp_kpt_path)
                ret_kpt = fk.read()
                fk.close()
            else:
                gamma_param = {"k_points":[1,1,1,0,0,0]}
                ret_kpt = make_abacus_scf_kpt(gamma_param)
                warnings.warn("Cannot find k-points information, gamma_only will be generated.")

    cwd = os.getcwd()
    for ii in fp_tasks:
        os.chdir(ii)
        sys_data = dpdata.System('POSCAR').data
        if 'mass_map' in jdata:
            sys_data['atom_masses'] = jdata['mass_map']
        with open('INPUT', 'w') as fp:
            fp.write(ret_input)
        if 'kspacing' not in fp_params.keys():
            with open("KPT", "w") as fp:
                fp.write(ret_kpt)
        ret_stru = make_abacus_scf_stru(sys_data, fp_pp_files, fp_orb_files, fp_dpks_descriptor, fp_params)
        with open("STRU", "w") as fp:
            fp.write(ret_stru)

        os.chdir(cwd)
    # link pp files
    _link_fp_vasp_pp(iter_index, jdata)
    if 'basis_type' in fp_params:
            if fp_params['basis_type'] == 'lcao':
                _link_fp_abacus_orb_descript(iter_index, jdata)


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

def make_fp_amber_diff(iter_index: int, jdata: dict):
    """Run amber twice to calculate high-level and low-level potential,
    and then generate difference between them.

    Besides AMBER, one needs to install `dpamber` package, which is avaiable at
    https://github.com/njzjz/dpamber

    Currently, it should be used with the AMBER model_devi driver.

    Parameters
    ----------
    iter_index : int
        iter index
    jdata : dict
        Run parameters. The following parameters are used in this method:
            mdin_prefix : str 
                The path prefix to AMBER mdin files
            qm_region : list[str]
                AMBER mask of the QM region. Each mask maps to a system.
            qm_charge : list[int]
                Charge of the QM region. Each charge maps to a system.
            high_level : str
                high level method
            low_level : str
                low level method
            fp_params : dict
                This parameters includes:
                    high_level_mdin : str
                        High-level AMBER mdin file. %qm_theory%, %qm_region%,
                        and %qm_charge% will be replace.
                    low_level_mdin : str
                        Low-level AMBER mdin file. %qm_theory%, %qm_region%,
                        and %qm_charge% will be replace.
            parm7_prefix : str
                The path prefix to AMBER PARM7 files
            parm7 : list[str]
                List of paths to AMBER PARM7 files. Each file maps to a system.
    
    References
    ----------
    .. [1] Development of Range-Corrected Deep Learning Potentials for Fast, Accurate Quantum
       Mechanical/Molecular Mechanical Simulations of Chemical Reactions in Solution, 
       Jinzhe Zeng, Timothy J. Giese, Şölen Ekesan, and Darrin M. York, Journal of Chemical
       Theory and Computation 2021 17 (11), 6993-7009 
    """
    # make config
    fp_tasks = _make_fp_vasp_configs(iter_index, jdata)
    # make amber input
    cwd = os.getcwd()
    # link two mdin files and param7
    os.chdir(os.path.join(fp_tasks[0], ".."))
    mdin_prefix = jdata.get('mdin_prefix', '')
    low_level_mdin = jdata['fp_params']['low_level_mdin']
    low_level_mdin = os.path.join(mdin_prefix, low_level_mdin)
    high_level_mdin = jdata['fp_params']['high_level_mdin']
    high_level_mdin = os.path.join(mdin_prefix, high_level_mdin)
    with open(low_level_mdin) as f:
        low_level_mdin_str = f.read()
    with open(high_level_mdin) as f:
        high_level_mdin_str = f.read()

    qm_region = jdata['qm_region']
    high_level = jdata['high_level']
    low_level = jdata['low_level']
    qm_charge = jdata['qm_charge']
    # qm_theory qm_region qm_charge
    for ii, _ in enumerate(qm_region):
        mdin_new_str = low_level_mdin_str.replace("%qm_theory%", low_level) \
                               .replace("%qm_region%", qm_region[ii]) \
                               .replace("%qm_charge%", str(qm_charge[ii]))
        with open('low_level%d.mdin'%ii, 'w') as f:
            f.write(mdin_new_str)

        mdin_new_str = high_level_mdin_str.replace("%qm_theory%", high_level) \
                               .replace("%qm_region%", qm_region[ii]) \
                               .replace("%qm_charge%", str(qm_charge[ii]))
        with open('high_level%d.mdin'%ii, 'w') as f:
            f.write(mdin_new_str)

    parm7 = jdata['parm7']
    parm7_prefix = jdata.get("parm7_prefix", "")
    parm7 = [os.path.join(parm7_prefix, pp) for pp in parm7]
    for ii, pp in enumerate(parm7):
        os.symlink(pp, "qmmm%d.parm7"%ii)
    
    rst7_prefix = jdata.get("sys_configs_prefix", "")
    for ii, ss in enumerate(jdata['sys_configs']):
        os.symlink(os.path.join(rst7_prefix, ss[0]), "init%d.rst7"%ii)

    with open("qm_region", 'w') as f:
        f.write("\n".join(qm_region))
    os.chdir(cwd)

def make_fp (iter_index,
             jdata,
             mdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        make_fp_vasp(iter_index, jdata)
    elif fp_style == "pwscf" :
        make_fp_pwscf(iter_index, jdata)
    elif fp_style == "abacus" :
        make_fp_abacus_scf(iter_index, jdata)
    elif fp_style == "siesta" :
        make_fp_siesta(iter_index, jdata)
    elif fp_style == "gaussian" :
        make_fp_gaussian(iter_index, jdata)
    elif fp_style == "cp2k" :
        make_fp_cp2k(iter_index, jdata)
    elif fp_style == "pwmat" :
        make_fp_pwmat(iter_index, jdata)
    elif fp_style == "amber/diff":
        make_fp_amber_diff(iter_index, jdata)
    else :
        raise RuntimeError ("unsupported fp style")
    # Copy user defined forward_files
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    symlink_user_forward_files(mdata=mdata, task_type="fp", work_path=work_path)

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

def _abacus_scf_check_fin(ii) :
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
    
    fp_style = jdata['fp_style']
    if fp_style == 'amber/diff':
        # firstly get sys_idx
        fp_command = (
            "TASK=$(basename $(pwd)) && "
            "SYS1=${TASK:5:3} && "
            "SYS=$((10#$SYS1)) && "
            'QM_REGION=$(awk "NR==$SYS+1" ../qm_region) &&'
        ) + fp_command + (
            " -O -p ../qmmm$SYS.parm7 -c ../init$SYS.rst7 -i ../low_level$SYS.mdin -o low_level.mdout -r low_level.rst7 "
            "-x low_level.nc -y rc.nc -frc low_level.mdfrc -inf low_level.mdinfo && "
        ) + fp_command + (
            " -O -p ../qmmm$SYS.parm7 -c ../init$SYS.rst7 -i ../high_level$SYS.mdin -o high_level.mdout -r high_level.rst7 "
            "-x high_level.nc -y rc.nc -frc high_level.mdfrc -inf high_level.mdinfo && "
        ) + (
            "dpamber corr --cutoff %f --parm7_file ../qmmm$SYS.parm7 --nc rc.nc --hl high_level --ll low_level --qm_region \"$QM_REGION\"") % (
               jdata['cutoff'],
        )


    fp_run_tasks = fp_tasks
    # for ii in fp_tasks :
    #     if not check_fin(ii) :
    #         fp_run_tasks.append(ii)
    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]

    user_forward_files = mdata.get("fp" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files += mdata.get("fp" + "_user_backward_files", [])
    
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
    fp_pp_files = jdata.get('fp_pp_files', [])

    if fp_style == "vasp" :
        forward_files = ['POSCAR', 'INCAR', 'POTCAR','KPOINTS']
        backward_files = ['fp.log','OUTCAR','vasprun.xml']
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
    elif fp_style == "abacus":
        fp_params = {}
        if 'user_fp_params' in jdata.keys() :
            fp_params = jdata['user_fp_params']
        elif 'fp_incar' in jdata.keys():
            fp_input_path = jdata['fp_incar']
            assert(os.path.exists(fp_input_path))
            fp_input_path = os.path.abspath(fp_input_path)
            fp_params = get_abacus_input_parameters(fp_input_path)
        forward_files = ["INPUT", "STRU"]
        if 'kspacing' not in fp_params.keys():
            forward_files = ["INPUT","STRU","KPT"]
        forward_files += fp_pp_files
        if "fp_orb_files" in jdata:
            forward_files += jdata["fp_orb_files"]
        if "fp_dpks_descriptor" in jdata:
            forward_files.append(jdata["fp_dpks_descriptor"])
        if "user_fp_params" in jdata:
            if "deepks_model" in jdata["user_fp_params"]:
                forward_files.append(jdata["user_fp_params"]["deepks_model"])
        backward_files = ["output", "OUT.ABACUS"]
        run_fp_inner(iter_index, jdata, mdata,  forward_files, backward_files, _abacus_scf_check_fin, log_file = 'output')
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
    elif fp_style == 'amber/diff':
        forward_files = ['rc.nc']
        backward_files = [
            'low_level.mdfrc', 'low_level.mdout',
            'high_level.mdfrc', 'high_level.mdout',
            'output', 'dataset'
        ]
        forward_common_files = ['low_level*.mdin', 'high_level*.mdin', 'qmmm*.parm7', 'qm_region', 'init*.rst7']
        run_fp_inner(iter_index, jdata, mdata, forward_files, backward_files, None, log_file = 'output',
                     forward_common_files=forward_common_files)
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
    ntask = len(fp_tasks)
    nfail = 0

    # check fail according to the number of collected data
    sys_data = glob.glob(os.path.join(work_path, "data.*"))
    sys_data.sort()
    nframe = 0
    for ii in sys_data :
        sys_paths = expand_sys_str(ii)
        for single_sys in sys_paths:
            sys = dpdata.LabeledSystem(os.path.join(single_sys), fmt = 'deepmd/npy')
            nframe += len(sys)
    nfail = ntask - nframe

    rfail = float(nfail) / float(ntask)
    dlog.info("failed tasks: %6d in %6d  %6.2f %% " % (nfail, ntask, rfail * 100.))
    if rfail > ratio_failed:
       raise RuntimeError("find too many unsuccessfully terminated jobs")


def post_fp_vasp (iter_index,
                  jdata,
                  rfailed=None):

    ratio_failed =  rfailed if rfailed else jdata.get('ratio_failed',0.05)
    model_devi_engine = jdata.get('model_devi_engine', "lammps")
    if model_devi_engine != 'calypso':
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
            except Exception:
                dlog.info('Try to parse from vasprun.xml')
                try:
                   _sys = dpdata.LabeledSystem(oo.replace('OUTCAR','vasprun.xml'), type_map = jdata['type_map'])
                except Exception:
                   _sys = dpdata.LabeledSystem()
                   dlog.info('Failed fp path: %s'%oo.replace('OUTCAR',''))
            if len(_sys) == 1:
                if all_sys is None:
                    all_sys = _sys
                else:
                    all_sys.append(_sys)
                # save ele_temp, if any
                if(os.path.exists(oo.replace('OUTCAR', 'job.json')) ): 
                    with open(oo.replace('OUTCAR', 'job.json')) as fp:
                        job_data = json.load(fp)
                    if 'ele_temp' in job_data:
                        assert(use_ele_temp)
                        ele_temp = job_data['ele_temp']
                        all_te.append(ele_temp)
            elif len(_sys) >= 2:
                raise RuntimeError("The vasp parameter NSW should be set as 1") 
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

    if(tcount == 0) :
        rfail = 0.0
        dlog.info("failed frame: %6d in %6d " % (icount, tcount))
    else :
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

def post_fp_abacus_scf (iter_index,
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
                      jdata,
                      rfailed=None):
                      
    ratio_failed =  rfailed if rfailed else jdata.get('ratio_failed',0.10)
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
    # tcount: num of all fp tasks
    tcount = 0
    # icount: num of converged fp tasks
    icount = 0
    for ss in system_index :
        sys_output = glob.glob(os.path.join(work_path, "task.%s.*/output"%ss))
        sys_output.sort()
        tcount += len(sys_output)
        all_sys = None
        for oo in sys_output :
            _sys = dpdata.LabeledSystem(oo, fmt = 'cp2k/output')
            #_sys.check_type_map(type_map = jdata['type_map'])
            if all_sys is None:
                all_sys = _sys
            else:
                all_sys.append(_sys)


        icount += len(all_sys)
        if all_sys is not None:
            sys_data_path = os.path.join(work_path, 'data.%s'%ss)
            all_sys.to_deepmd_raw(sys_data_path)
            all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output))

    if(tcount == 0) :
        rfail = 0.0
        dlog.info("failed frame: %6d in %6d " % (tcount - icount, tcount))
    else :
        rfail=float(tcount - icount)/float(tcount)
        dlog.info("failed frame: %6d in %6d  %6.2f %% " % (tcount - icount, tcount, rfail * 100.))

    if rfail>ratio_failed:
       raise RuntimeError("find too many unsuccessfully terminated jobs. Too many FP tasks are not converged. Please check your files in directories \'iter.*.*/02.fp/task.*.*/.\'")


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


def post_fp_amber_diff(iter_index, jdata):
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

    for ss in system_index :
        sys_output = glob.glob(os.path.join(work_path, "task.%s.*"%ss))
        sys_output.sort()
        all_sys=dpdata.MultiSystems(type_map=jdata['type_map'])
        for oo in sys_output :
            sys=dpdata.MultiSystems(type_map=jdata['type_map']).from_deepmd_npy(os.path.join(oo, 'dataset'))
            all_sys.append(sys)
        sys_data_path = os.path.join(work_path, 'data.%s'%ss)
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size = len(sys_output), prec=np.float64)

def post_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']
    if fp_style == "vasp" :
        post_fp_vasp(iter_index, jdata)
    elif fp_style == "pwscf" :
        post_fp_pwscf(iter_index, jdata)
    elif fp_style == "abacus":
        post_fp_abacus_scf(iter_index, jdata)
    elif fp_style == "siesta":
        post_fp_siesta(iter_index, jdata)
    elif fp_style == 'gaussian' :
        post_fp_gaussian(iter_index, jdata)
    elif fp_style == 'cp2k' :
        post_fp_cp2k(iter_index, jdata)
    elif fp_style == 'pwmat' :
        post_fp_pwmat(iter_index, jdata)
    elif fp_style == 'amber/diff':
        post_fp_amber_diff(iter_index, jdata)
    else :
        raise RuntimeError ("unsupported fp style")
    post_fp_check_fail(iter_index, jdata)
    # clean traj
    clean_traj = True
    if 'model_devi_clean_traj' in jdata :
        clean_traj = jdata['model_devi_clean_traj']
    modd_path =  None
    if isinstance(clean_traj, bool):
        iter_name = make_iter_name(iter_index)
        if clean_traj:
            modd_path = os.path.join(iter_name, model_devi_name)
    elif isinstance(clean_traj, int):
        clean_index = iter_index - clean_traj
        if clean_index >= 0:
            modd_path = os.path.join(make_iter_name(clean_index), model_devi_name)
    if modd_path is not None:
        md_trajs = glob.glob(os.path.join(modd_path, 'task*/traj'))
        for ii in md_trajs:
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
    except Exception:
       with open (param_file, 'r') as fp :
           jdata = json.load (fp)
       with open (machine_file, 'r') as fp:
           mdata = json.load (fp)

    jdata_arginfo = run_jdata_arginfo()
    jdata = normalize(jdata_arginfo, jdata, strict_check=False)

    update_mass_map(jdata)
        
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
    # Convert mdata
    mdata = convert_mdata(mdata)
    max_tasks = 10000
    numb_task = 9
    record = "record.dpgen"
    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec :
                iter_rec = [int(x) for x in line.split()]
        if len(iter_rec) == 0: 
            raise ValueError("There should not be blank lines in record.dpgen.")
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
                run_model_devi (ii, jdata, mdata)

            elif jj == 5 :
                log_iter ("post_model_devi", ii, jj)
                post_model_devi (ii, jdata, mdata)
            elif jj == 6 :
                log_iter ("make_fp", ii, jj)
                make_fp (ii, jdata, mdata)
            elif jj == 7 :
                log_iter ("run_fp", ii, jj)
                run_fp (ii, jdata, mdata)
            elif jj == 8 :
                log_iter ("post_fp", ii, jj)
                post_fp (ii, jdata)
            else :
                raise RuntimeError ("unknown task %d, something wrong" % jj)
            record_iter (record, ii, jj)


def get_atomic_masses(atom):
    element_names = ['Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon', 'Nitrogen',
                     'Oxygen', 'Fluorine', 'Neon', 'Sodium', 'Magnesium', 'Aluminium', 'Silicon',
                     'Phosphorus', 'Sulfur', 'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
                     'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron', 'Cobalt', 'Nickel',
                     'Copper', 'Zinc', 'Gallium', 'Germanium', 'Arsenic', 'Selenium', 'Bromine',
                     'Krypton', 'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium',
                     'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
                     'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium', 'Iodine', 'Xenon',
                     'Caesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium',
                     'Promethium', 'Samarium', 'Europium', 'Gadolinium', 'Terbium', 'Dysprosium',
                     'Holmium', 'Erbium', 'Thulium', 'Ytterbium', 'Lutetium', 'Hafnium', 'Tantalum',
                     'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
                     'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine', 'Radon', 'Francium',
                     'Radium', 'Actinium', 'Thorium', 'Protactinium', 'Uranium', 'Neptunium',
                     'Plutonium', 'Americium', 'Curium', 'Berkelium', 'Californium', 'Einsteinium',
                     'Fermium', 'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium', 'Dubnium',
                     'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium', 'Darmastadtium', 'Roentgenium',
                     'Copernicium', 'Nihonium', 'Flerovium', 'Moscovium', 'Livermorium', 'Tennessine',
                     'Oganesson']
    chemical_symbols  = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
                         'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
                         'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
                         'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                         'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
                         'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',
                         'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                         'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
                         'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
                         'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    atomic_number  = [ i+1 for i in range(len(chemical_symbols)) ]

    # NIST Standard Reference Database 144
    # URL: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii&isotype=all
    atomic_masses_common = [1.00782503223, 4.00260325413, 7.0160034366, 9.012183065, 11.00930536,
                            12.0, 14.00307400443, 15.99491461957, 18.99840316273, 19.9924401762,
                            22.989769282, 23.985041697, 26.98153853, 27.97692653465, 30.97376199842,
                            31.9720711744, 34.968852682, 39.9623831237, 38.9637064864, 39.962590863,
                            44.95590828, 47.94794198, 50.94395704, 51.94050623, 54.93804391,
                            55.93493633, 58.93319429, 57.93534241, 62.92959772, 63.92914201,
                            68.9255735, 73.921177761, 74.92159457, 79.9165218, 78.9183376, 83.9114977282,
                            84.9117897379, 87.9056125, 88.9058403, 89.9046977, 92.906373, 97.90540482,
                            96.9063667, 101.9043441, 102.905498, 105.9034804, 106.9050916, 113.90336509,
                            114.903878776, 119.90220163, 120.903812, 129.906222748, 126.9044719,
                            131.9041550856, 132.905451961, 137.905247, 138.9063563, 139.9054431,
                            140.9076576, 141.907729, 144.9127559, 151.9197397, 152.921238, 157.9241123,
                            158.9253547, 163.9291819, 164.9303288, 165.9302995, 168.9342179, 173.9388664,
                            174.9407752, 179.946557, 180.9479958, 183.95093092, 186.9557501, 191.961477,
                            192.9629216, 194.9647917, 196.96656879, 201.9706434, 204.9744278, 207.9766525,
                            208.9803991, 208.9824308, 209.9871479, 222.0175782, 223.019736, 226.0254103,
                            227.0277523, 232.0380558, 231.0358842, 238.0507884, 237.0481736, 244.0642053,
                            243.0613813, 247.0703541, 247.0703073, 251.0795886, 252.08298, 257.0951061,
                            258.0984315, 259.10103, 262.10961, 267.12179, 268.12567, 271.13393, 272.13826,
                            270.13429, 276.15159, 281.16451, 280.16514, 285.17712, 284.17873, 289.19042, 
                            288.19274, 293.20449, 292.20746, 294.21392]
    # IUPAC Technical Report
    # doi:10.1515/pac-2015-0305
    atomic_masses_2013 = [1.00784, 4.002602, 6.938, 9.0121831, 10.806, 12.0096, 14.00643, 15.99903,
                          18.99840316, 20.1797, 22.98976928, 24.304, 26.9815385, 28.084, 30.973762,
                          32.059, 35.446, 39.948, 39.0983, 40.078, 44.955908, 47.867, 50.9415, 51.9961,
                          54.938044, 55.845, 58.933194, 58.6934, 63.546, 65.38, 69.723, 72.63, 74.921595,
                          78.971, 79.901, 83.798, 85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, None,
                          101.07, 102.9055, 106.42, 107.8682, 112.414, 114.818, 118.71, 121.76, 127.6,
                          126.90447, 131.293, 132.905452, 137.327, 138.90547, 140.116, 140.90766, 144.242,
                          None, 150.36, 151.964, 157.25, 158.92535, 162.5, 164.93033, 167.259, 168.93422,
                          173.054, 174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084,
                          196.966569, 200.592, 204.382, 207.2, 208.9804, None, None, None, None, None, None,
                          232.0377, 231.03588, 238.02891, None, None, None, None, None, None, None, None,
                          None, None, None, None, None, None, None, None, None, None, None, None, None,
                          None, None, None, None, None]
    # IUPAC Technical Report
    # doi:10.1515/pac-2019-0603
    atomic_masses_2021 = [1.00784, 4.002602, 6.938, 9.0121831, 10.806, 12.0096, 14.00643, 15.99903,
                          18.99840316, 20.1797, 22.98976928, 24.304, 26.9815384, 28.084, 30.973762,
                          32.059, 35.446, 39.792, 39.0983, 40.078, 44.955907, 47.867, 50.9415, 51.9961,
                          54.938043, 55.845, 58.933194, 58.6934, 63.546, 65.38, 69.723, 72.63, 74.921595,
                          78.971, 79.901, 83.798, 85.4678, 87.62, 88.905838, 91.224, 92.90637, 95.95,
                          None, 101.07, 102.90549, 106.42, 107.8682, 112.414, 114.818, 118.71, 121.76,
                          127.6, 126.90447, 131.293, 132.905452, 137.327, 138.90547, 140.116, 140.90766,
                          144.242, None, 150.36, 151.964, 157.25, 158.925354, 162.5, 164.930329, 167.259,
                          168.934219, 173.045, 174.9668, 178.486, 180.94788, 183.84, 186.207, 190.23,
                          192.217, 195.084, 196.96657, 200.592, 204.382, 206.14, 208.9804, None, None,
                          None, None, None, None, 232.0377, 231.03588, 238.02891, None, None, None,
                          None, None, None, None, None, None, None, None, None, None, None, None, None,
                          None, None, None, None, None, None, None, None, None, None]

    atomic_masses = [atomic_masses_common[n] if i is None else i for n,i in enumerate(atomic_masses_2021)]

    if atom in element_names:
        return atomic_masses[element_names.index(atom)]
    elif atom in chemical_symbols:
        return atomic_masses[chemical_symbols.index(atom)]
    elif atom in atomic_number:
        return atomic_masses[atomic_number.index(atom)]
    else:
        raise RuntimeError('unknown atomic identifier', atom, 'if one want to use isotopes, or non-standard element names, chemical symbols, or atomic number in the type_map list, please customize the mass_map list instead of using "auto".')


def update_mass_map(jdata):
    if jdata['mass_map'] == 'auto':
        jdata['mass_map'] = [get_atomic_masses(i) for i in jdata['type_map']]
        
        
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

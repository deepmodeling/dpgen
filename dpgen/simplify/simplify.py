"""Simplify dataset (minimize the dataset size).

Init:
pick up init data from dataset randomly

Iter:
00: train models (same as generator)
01: calculate model deviations of the rest dataset, pick up data with proper model deviaiton 
02: fp (optional, if the original dataset do not have fp data, same as generator)
"""
import logging
import warnings
import queue
import os
import json
import argparse
import pickle
import glob
import fnmatch
import dpdata
import numpy as np
from typing import Union, List

from dpgen import dlog
from dpgen import SHORT_CMD
from dpgen.util import sepline, expand_sys_str, normalize
from distutils.version import LooseVersion
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, make_dispatcher, make_submission
from dpgen.generator.run import make_train, run_train, post_train, run_fp, post_fp, fp_name, model_devi_name, train_name, train_task_fmt, sys_link_fp_vasp_pp, make_fp_vasp_incar, make_fp_vasp_kp, make_fp_vasp_cp_cvasp, data_system_fmt, model_devi_task_fmt, fp_task_fmt
# TODO: maybe the following functions can be moved to dpgen.util
from dpgen.generator.lib.utils import log_iter, make_iter_name, create_path, record_iter
from dpgen.generator.lib.gaussian import make_gaussian_input
from dpgen.remote.decide_machine import  convert_mdata
from .arginfo import simplify_jdata_arginfo


picked_data_name = "data.picked"
rest_data_name = "data.rest"
accurate_data_name = "data.accurate"
detail_file_name_prefix = "details"
sys_name_fmt = 'sys.' + data_system_fmt
sys_name_pattern = 'sys.[0-9]*[0-9]'


def get_system_cls(jdata):
    if jdata.get("labeled", False):
        return dpdata.LabeledSystem
    return dpdata.System


def get_multi_system(path: Union[str, List[str]], jdata: dict) -> dpdata.MultiSystems:
    """Get MultiSystems from a path or list of paths.

    Both NumPy and HDF5 formats are supported. For details
    of two formats, refer to DeePMD-kit documentation.

    If `labeled` in jdata is True, returns MultiSystems with LabeledSystem.
    Otherwise, returns MultiSystems with System.
    
    Parameters
    ----------
    path : str or list of str
        path or list of paths to the dataset
    jdata : dict
        parameters which may contain `labeled` key

    Returns
    -------
    dpdata.MultiSystems
        MultiSystems with LabeledSystem or System
    """
    system = get_system_cls(jdata)
    if not isinstance(path, (list, tuple)):
        path = [path]
    system_paths = []
    for pp in path:
        system_paths.extend(expand_sys_str(pp))
    systems = dpdata.MultiSystems(
        *[system(s, fmt=('deepmd/npy' if "#" not in s else 'deepmd/hdf5')) for s in system_paths],
        type_map=jdata['type_map'],
    )
    return systems


def init_model(iter_index, jdata, mdata):
    training_init_model = jdata.get('training_init_model', False)
    if not training_init_model:
        return
    iter0_models = []
    training_iter0_model = jdata.get('training_iter0_model_path', [])
    if type(training_iter0_model) == str:
        training_iter0_model = [training_iter0_model]
    for ii in training_iter0_model:            
        model_is = glob.glob(ii)
        model_is.sort()
        iter0_models += [os.path.abspath(ii) for ii in model_is]
    numb_models = jdata['numb_models']
    assert(numb_models == len(iter0_models)), "training_iter0_model_path should be provided, and the number of models should be equal to %d" % numb_models
    work_path = os.path.join(make_iter_name(iter_index), train_name)
    create_path(work_path)
    cwd = os.getcwd()
    for ii in range(len(iter0_models)):
        train_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(train_path)
        os.chdir(train_path)
        ckpt_files = glob.glob(os.path.join(iter0_models[ii], 'model.ckpt*'))
        for jj in ckpt_files:
            os.symlink(jj, os.path.basename(jj))
        os.chdir(cwd)


def init_pick(iter_index, jdata, mdata):
    """pick up init data from dataset randomly"""
    pick_data = jdata['pick_data']
    init_pick_number = jdata['init_pick_number']
    # use MultiSystems with System
    # TODO: support System and LabeledSystem
    # TODO: support other format
    systems = get_multi_system(pick_data, jdata)
    # label the system
    labels = []
    items = systems.systems.items()
    for key, system in items:
        labels.extend([(key, j) for j in range(len(system))])

    # random pick
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    create_path(work_path)
    idx = np.arange(len(labels))
    np.random.shuffle(idx)
    pick_idx = idx[:init_pick_number]
    rest_idx = idx[init_pick_number:]

    # dump the init data
    sys_data_path = os.path.join(work_path, picked_data_name)
    _init_dump_selected_frames(systems, labels, pick_idx, sys_data_path, jdata)

    # dump the rest data
    sys_data_path = os.path.join(work_path, rest_data_name)
    _init_dump_selected_frames(systems, labels, rest_idx, sys_data_path, jdata)


def _init_dump_selected_frames(systems, labels, selc_idx, sys_data_path, jdata):
    selc_systems = dpdata.MultiSystems(type_map=jdata['type_map'])
    for j in selc_idx:
        sys_name, sys_id = labels[j]
        selc_systems.append(systems[sys_name][sys_id])
    selc_systems.to_deepmd_raw(sys_data_path)
    selc_systems.to_deepmd_npy(sys_data_path, set_size=selc_idx.size)


def make_model_devi(iter_index, jdata, mdata):
    """calculate the model deviation of the rest idx"""
    pick_data = jdata['pick_data']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    create_path(work_path)
    # link the model
    train_path = os.path.join(iter_name, train_name)
    train_path = os.path.abspath(train_path)
    models = glob.glob(os.path.join(train_path, "graph*pb"))
    for mm in models:
        model_name = os.path.basename(mm)
        os.symlink(mm, os.path.join(work_path, model_name))
    # link the last rest data
    last_iter_name = make_iter_name(iter_index-1)
    rest_data_path = os.path.join(last_iter_name, model_devi_name, rest_data_name)
    if not os.path.exists(rest_data_path):
        return False
    os.symlink(os.path.abspath(rest_data_path), os.path.join(work_path, rest_data_name + ".old"))
    return True


def run_model_devi(iter_index, jdata, mdata):
    """submit dp test tasks"""
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    # generate command
    commands = []
    run_tasks = ["."]
    # get models
    models = glob.glob(os.path.join(work_path, "graph*pb"))
    model_names = [os.path.basename(ii) for ii in models]
    task_model_list = []
    for ii in model_names:
        task_model_list.append(os.path.join('.', ii))
    # models
    commands = []
    detail_file_name = detail_file_name_prefix
    command = "{dp} model-devi -m {model} -s {system} -o {detail_file}".format(
        dp=mdata.get('model_devi_command', 'dp'),
        model=" ".join(task_model_list),
        system=rest_data_name + ".old",
        detail_file=detail_file_name,
    )
    commands = [command]
    # submit
    model_devi_group_size = mdata.get('model_devi_group_size', 1)

    forward_files = [rest_data_name + ".old"]
    backward_files = [detail_file_name]

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


def post_model_devi(iter_index, jdata, mdata):
    """calculate the model deviation"""
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)

    f_trust_lo = jdata['model_devi_f_trust_lo']
    f_trust_hi = jdata['model_devi_f_trust_hi']

    type_map = jdata.get("type_map", [])
    sys_accurate = dpdata.MultiSystems(type_map=type_map)
    sys_candinate = dpdata.MultiSystems(type_map=type_map)
    sys_failed = dpdata.MultiSystems(type_map=type_map)
    
    labeled = jdata.get("labeled", False)
    sys_entire = dpdata.MultiSystems(type_map = type_map).from_deepmd_npy(os.path.join(work_path, rest_data_name + ".old"), labeled=labeled)
    
    detail_file_name = detail_file_name_prefix
    with open(os.path.join(work_path, detail_file_name)) as f:
        for line in f:
            if line.startswith("# data.rest.old"):
                name = (line.split()[1]).split("/")[-1]
            elif line.startswith("#"):
                pass
            else:
                idx = int(line.split()[0])
                f_devi = float(line.split()[4])
                subsys = sys_entire[name][idx]
                if f_trust_lo <= f_devi < f_trust_hi:
                    sys_candinate.append(subsys)
                elif f_devi >= f_trust_hi:
                    sys_failed.append(subsys)
                elif f_devi < f_trust_lo:
                    sys_accurate.append(subsys)
                else:
                    raise RuntimeError('reach a place that should NOT be reached...')

    counter = {"candidate": sys_candinate.get_nframes(), "accurate": sys_accurate.get_nframes(), "failed": sys_failed.get_nframes()}
    fp_sum = sum(counter.values())
    for cc_key, cc_value in counter.items():
        dlog.info("{0:9s} : {1:6d} in {2:6d} {3:6.2f} %".format(cc_key, cc_value, fp_sum, cc_value/fp_sum*100))
    
    if counter['candidate'] == 0 and counter['failed'] > 0:
        raise RuntimeError('no candidate but still have failed cases, stop. You may want to refine the training or to increase the trust level hi')

    # label the candidate system
    labels = []
    items = sys_candinate.systems.items()

    for key, system in items:
        labels.extend([(key, j) for j in range(len(system))])
    # candinate: pick up randomly
    iter_pick_number = jdata['iter_pick_number']
    idx = np.arange(counter['candidate'])
    assert(len(idx) == len(labels))
    np.random.shuffle(idx)
    pick_idx = idx[:iter_pick_number]
    rest_idx = idx[iter_pick_number:]
    if(counter['candidate'] == 0) : 
        dlog.info("no candidate")
    else :
        dlog.info("total candidate {0:6d}   picked {1:6d} ({2:6.2f} %) rest {3:6d} ({4:6.2f} % )".format\
              (counter['candidate'], len(pick_idx), float(len(pick_idx))/counter['candidate']*100., len(rest_idx), float(len(rest_idx))/counter['candidate']*100.))

    # dump the picked candinate data
    picked_systems = dpdata.MultiSystems(type_map = type_map)
    for j in pick_idx:
        sys_name, sys_id = labels[j]
        picked_systems.append(sys_candinate[sys_name][sys_id])
    sys_data_path = os.path.join(work_path, picked_data_name)
    picked_systems.to_deepmd_raw(sys_data_path)
    picked_systems.to_deepmd_npy(sys_data_path, set_size=iter_pick_number)


    # dump the rest data (not picked candinate data and failed data)
    rest_systems = dpdata.MultiSystems(type_map = type_map)
    for j in rest_idx:
        sys_name, sys_id = labels[j]
        rest_systems.append(sys_candinate[sys_name][sys_id])
    rest_systems += sys_failed
    sys_data_path = os.path.join(work_path, rest_data_name)
    rest_systems.to_deepmd_raw(sys_data_path)
    if rest_idx.size:
        rest_systems.to_deepmd_npy(sys_data_path, set_size=rest_idx.size)


    # dump the accurate data -- to another directory
    sys_data_path = os.path.join(work_path, accurate_data_name)
    sys_accurate.to_deepmd_raw(sys_data_path)
    sys_accurate.to_deepmd_npy(sys_data_path, set_size=sys_accurate.get_nframes())


def make_fp_labeled(iter_index, jdata):    
    dlog.info("already labeled, skip make_fp and link data directly")
    pick_data = jdata['pick_data']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    picked_data_path = os.path.join(iter_name, model_devi_name, picked_data_name)
    os.symlink(os.path.abspath(picked_data_path), os.path.abspath(
        os.path.join(work_path, "task." + data_system_fmt % 0)))
    os.symlink(os.path.abspath(picked_data_path), os.path.abspath(
        os.path.join(work_path, "data." + data_system_fmt % 0)))


def make_fp_configs(iter_index, jdata):
    pick_data = jdata['pick_data']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    picked_data_path = os.path.join(iter_name, model_devi_name, picked_data_name)
    systems = get_multi_system(picked_data_path, jdata)
    ii = 0
    jj = 0
    for system in systems:
        for subsys in system:
            task_name = "task." + fp_task_fmt % (ii, jj)
            task_path = os.path.join(work_path, task_name)
            create_path(task_path)
            subsys.to('vasp/poscar', os.path.join(task_path, 'POSCAR'))
            jj += 1
        ii += 1


def make_fp_gaussian(iter_index, jdata):
    work_path = os.path.join(make_iter_name(iter_index), fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    cwd = os.getcwd()
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


def make_fp_vasp(iter_index, jdata):
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


def make_fp_calculation(iter_index, jdata):
    fp_style = jdata['fp_style']
    if fp_style == 'vasp':
        make_fp_vasp(iter_index, jdata)
    elif fp_style == 'gaussian':
        make_fp_gaussian(iter_index, jdata)
    else :
        raise RuntimeError('unsupported fp_style ' + fp_style)


def make_fp(iter_index, jdata, mdata):
    labeled = jdata.get("labeled", False)
    if labeled:
        make_fp_labeled(iter_index, jdata)
    else:
        make_fp_configs(iter_index, jdata)
        make_fp_calculation(iter_index, jdata)


def run_iter(param_file, machine_file):
    """ init (iter 0): init_pick

    tasks (iter > 0):
    00 make_train (same as generator)
    01 run_train (same as generator)
    02 post_train (same as generator)
    03 make_model_devi
    04 run_model_devi
    05 post_model_devi
    06 make_fp
    07 run_fp (same as generator)
    08 post_fp (same as generator)
    """
    # TODO: function of handling input json should be combined as one function
    try:
        import ruamel
        from monty.serialization import loadfn, dumpfn
        warnings.simplefilter(
            'ignore', ruamel.yaml.error.MantissaNoDotYAML1_1Warning)
        jdata = loadfn(param_file)
        mdata = loadfn(machine_file)
    except Exception:
        with open(param_file, 'r') as fp:
            jdata = json.load(fp)
        with open(machine_file, 'r') as fp:
            mdata = json.load(fp)

    jdata_arginfo = simplify_jdata_arginfo()
    jdata = normalize(jdata_arginfo, jdata)

    if mdata.get('handlers', None):
        if mdata['handlers'].get('smtp', None):
            que = queue.Queue(-1)
            queue_handler = logging.handlers.QueueHandler(que)
            smtp_handler = logging.handlers.SMTPHandler(
                **mdata['handlers']['smtp'])
            listener = logging.handlers.QueueListener(que, smtp_handler)
            dlog.addHandler(queue_handler)
            listener.start()
            
    mdata = convert_mdata(mdata)
    max_tasks = 10000
    numb_task = 9
    record = "record.dpgen"
    iter_rec = [0, -1]
    if os.path.isfile(record):
        with open(record) as frec:
            for line in frec:
                iter_rec = [int(x) for x in line.split()]
        dlog.info("continue from iter %03d task %02d" %
                  (iter_rec[0], iter_rec[1]))

    cont = True
    ii = -1
    while cont:
        ii += 1
        iter_name = make_iter_name(ii)
        sepline(iter_name, '=')
        for jj in range(numb_task):
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1]:
                continue
            task_name = "task %02d" % jj
            sepline("{} {}".format(iter_name, task_name), '-')
            jdata['model_devi_jobs'] = [{} for _ in range(ii+1)]
            if ii == 0 and jj < 6:
                if jj == 0:
                    log_iter("init_pick", ii, jj)
                    init_model(ii, jdata, mdata)
                    init_pick(ii, jdata, mdata)
                dlog.info("first iter, skip step 1-5")
            elif jj == 0:
                log_iter("make_train", ii, jj)
                make_train(ii, jdata, mdata)
            elif jj == 1:
                log_iter("run_train", ii, jj)
                #disp = make_dispatcher(mdata['train_machine'])
                run_train(ii, jdata, mdata)
            elif jj == 2:
                log_iter("post_train", ii, jj)
                post_train(ii, jdata, mdata)
            elif jj == 3:
                log_iter("make_model_devi", ii, jj)
                cont = make_model_devi(ii, jdata, mdata)
                if not cont or ii >= jdata.get("stop_iter", ii+1):
                    break
            elif jj == 4:
                log_iter("run_model_devi", ii, jj)
                #disp = make_dispatcher(mdata['model_devi_machine'])
                run_model_devi(ii, jdata, mdata)
            elif jj == 5:
                log_iter("post_model_devi", ii, jj)
                post_model_devi(ii, jdata, mdata)
            elif jj == 6:
                log_iter("make_fp", ii, jj)
                make_fp(ii, jdata, mdata)
            elif jj == 7:
                log_iter("run_fp", ii, jj)
                if jdata.get("labeled", False):
                    dlog.info("already have labeled data, skip run_fp")
                else:
                    #disp = make_dispatcher(mdata['fp_machine'])
                    run_fp(ii, jdata, mdata)
            elif jj == 8:
                log_iter("post_fp", ii, jj)
                if jdata.get("labeled", False):
                    dlog.info("already have labeled data, skip post_fp")
                else:
                    post_fp(ii, jdata)
            else:
                raise RuntimeError("unknown task %d, something wrong" % jj)
            record_iter(record, ii, jj)


def gen_simplify(args):
    if args.PARAM and args.MACHINE:
        if args.debug:
            dlog.setLevel(logging.DEBUG)
        dlog.info("start simplifying")
        run_iter(args.PARAM, args.MACHINE)
        dlog.info("finished")

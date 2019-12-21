"""Simplify dataset (minimize the dataset size).

Init:
pick up init data from dataset randomly

Iter:
00: train models (same as generator)
01: calculate model deviations of the rest dataset, pick up data with proper model deviaiton 
02: fp (optional, if the original dataset do not have fp data, same as generator)
"""
import logging
import queue
import os
import json
import argparse
import pickle
import glob

import dpdata
import numpy as np

from dpgen import dlog
from dpgen import SHORT_CMD
from dpgen.util import sepline
from dpgen.remote.decide_machine import decide_train_machine
from dpgen.dispatcher.Dispatcher import Dispatcher, make_dispatcher
from dpgen.generator.run import make_train, run_train, post_train, run_fp, post_fp, fp_name, model_devi_name, train_name
# TODO: maybe the following functions can be moved to dpgen.util
from dpgen.generator.lib.utils import log_iter, make_iter_name, create_path, record_iter
from dpgen.remote.decide_machine import decide_train_machine, decide_fp_machine, decide_model_devi_machine
from dpgen.generator.lib.gaussian import make_gaussian_input


picked_data_name = "data.picked"
rest_data_name = "data.rest"
accurate_data_name = "data.accurate"
detail_file_name_prefix = "details"


def get_system_cls(jdata):
    if jdata.get("labeled", False):
        return dpdata.LabeledSystem
    return dpdata.System


def get_systems(path, jdata):
    system = get_system_cls(jdata)
    systems = dpdata.MultiSystems(
        *[system(os.path.join(path, s), fmt='deepmd/npy') for s in os.listdir(path)])
    return systems


def init_pick(iter_index, jdata, mdata):
    """pick up init data from dataset randomly"""
    pick_data = jdata['pick_data']
    init_pick_number = jdata['init_pick_number']
    # use MultiSystems with System
    # TODO: support System and LabeledSystem
    # TODO: support other format
    systems = get_systems(pick_data, jdata)
    # label the system
    labels = []
    for key, system in systems.systems.items():
        labels.extend([(key, j) for j in range(len(system))])

    # random pick
    iter_name = make_iter_name(iter_index)
    create_path(iter_name)
    work_path = os.path.join(iter_name, model_devi_name)
    create_path(work_path)
    idx = np.arange(len(labels))
    np.random.shuffle(idx)
    pick_idx = idx[:init_pick_number]
    rest_idx = idx[init_pick_number:]

    # dump the init data
    picked_systems = dpdata.MultiSystems()
    for j in pick_idx:
        sys_name, sys_id = labels[j]
        picked_systems.append(systems[sys_name][sys_id])
    sys_data_path = os.path.join(work_path, picked_data_name)

    picked_systems.to_deepmd_raw(sys_data_path)
    picked_systems.to_deepmd_npy(sys_data_path, set_size=init_pick_number)

    # dump the rest data
    rest_systems = dpdata.MultiSystems()
    for j in rest_idx:
        sys_name, sys_id = labels[j]
        rest_systems.append(systems[sys_name][sys_id])
    sys_data_path = os.path.join(work_path, rest_data_name)
    rest_systems.to_deepmd_raw(sys_data_path)
    rest_systems.to_deepmd_npy(sys_data_path, set_size=rest_idx.size)


def make_model_devi(iter_index, jdata, mdata):
    """calculate the model deviation of the rest idx"""
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
    for jj, subsystem in enumerate(os.listdir(rest_data_path)):
        task_name = "task.%03d.%06d" % (0, jj)
        task_path = os.path.join(work_path, task_name)
        create_path(task_path)
        os.symlink(os.path.abspath(os.path.join(rest_data_path, subsystem)),
                   os.path.abspath(os.path.join(task_path, rest_data_name)))
    return True


def run_model_devi(iter_index, jdata, mdata, dispatcher):
    """submit dp test tasks"""
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    # generate command
    commands = []
    tasks = glob.glob(os.path.join(work_path, "task.*"))
    run_tasks = [os.path.basename(ii) for ii in tasks]
    # get models
    models = glob.glob(os.path.join(work_path, "graph*pb"))
    model_names = [os.path.basename(ii) for ii in models]
    task_model_list = []
    for ii in model_names:
        task_model_list.append(os.path.join('..', ii))
    # get max data size
    data_size = max([len(dpdata.System(os.path.join(
        task, rest_data_name), fmt="deepmd/npy")) for task in tasks])
    # models
    commands = []
    detail_file_names = []
    for ii, mm in enumerate(task_model_list):
        detail_file_name = "{prefix}.{ii}".format(
            prefix=detail_file_name_prefix,
            ii=ii,
        )
        # TODO: support 0.x?
        command = "{python} -m deepmd test -m {model} -s {system} -n {numb_test} -d {detail_file}".format(
            python=mdata['python_test_path'],
            model=mm,
            system=rest_data_name,
            numb_test=data_size,
            detail_file=detail_file_name,
        )
        commands.append(command)
        detail_file_names.append(detail_file_name)
    # submit
    try:
        model_devi_group_size = mdata['model_devi_group_size']
    except:
        model_devi_group_size = 1

    forward_files = [rest_data_name]
    backward_files = sum([[pf+".e.out", pf+".f.out", pf+".v.out"] for pf in detail_file_names], [])

    dispatcher.run_jobs(mdata['model_devi_resources'],
                        commands,
                        work_path,
                        run_tasks,
                        model_devi_group_size,
                        model_names,
                        forward_files,
                        backward_files,
                        outlog='model_devi.log',
                        errlog='model_devi.log')


def post_model_devi(iter_index, jdata, mdata):
    """calculate the model deviation"""
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    tasks = glob.glob(os.path.join(work_path, "task.*"))

    e_trust_lo = jdata['e_trust_lo']
    e_trust_hi = jdata['e_trust_hi']
    f_trust_lo = jdata['f_trust_lo']
    f_trust_hi = jdata['f_trust_hi']

    sys_accurate = dpdata.MultiSystems()
    sys_candinate = dpdata.MultiSystems()
    sys_failed = dpdata.MultiSystems()

    for task in tasks:
        # e.out
        details_e = glob.glob(os.path.join(task, "{}.*.e.out".format(detail_file_name_prefix)))
        e_all = np.array([np.loadtxt(detail_e, ndmin=2)[:, 1] for detail_e in details_e])
        e_std = np.std(e_all, axis=0)
        n_frame = e_std.size
        
        # f.out
        details_f = glob.glob(os.path.join(task, "{}.*.f.out".format(detail_file_name_prefix)))
        f_all = np.array([np.loadtxt(detail_f, ndmin=2)[:, 3:6].reshape((n_frame, -1, 3)) for detail_f in details_f])
        # (n_model, n_frame, n_atom, 3)
        f_std = np.std(f_all, axis=0)
        # (n_frame, n_atom, 3)
        f_std = np.linalg.norm(f_std, axis=2)
        # (n_frame, n_atom)
        f_std = np.max(f_std, axis=1)
        # (n_frame,)

        system_cls = get_system_cls(jdata)
        for subsys, e_devi, f_devi in zip(system_cls(os.path.join(task, rest_data_name), fmt='deepmd/npy'), e_std, f_std):
            if (e_devi < e_trust_hi and e_devi >= e_trust_lo) or (f_devi < f_trust_hi and f_devi >= f_trust_lo) :
                sys_candinate.append(subsys)
            elif (e_devi >= e_trust_hi ) or (f_devi >= f_trust_hi ):
                sys_failed.append(subsys)
            elif (e_devi < e_trust_lo and f_devi < f_trust_lo ):
                sys_accurate.append(subsys)
    counter = {"candidate": sys_candinate.get_nframes(), "accurate": sys_accurate.get_nframes(), "failed": sys_failed.get_nframes()}
    fp_sum = sum(counter.values())
    for cc_key, cc_value in counter.items():
        dlog.info("{0:9s} : {1:6d} in {2:6d} {3:6.2f} %".format(cc_key, cc_value, fp_sum, cc_value/fp_sum*100))
    
    # label the candidate system
    labels = []
    for key, system in sys_candinate.systems.items():
        labels.extend([(key, j) for j in range(len(system))])
    # candinate: pick up randomly
    iter_pick_number = jdata['iter_pick_number']
    idx = np.arange(counter['candidate'])
    np.random.shuffle(idx)
    pick_idx = idx[:iter_pick_number]
    rest_idx = idx[iter_pick_number:]

    # dump the picked candinate data
    picked_systems = dpdata.MultiSystems()
    for j in pick_idx:
        sys_name, sys_id = labels[j]
        picked_systems.append(sys_candinate[sys_name][sys_id])
    sys_data_path = os.path.join(work_path, picked_data_name)

    picked_systems.to_deepmd_raw(sys_data_path)
    picked_systems.to_deepmd_npy(sys_data_path, set_size=iter_pick_number)

    # dump the rest data (not picked candinate data and failed data)
    rest_systems = dpdata.MultiSystems()
    for j in rest_idx:
        sys_name, sys_id = labels[j]
        rest_systems.append(sys_candinate[sys_name][sys_id])
    rest_systems += sys_failed
    sys_data_path = os.path.join(work_path, rest_data_name)
    rest_systems.to_deepmd_raw(sys_data_path)
    rest_systems.to_deepmd_npy(sys_data_path, set_size=rest_idx.size)

    # dump the accurate data -- to another directory
    sys_data_path = os.path.join(work_path, accurate_data_name)
    sys_accurate.to_deepmd_raw(sys_data_path)
    sys_accurate.to_deepmd_npy(sys_data_path, set_size=rest_idx.size)


def make_fp(iter_index, jdata, mdata):
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    picked_data_path = os.path.join(iter_name, model_devi_name, picked_data_name)
    if jdata.get("labeled", False):
        dlog.info("already labeled, skip make_fp and link data directly")
        os.symlink(os.path.abspath(picked_data_path), os.path.abspath(
            os.path.join(work_path, "task.%03d" % 0)))
        os.symlink(os.path.abspath(picked_data_path), os.path.abspath(
            os.path.join(work_path, "data.%03d" % 0)))
        return
    systems = get_systems(picked_data_path, jdata)
    fp_style = jdata['fp_style']
    if 'user_fp_params' in jdata.keys() :
        fp_params = jdata['user_fp_params']
    else:
        fp_params = jdata['fp_params']
    jj = 0
    for system in systems:
        for subsys in system:
            sys_data = subsys.data
            task_name = "task.%03d.%06d" % (0, jj)
            task_path = os.path.join(work_path, task_name)
            create_path(task_path)
            if fp_style == "gaussian" :
                ret = make_gaussian_input(sys_data, fp_params)
                with open(os.path.join(task_path, 'input'), 'w') as fp:
                    fp.write(ret)
            else :
                # TODO: support other formats
                raise RuntimeError ("unsupported fp style")
            jj += 1


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
    except:
        with open(param_file, 'r') as fp:
            jdata = json.load(fp)
        with open(machine_file, 'r') as fp:
            mdata = json.load(fp)

    if jdata.get('pretty_print', False):
        fparam = SHORT_CMD+'_' + \
            param_file.split('.')[0]+'.'+jdata.get('pretty_format', 'json')
        dumpfn(jdata, fparam, indent=4)
        fmachine = SHORT_CMD+'_' + \
            machine_file.split('.')[0]+'.'+jdata.get('pretty_format', 'json')
        dumpfn(mdata, fmachine, indent=4)

    if mdata.get('handlers', None):
        if mdata['handlers'].get('smtp', None):
            que = queue.Queue(-1)
            queue_handler = logging.handlers.QueueHandler(que)
            smtp_handler = logging.handlers.SMTPHandler(
                **mdata['handlers']['smtp'])
            listener = logging.handlers.QueueListener(que, smtp_handler)
            dlog.addHandler(queue_handler)
            listener.start()

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
                    init_pick(ii, jdata, mdata)
                dlog.info("first iter, skip step 1-5")
            elif jj == 0:
                log_iter("make_train", ii, jj)
                make_train(ii, jdata, mdata)
            elif jj == 1:
                log_iter("run_train", ii, jj)
                mdata = decide_train_machine(mdata)
                disp = make_dispatcher(mdata['train_machine'])
                run_train(ii, jdata, mdata, disp)
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
                mdata = decide_model_devi_machine(mdata)
                disp = make_dispatcher(mdata['model_devi_machine'])
                run_model_devi(ii, jdata, mdata, disp)
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
                    mdata = decide_fp_machine(mdata)
                    disp = make_dispatcher(mdata['fp_machine'])
                    run_fp(ii, jdata, mdata, disp)
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

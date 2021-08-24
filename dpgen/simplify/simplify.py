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
import fnmatch
import dpdata
import numpy as np

from dpgen import dlog
from dpgen import SHORT_CMD
from dpgen.util import sepline
from dpgen.dispatcher.Dispatcher import Dispatcher, make_dispatcher
from dpgen.generator.run import make_train, run_train, post_train, run_fp, post_fp, fp_name, model_devi_name, train_name, train_task_fmt, sys_link_fp_vasp_pp, make_fp_vasp_incar, make_fp_vasp_kp, make_fp_vasp_cp_cvasp, data_system_fmt, model_devi_task_fmt, fp_task_fmt
# TODO: maybe the following functions can be moved to dpgen.util
from dpgen.generator.lib.utils import log_iter, make_iter_name, create_path, record_iter
from dpgen.generator.lib.gaussian import make_gaussian_input


picked_data_name = "data.picked"
rest_data_name = "data.rest"
accurate_data_name = "data.accurate"
detail_file_name_prefix = "details"
sys_name_fmt = 'sys.' + data_system_fmt
sys_name_pattern = 'sys.[0-9]*[0-9]'

def expand_sys_str(root_dir):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir, followlinks=True):
        for filename in fnmatch.filter(filenames, 'type.raw'):
            matches.append(root)
    matches.sort()
    dirnames = [os.path.basename(ii) for ii in matches]
    if (len(list(set(dirnames))) != len(matches)) :
        raise RuntimeError('duplicated system name: it is highly recommend to place all systems in the same level of directory and has different names')
    return matches


def get_system_cls(jdata):
    if jdata.get("labeled", False):
        return dpdata.LabeledSystem
    return dpdata.System


def get_multi_system(path, jdata):
    system = get_system_cls(jdata)
    systems = dpdata.MultiSystems(
        *[system(os.path.join(path, s), fmt='deepmd/npy') for s in os.listdir(path)])
    return systems


def get_systems(path, jdata):
    system_cls = get_system_cls(jdata)
    system_paths = expand_sys_str(path)    
    systems = {}
    for ii in system_paths:
        systems[os.path.basename(ii)] = system_cls(ii, fmt='deepmd/npy')
    return systems


def get_system_idx(path):
    system_paths = expand_sys_str(path)    
    sys_idx_map = {}
    for idx,ii in enumerate(system_paths):
        sys_idx_map[os.path.basename(ii)] = idx
    return sys_idx_map


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
    use_clusters = jdata.get('use_clusters', False)
    # use MultiSystems with System
    # TODO: support System and LabeledSystem
    # TODO: support other format
    if use_clusters:
        systems = get_multi_system(pick_data, jdata)
    else:
        systems = get_systems(pick_data, jdata)
    # label the system
    labels = []
    if use_clusters:
        items = systems.systems.items()
    else:
        items = systems.items()
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


def _add_system(systems, key, system):
    if key in systems.keys():
        systems[key].append(system)
    else:
        systems[key] = system
    return systems


def _init_dump_selected_frames(systems, labels, selc_idx, sys_data_path, jdata):
    pick_data = jdata['pick_data']
    use_clusters = jdata.get('use_clusters', False)
    if use_clusters:
        selc_systems = dpdata.MultiSystems()
        for j in selc_idx:
            sys_name, sys_id = labels[j]
            selc_systems.append(systems[sys_name][sys_id])
        selc_systems.to_deepmd_raw(sys_data_path)
        selc_systems.to_deepmd_npy(sys_data_path, set_size=selc_idx.size)
    else:
        selc_systems = {}
        for j in selc_idx:
            sys_name, sys_id = labels[j]
            selc_systems = _add_system(selc_systems, sys_name, systems[sys_name][sys_id])
        sys_idx_map = get_system_idx(pick_data)
        for kk in selc_systems.keys():
            sub_path = os.path.join(sys_data_path, sys_name_fmt % sys_idx_map[kk])
            selc_systems[kk].to_deepmd_raw(sub_path)
            selc_systems[kk].to_deepmd_npy(sub_path, set_size=selc_idx.size)
        with open(os.path.join(sys_data_path, 'sys_idx_map.json'), 'w') as fp:
            json.dump(sys_idx_map, fp, indent=4)

def _dump_system_dict(systems, path):
    for kk in systems:
        sub_path = os.path.join(path, sys_name_fmt % (int(kk)))
        systems[kk].to_deepmd_raw(sub_path)
        systems[kk].to_deepmd_npy(sub_path, set_size=systems[kk].get_nframes())


def make_model_devi(iter_index, jdata, mdata):
    """calculate the model deviation of the rest idx"""
    pick_data = jdata['pick_data']
    use_clusters = jdata.get('use_clusters', False)
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
    if use_clusters:
        for jj, subsystem in enumerate(os.listdir(rest_data_path)):
            task_name = "task." + model_devi_task_fmt % (0, jj)
            task_path = os.path.join(work_path, task_name)
            create_path(task_path)
            os.symlink(os.path.abspath(os.path.join(rest_data_path, subsystem)),
                       os.path.abspath(os.path.join(task_path, rest_data_name)))
    else:
        rest_data_path = os.path.abspath(rest_data_path)
        sys_path = glob.glob(os.path.join(rest_data_path, sys_name_pattern))
        cwd = os.getcwd()
        for ii in sys_path:
            task_name = "task." + model_devi_task_fmt % (int(os.path.basename(ii).split('.')[1]), 0)
            task_path = os.path.join(work_path, task_name)
            create_path(task_path)            
            os.chdir(task_path)
            os.symlink(os.path.relpath(ii), rest_data_name)
            os.chdir(cwd)
        os.chdir(cwd)
    return True


def run_model_devi(iter_index, jdata, mdata):
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

    dispatcher = make_dispatcher(mdata['model_devi_machine'], mdata['model_devi_resources'], work_path, run_tasks, model_devi_group_size)
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
    use_clusters = jdata.get('use_clusters', False)
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    tasks = glob.glob(os.path.join(work_path, "task.*"))
    tasks.sort()

    e_trust_lo = jdata['e_trust_lo']
    e_trust_hi = jdata['e_trust_hi']
    f_trust_lo = jdata['f_trust_lo']
    f_trust_hi = jdata['f_trust_hi']

    if use_clusters:
        sys_accurate = dpdata.MultiSystems()
        sys_candinate = dpdata.MultiSystems()
        sys_failed = dpdata.MultiSystems()
    else:
        sys_accurate = {}
        sys_candinate = {}
        sys_failed = {}
        all_names = set()

    for task in tasks:
        if not use_clusters:
            sys_name = os.path.basename(task).split('.')[1]
            all_names.add(sys_name)
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
                if use_clusters:
                    sys_candinate.append(subsys)
                else:
                    sys_candinate = _add_system(sys_candinate, sys_name, subsys)
            elif (e_devi >= e_trust_hi ) or (f_devi >= f_trust_hi ):
                if use_clusters:
                    sys_failed.append(subsys)
                else:
                    sys_failed = _add_system(sys_failed, sys_name, subsys)
            elif (e_devi < e_trust_lo and f_devi < f_trust_lo ):
                if use_clusters:
                    sys_accurate.append(subsys)
                else:
                    sys_accurate = _add_system(sys_accurate, sys_name, subsys)
            else:
                raise RuntimeError('reach a place that should NOT be reached...')
    if use_clusters:
        counter = {"candidate": sys_candinate.get_nframes(), "accurate": sys_accurate.get_nframes(), "failed": sys_failed.get_nframes()}
        fp_sum = sum(counter.values())
        for cc_key, cc_value in counter.items():
            dlog.info("{0:9s} : {1:6d} in {2:6d} {3:6.2f} %".format(cc_key, cc_value, fp_sum, cc_value/fp_sum*100))
    else:
        all_names = list(all_names)
        all_names.sort()
        counter = {"candidate": 0, "accurate": 0, "failed": 0}
        for kk in all_names:
            sys_counter = {"candidate": 0, "accurate": 0, "failed": 0}
            if kk in sys_candinate.keys():
                sys_counter['candidate'] += sys_candinate[kk].get_nframes()
            if kk in sys_accurate.keys():
                sys_counter['accurate'] += sys_accurate[kk].get_nframes()
            if kk in sys_failed.keys():
                sys_counter['failed'] += sys_failed[kk].get_nframes()
            fp_sum = sum(sys_counter.values())
            for cc_key, cc_value in sys_counter.items():
                if fp_sum != 0:
                    dlog.info("sys{0:s} {1:9s} : {2:6d} in {3:6d} {4:6.2f} %".format(kk, cc_key, cc_value, fp_sum, cc_value/fp_sum*100))
                else:
                    dlog.info("sys{0:s} {1:9s} : {2:6d} in {3:6d} {4:6.2f} %".format(kk, cc_key, cc_value, fp_sum, 0*100))
            for ii in ['candidate', 'accurate', 'failed']:
                counter[ii] += sys_counter[ii]
    
    if counter['candidate'] == 0 and counter['failed'] > 0:
        raise RuntimeError('no candidate but still have failed cases, stop. You may want to refine the training or to increase the trust level hi')

    # label the candidate system
    labels = []
    if use_clusters:
        items = sys_candinate.systems.items()
    else:
        items = sys_candinate.items()
    for key, system in items:
        labels.extend([(key, j) for j in range(len(system))])
    # candinate: pick up randomly
    iter_pick_number = jdata['iter_pick_number']
    idx = np.arange(counter['candidate'])
    assert(len(idx) == len(labels))
    np.random.shuffle(idx)
    pick_idx = idx[:iter_pick_number]
    rest_idx = idx[iter_pick_number:]
    dlog.info("total candidate {0:6d}   picked {1:6d} ({2:6.2f} %) rest {3:6d} ({4:6.2f} % )".format\
              (counter['candidate'], len(pick_idx), float(len(pick_idx))/counter['candidate']*100., len(rest_idx), float(len(rest_idx))/counter['candidate']*100.))

    # dump the picked candinate data
    if use_clusters:
        picked_systems = dpdata.MultiSystems()
        for j in pick_idx:
            sys_name, sys_id = labels[j]
            picked_systems.append(sys_candinate[sys_name][sys_id])
        sys_data_path = os.path.join(work_path, picked_data_name)
        picked_systems.to_deepmd_raw(sys_data_path)
        picked_systems.to_deepmd_npy(sys_data_path, set_size=iter_pick_number)
    else:
        selc_systems = {}
        for j in pick_idx:
            sys_name, sys_id = labels[j]
            selc_systems = _add_system(selc_systems, sys_name, sys_candinate[sys_name][sys_id])
        sys_data_path = os.path.join(work_path, picked_data_name)
        _dump_system_dict(selc_systems, sys_data_path)

    # dump the rest data (not picked candinate data and failed data)
    if use_clusters:
        rest_systems = dpdata.MultiSystems()
        for j in rest_idx:
            sys_name, sys_id = labels[j]
            rest_systems.append(sys_candinate[sys_name][sys_id])
        rest_systems += sys_failed
        sys_data_path = os.path.join(work_path, rest_data_name)
        rest_systems.to_deepmd_raw(sys_data_path)
        rest_systems.to_deepmd_npy(sys_data_path, set_size=rest_idx.size)
    else:
        selc_systems = {}
        for j in rest_idx:
            sys_name, sys_id = labels[j]
            selc_systems = _add_system(selc_systems, sys_name, sys_candinate[sys_name][sys_id])
        for kk in sys_failed.keys():
            selc_systems = _add_system(selc_systems, kk, sys_failed[kk])        
        sys_data_path = os.path.join(work_path, rest_data_name)
        _dump_system_dict(selc_systems, sys_data_path)

    # dump the accurate data -- to another directory
    if use_clusters:
        sys_data_path = os.path.join(work_path, accurate_data_name)
        sys_accurate.to_deepmd_raw(sys_data_path)
        sys_accurate.to_deepmd_npy(sys_data_path, set_size=sys_accurate.get_nframes())
    else:
        sys_data_path = os.path.join(work_path, accurate_data_name)
        _dump_system_dict(sys_accurate, sys_data_path)


def make_fp_labeled(iter_index, jdata):    
    dlog.info("already labeled, skip make_fp and link data directly")
    pick_data = jdata['pick_data']
    use_clusters = jdata.get('use_clusters', False)
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    picked_data_path = os.path.join(iter_name, model_devi_name, picked_data_name)
    if use_clusters:
        os.symlink(os.path.abspath(picked_data_path), os.path.abspath(
            os.path.join(work_path, "task." + data_system_fmt % 0)))
        os.symlink(os.path.abspath(picked_data_path), os.path.abspath(
            os.path.join(work_path, "data." + data_system_fmt % 0)))
    else:
        picked_data_path = os.path.abspath(picked_data_path)
        sys_path = glob.glob(os.path.join(picked_data_path, sys_name_pattern))
        cwd = os.getcwd()
        os.chdir(work_path)
        for ii in sys_path:
            sys_idx = os.path.basename(ii).split('.')[1]
            data_dir = 'data.' + data_system_fmt % int(sys_idx)
            task_dir = 'task.' + data_system_fmt % int(sys_idx)
            os.symlink(os.path.relpath(ii), data_dir)
            os.symlink(os.path.relpath(ii), task_dir)
        os.chdir(cwd)


def make_fp_configs(iter_index, jdata):
    pick_data = jdata['pick_data']
    use_clusters = jdata.get('use_clusters', False)
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    picked_data_path = os.path.join(iter_name, model_devi_name, picked_data_name)
    if use_clusters:
        systems = get_multi_system(picked_data_path, jdata)
        jj = 0
        for system in systems:
            for subsys in system:
                task_name = "task." + fp_task_fmt % (0, jj)
                task_path = os.path.join(work_path, task_name)
                create_path(task_path)
                subsys.to('vasp/poscar', os.path.join(task_path, 'POSCAR'))
                jj += 1
    else:
        picked_data_path = os.path.abspath(picked_data_path)
        sys_path = glob.glob(os.path.join(picked_data_path, sys_name_pattern))
        for ii in sys_path:
            tmp_sys = dpdata.System(ii, fmt = 'deepmd/npy')
            sys_idx = os.path.basename(ii).split('.')[1]
            jj = 0
            for ss in tmp_sys:
                task_name = "task." + fp_task_fmt % (int(sys_idx), jj)
                task_path = os.path.join(work_path, task_name)
                create_path(task_path)
                ss.to('vasp/poscar', os.path.join(task_path, 'POSCAR'))
                job = {}
                with open(os.path.join(task_path, 'job.json'), 'w') as fp:
                    json.dump(job, fp, indent=4)
                jj += 1


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

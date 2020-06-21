import os
from dpgen.auto_test.EOS import EOS
from dpgen.auto_test.Elastic import Elastic
from dpgen.auto_test.Vacancy import Vacancy
from dpgen.auto_test.Interstitial import Interstitial
from dpgen.auto_test.Surface import Surface
from dpgen.auto_test.common_task import make_task,make_task_trans_files

import dpgen.auto_test.lib.crys as crys
import glob, warnings, json
from dpgen.remote.decide_machine import decide_fp_machine, decide_model_devi_machine
import dpgen.auto_test.lib.util as util
from dpgen.dispatcher.Dispatcher import make_dispatcher

lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']

def make_property_instance(paramters):
    """
    Make an instance of Property
    """
    prop_type = paramters['type']
    if prop_type == 'eos':
        return EOS(paramters)
    elif prop_type == 'elastic':
        return Elastic(paramters)
    elif prop_type == 'vacancy':
        return Vacancy(paramters)
    elif prop_type == 'interstitial':
        return Interstitial(paramters)
    elif prop_type == 'surface':
        return Surface(paramters)
    else:
        raise RuntimeError(f'unknown property type {prop_type}')


def make_property(confs,
                  inter_param,
                  property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    conf_dirs = glob.glob(confs)
    conf_dirs.sort()
    for ii in conf_dirs:
        for jj in property_list:
            if 'init_from_suffix' and 'output_suffix' in jj:
                do_refine = True
                suffix = jj['output_suffix']
            elif 'reprod-opt' in jj and jj['reprod-opt']:
                do_refine = False
                suffix = 'reprod'
            else:
                do_refine = False
                suffix = '00'
            # generate working directory like mp-xxx/eos_00 if jj['type'] == 'eos'
            # handel the exception that the working directory exists
            # ...

            # determine the suffix: from scratch or refine
            # ...

            property_type = jj['type']
            path_to_equi = os.path.join(ii, 'relaxation')
            path_to_work = os.path.join(ii, property_type + '_' + suffix)

            if os.path.exists(path_to_work):
                warnings.warn('%s already exists' % path_to_work)
            else:
                os.makedirs(path_to_work)

            prop = make_property_instance(jj)
            task_list = prop.make_confs(path_to_work, path_to_equi, do_refine)

            for kk in task_list:
                poscar = os.path.join(kk, 'POSCAR')
                inter = make_task(inter_param, poscar)
                inter.make_potential_files(kk)
                inter.make_input_file(kk, prop.task_type, prop.task.pararm)


def run_property(confs,
                 inter_param,
                 property_list,
                 mdata):
    # find all POSCARs and their name like mp-xxx
    # ...
    conf_dirs = glob.glob(confs)
    conf_dirs.sort()
    task_list = []
    work_path_list = []
    for ii in conf_dirs:
        for jj in property_list:
            # determine the suffix: from scratch or refine
            # ...
            if 'init_from_suffix' and 'output_suffix' in jj:
                suffix = jj['output_suffix']
            elif 'reprod-opt' in jj and jj['reprod-opt']:
                suffix = 'reprod'
            else:
                suffix = 0

            property_type = jj['type']
            path_to_work = os.path.join(ii, property_type + '_' + suffix)

            work_path_list.append(path_to_work)
            tmp_task_list = glob.glob(os.path.join(path_to_work, 'task.[0-9]*[0-9]'))
            tmp_task_list.sort()
            task_list.append(tmp_task_list)

    # dispatch the tasks
    forward_files, forward_common_files, backward_files = make_task_trans_files(inter_param)
    #    backward_files += logs
    # ...
    inter_type = inter_param['type']
    # vasp
    if inter_type == "vasp":
        mdata = decide_fp_machine(mdata)
    elif inter_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)
    else:
        raise RuntimeError("unknown task %s, something wrong" % inter_type)

    for ii in range(len(work_path_list)):
        work_path = work_path_list[ii]
        all_task = task_list[ii]
        run_tasks = util.collect_task(all_task, inter_type)
        if len(run_tasks) == 0:
            return
        else:
            run_tasks = [os.path.basename(ii) for ii in all_task]
            machine, resources, command, group_size = util.get_machine_info(mdata, inter_type)
            disp = make_dispatcher(machine, resources, work_path, run_tasks, group_size)
            disp.run_jobs(resources,
                          command,
                          work_path,
                          run_tasks,
                          group_size,
                          forward_common_files,
                          forward_files,
                          backward_files,
                          outlog=inter_type + '.out',
                          errlog=inter_type + '.err')


def post_property(confs,
                  #                  inter_param,
                  property_list):
    # find all POSCARs and their name like mp-xxx
    # ...
    #    task_list = []
    conf_dirs = glob.glob(confs)
    conf_dirs.sort()
    for ii in conf_dirs:
        for jj in property_list:
            # determine the suffix: from scratch or refine
            # ...
            if 'init_from_suffix' and 'output_suffix' in jj:
                suffix = jj['output_suffix']
            elif 'reprod-opt' in jj and jj['reprod-opt']:
                suffix = 'reprod'
            else:
                suffix = 0
            property_type = jj['type']
            path_to_work = os.path.join(ii, property_type + '_' + suffix)
            prop = make_property_instance(jj)
            prop.compute('result.json', 'result.out', path_to_work)

import os
from VASP import VASP
from DEEPMD_LMP import DEEPMD_LMP
from MEAM_LMP import MEAM_LMP
from EAM_FS_LMP import EAM_FS_LMP
from EAM_ALLOY_LMP import EAM_ALLOY_LMP
from EOS import EOS
from Elastic import Elastic
from Vacancy import Vacancy
from Interstitial import Interstitial
from Surface import Surface
import dpgen.auto_test.lib.crys as crys
import glob, warnings, json
from dpgen.remote.decide_machine import decide_fp_machine, decide_model_devi_machine
import dpgen.auto_test.lib.util as util
from dpgen.dispatcher.Dispatcher import make_dispatcher

lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']


def make_task(inter_parameter,
              path_to_poscar):
    """
    Make an instance of Task
    """
    inter_type = inter_parameter['type']
    if inter_type == 'vasp':
        return VASP(inter_parameter, path_to_poscar)
    elif inter_type == 'deepmd':
        return DEEPMD_LMP(inter_parameter, path_to_poscar)
    elif inter_type == 'meam':
        return MEAM_LMP(inter_parameter, path_to_poscar)
    elif inter_type == 'eam_fs':
        return EAM_FS_LMP(inter_parameter, path_to_poscar)
    elif inter_type == 'eam_alloy':
        return EAM_ALLOY_LMP(inter_parameter, path_to_poscar)
    else:
        raise RuntimeError(f'unknown interaction {inter_type}')


def make_task_trans_files(inter_parameter):
    """
    Make the forward and backward file of an Task
    """
    inter_type = inter_parameter['type']
    if inter_type == 'vasp':
        return VASP.forward_files, VASP.forward_common_files, VASP.backward_files
    elif inter_type == 'deepmd':
        return DEEPMD_LMP.forward_files, DEEPMD_LMP.forward_common_files, DEEPMD_LMP.backward_files
    elif inter_type == 'meam':
        return MEAM_LMP.forward_files, MEAM_LMP.forward_common_files, MEAM_LMP.backward_files
    elif inter_type == 'eam_fs':
        return EAM_FS_LMP.forward_files, EAM_FS_LMP.forward_common_files, EAM_FS_LMP.backward_files
    elif inter_type == 'eam_alloy':
        return EAM_ALLOY_LMP.forward_files, EAM_ALLOY_LMP.forward_common_files, EAM_ALLOY_LMP.backward_files
    else:
        raise RuntimeError(f'unknown interaction {inter_type}')


def make_property(paramters):
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


def make_equi(confs,
              inter_param,
              relax_param):
    # find all POSCARs and their name like mp-xxx
    # ...
    ele_list = [key for key in inter_param['type_map'].keys()]
    conf_dirs = glob.glob(confs)
    conf_dirs.sort()

    # generate a list of task names like mp-xxx/relaxation
    # ...
    cwd = os.getcwd()
    # generate poscar for single element crystal
    if len(ele_list) == 1:
        for ii in conf_dirs:
            os.chdir(ii)
            crys_type = ii[3:]
            if crys_type == 'fcc':
                if not os.path.exists('POSCAR'):
                    crys.fcc(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'hcp':
                if not os.path.exists('POSCAR'):
                    crys.hcp(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'dhcp':
                if not os.path.exists('POSCAR'):
                    crys.dhcp(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'bcc':
                if not os.path.exists('POSCAR'):
                    crys.bcc(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'diamond':
                if not os.path.exists('POSCAR'):
                    crys.diamond(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'sc':
                if not os.path.exists('POSCAR'):
                    crys.sc(ele_list[0]).to('POSCAR', 'POSCAR')
            os.chdir(cwd)
    task_dirs = []
    # make task directories like mp-xxx/relaxation
    # if mp-xxx/exists then print a warning and exit.
    # ...
    for ii in conf_dirs:
        poscar = os.path.abspath(os.path.join(ii, 'POSCAR'))
        if not os.path.exists(poscar):
            raise FileNotFoundError('no configuration for autotest')
        relax_dirs = os.path.abspath(os.path.join(ii, 'relaxation'))
        if os.path.exists(relax_dirs):
            warnings.warn('%s already exists' % relax_dirs)
        else:
            os.makedirs(relax_dirs)
            task_dirs.append(relax_dirs)
            os.chdir(relax_dirs)
            # copy POSCARs to mp-xxx/relaxation
            # ...
            os.symlink(os.path.relpath(poscar), 'POSCAR')
            os.chdir(cwd)
    task_dirs.sort()
    # generate task files
    for ii in task_dirs:
        poscar = os.path.join(ii, 'POSCAR')
        inter = make_task(inter_param, poscar)
        inter.make_potential_files(ii)
        inter.make_input_file(ii, 'relaxation', relax_param)


def run_equi(confs,
             inter_param,
             mdata):
    # find all POSCARs and their name like mp-xxx
    # ...
    conf_dirs = glob.glob(confs)
    conf_dirs.sort()
    # generate a list of task names like mp-xxx/relaxation
    # ...
    work_path_list = []
    for ii in conf_dirs:
        work_path_list.append(os.path.join(ii, 'relaxation'))
    all_task = []
    for ii in work_path_list:
        all_task.append(os.path.join(ii, '.'))

    inter_type = inter_param['type']
    # vasp
    if inter_type == "vasp":
        mdata = decide_fp_machine(mdata)
    elif inter_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)
    else:
        raise RuntimeError("unknown task %s, something wrong" % task_type)

    # dispatch the tasks
    forward_files, forward_common_files, backward_files = make_task_trans_files(inter_param)
    #    backward_files += logs
    # ...
    run_tasks = util.collect_task(all_task, inter_type)
    if len(run_tasks) == 0:
        return
    else:
        run_tasks = [os.path.basename(ii) for ii in all_task]
        machine, resources, command, group_size = util.get_machine_info(mdata, inter_type)
        for ii in range(len(work_path_list)):
            work_path = work_path_list[ii]
            disp = make_dispatcher(machine, resources, work_path, run_tasks[ii], group_size)
            disp.run_jobs(resources,
                          command,
                          work_path,
                          run_tasks[ii],
                          group_size,
                          forward_common_files,
                          forward_files,
                          backward_files,
                          outlog=inter_type + '.out',
                          errlog=inter_type + '.err')


def post_equi(confs, inter_param):
    # find all POSCARs and their name like mp-xxx
    # ...
    conf_dirs = glob.glob(confs)
    conf_dirs.sort()
    task_dirs = []
    for ii in conf_dirs:
        task_dirs.append(os.path.abspath(os.path.join(ii, 'relaxation')))
    task_dirs.sort()

    # generate a list of task names like mp-xxx/relaxation
    # ...

    # dump the relaxation result.
    for ii in task_dirs:
        poscar = os.path.join(ii, 'POSCAR')
        inter = make_task(inter_param, poscar)
        res = inter.compute(ii)
        with open(os.path.join(ii, 'result.json'), 'w') as fp:
            json.dump(res, fp, indent=4)


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
            else:
                do_refine = False
                suffix = 0
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

            prop = make_property(jj)
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
    task_type = inter_param['type']
    # vasp
    if task_type == "vasp":
        mdata = decide_fp_machine(mdata)
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)
    else:
        raise RuntimeError("unknown task %s, something wrong" % task_type)

    for ii in range(len(work_path_list)):
        work_path = work_path_list[ii]
        all_task = task_list[ii]
        run_tasks = util.collect_task(all_task, task_type)
        if len(run_tasks) == 0:
            return
        else:
            run_tasks = [os.path.basename(ii) for ii in all_task]
            machine, resources, command, group_size = util.get_machine_info(mdata, task_type)
            disp = make_dispatcher(machine, resources, work_path, run_tasks, group_size)
            disp.run_jobs(resources,
                          command,
                          work_path,
                          run_tasks,
                          group_size,
                          forward_common_files,
                          forward_files,
                          backward_files,
                          outlog=task_type + '.out',
                          errlog=task_type + '.err')


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
            else:
                suffix = 0
            property_type = jj['type']
            path_to_work = os.path.join(ii, property_type + '_' + suffix)
            prop = make_property(jj)
            prop.compute('result.json', 'result.out', path_to_work)

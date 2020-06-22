import os
import dpgen.auto_test.lib.crys as crys
import glob, warnings, json
import dpgen.auto_test.lib.util as util
from dpgen import dlog
from dpgen.dispatcher.Dispatcher import make_dispatcher
from dpgen.auto_test.calculator import make_calculator
from dpgen.remote.decide_machine import decide_fp_machine, decide_model_devi_machine
from dpgen.auto_test.mpdb import get_structure
from monty.serialization import loadfn, dumpfn

lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']


def make_equi(confs,
              inter_param,
              relax_param):
    # find all POSCARs and their name like mp-xxx
    # ...
    dlog.debug('debug info make equi')
    if 'type_map' in inter_param:
        ele_list = [key for key in inter_param['type_map'].keys()]
    else:
        ele_list = [key for key in inter_param['potcars'].keys()]
    dlog.debug("ele_list %s" % ':'.join(ele_list))
    conf_dirs = []
    for conf in confs:
        conf_dirs.extend(glob.glob(conf))
    conf_dirs.sort()

    # generate a list of task names like mp-xxx/relaxation
    # ...
    cwd = os.getcwd()
    # generate poscar for single element crystal
    if len(ele_list) == 1:
        for ii in conf_dirs:
            os.chdir(ii)
            crys_type = ii.split('/')[-1]
            dlog.debug('crys_type: %s' % crys_type)
            dlog.debug('pwd: %s' % os.getcwd())
            if crys_type == 'std-fcc':
                if not os.path.exists('POSCAR'):
                    crys.fcc(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'std-hcp':
                if not os.path.exists('POSCAR'):
                    crys.hcp(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'std-dhcp':
                if not os.path.exists('POSCAR'):
                    crys.dhcp(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'std-bcc':
                if not os.path.exists('POSCAR'):
                    crys.bcc(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'std-diamond':
                if not os.path.exists('POSCAR'):
                    crys.diamond(ele_list[0]).to('POSCAR', 'POSCAR')
            elif crys_type == 'std-sc':
                if not os.path.exists('POSCAR'):
                    crys.sc(ele_list[0]).to('POSCAR', 'POSCAR')

            os.chdir(cwd)
    task_dirs = []
    # make task directories like mp-xxx/relaxation
    # if mp-xxx/exists then print a warning and exit.
    # ...
    for ii in conf_dirs:
        crys_type = ii.split('/')[-1]
        dlog.debug('crys_type: %s' % crys_type)

        if 'mp-' in crys_type and not os.path.exists(os.path.join(ii, 'POSCAR')):
            get_structure(crys_type).to('POSCAR', os.path.join(ii, 'POSCAR'))

        poscar = os.path.abspath(os.path.join(ii, 'POSCAR'))
        if not os.path.exists(poscar):
            raise FileNotFoundError('no configuration for autotest')
        relax_dirs = os.path.abspath(os.path.join(ii, 'relaxation'))
        if os.path.exists(relax_dirs):
            dlog.warning('%s already exists' % relax_dirs)
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
        dlog.debug('task_dir %s' % ii)
        inter = make_calculator(inter_param, poscar)
        inter.make_potential_files(ii)
        inter.make_input_file(ii, 'relaxation', relax_param)


def run_equi(confs,
             inter_param,
             mdata):
    # find all POSCARs and their name like mp-xxx
    # ...
    conf_dirs = []
    for conf in confs:
        conf_dirs.extend(glob.glob(conf))
    conf_dirs.sort()
    # generate a list of task names like mp-xxx/relaxation
    # ...
    work_path_list = []
    for ii in conf_dirs:
        work_path_list.append(os.path.abspath(os.path.join(ii, 'relaxation')))
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
        raise RuntimeError("unknown task %s, something wrong" % inter_type)

    # dispatch the tasks
    # POSCAR here is useless
    virutual_calculator = make_calculator(inter_param, "POSCAR")
    forward_files = virutual_calculator.forward_files()
    forward_common_files = virutual_calculator.forward_common_files()
    backward_files = virutual_calculator.backward_files()
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
                          outlog='outlog',
                          errlog='errlog')


def post_equi(confs, inter_param):
    # find all POSCARs and their name like mp-xxx
    # ...
    conf_dirs = []
    for conf in confs:
        conf_dirs.extend(glob.glob(conf))
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
        inter = make_calculator(inter_param, poscar)
        res = inter.compute(ii, inter_param)

        dumpfn(res, os.path.join(ii, 'result.json'), indent=4)

import numpy as np
import requests
import os,re
from dpgen import dlog
from dpgen.auto_test.lib import vasp
from dpgen.auto_test.lib import lammps
from dpgen.auto_test.lib.utils import cmd_append_log

lammps_task_type=['deepmd','meam','eam_fs','eam_alloy']   # 06/13 revised

def voigt_to_stress(inpt) :
    ret = np.zeros((3,3))
    ret[0][0] = inpt[0]
    ret[1][1] = inpt[1]
    ret[2][2] = inpt[2]
    ret[0][1] = inpt[3]
    ret[0][2] = inpt[4]
    ret[1][2] = inpt[5]
    ret[2][0] = ret[0][2]
    ret[1][0] = ret[0][1]
    ret[2][1] = ret[1][2]
    return ret

def insert_data(task,task_type,username,file_name):
    assert task in ['eos','elastic','surf']
    assert task_type in ['vasp','deepmd']
    url='http://115.27.161.2:5000/insert_test_data?username=%s&expr_type=%s&data_type=%s' % (username,task_type,task)
    res = requests.post(url, data=open(file_name).read())
    print('Successful upload!')


def make_work_path(jdata,task,reprod_opt,static,user):

    task_type=jdata['task_type']
    conf_dir=jdata['conf_dir']
    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', task, conf_path)

    if task_type=="vasp":
        if user:
            work_path=os.path.join(task_path, 'vasp-user_incar')
            assert(os.path.isdir(work_path))
            return work_path
        if static:
            if 'scf_incar' in jdata.keys():
                task_type=task_type+'-static-scf_incar'
            else:
                kspacing = jdata['vasp_params']['kspacing']
                task_type=task_type+'-static-k%.2f' % (kspacing)
        else:
            if 'relax_incar' in jdata.keys():
                task_type=task_type+'-relax_incar'
            else:
                kspacing = jdata['vasp_params']['kspacing']
                task_type=task_type+'-k%.2f' % (kspacing)
    elif task_type in lammps_task_type:
        if static:
            task_type=task_type+'-static'
        elif reprod_opt :
            if 'relax_incar' in jdata.keys():
                task_type=task_type+'-reprod-relax_incar'
            else:
                kspacing = jdata['vasp_params']['kspacing']
                task_type=task_type+'-reprod-k%.2f'% (kspacing)

    work_path=os.path.join(task_path, task_type)
    assert(os.path.isdir(work_path))
    return work_path


def get_machine_info(mdata,task_type):
    if task_type=="vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        command = vasp_exec
        command = cmd_append_log(command, "log")
    elif task_type in lammps_task_type:
        model_devi_exec = mdata['model_devi_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        command = model_devi_exec + " -i in.lammps"
        command = cmd_append_log(command, "model_devi.log")
    return machine, resources, command, group_size

def collect_task(all_task,task_type):

    if task_type == 'vasp':
        output_file ='OUTCAR'
        check_finished = vasp.check_finished
    elif task_type in lammps_task_type:
        output_file = 'log.lammps'
        check_finished = lammps.check_finished

    run_tasks_ = []
    for ii in all_task:
        fres = os.path.join(ii, output_file)
        if os.path.isfile(fres) :
            if not check_finished(fres):
                run_tasks_.append(ii)
        else :
            run_tasks_.append(ii)

    run_tasks = [os.path.basename(ii) for ii in run_tasks_]
    return run_tasks

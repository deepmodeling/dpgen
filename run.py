#!/usr/bin/evn python3

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
import numpy as np
import subprocess as sp
from lib.utils import make_iter_name
from lib.utils import create_path
from lib.utils import copy_file_list
from lib.utils import replace
from lib.utils import cmd_append_log
from lib.utils import log_iter
from lib.utils import record_iter
from lib.utils import log_task
from lib.lammps import cvt_lammps_conf
from lib.lammps import make_lammps_input
from lib.vasp import make_vasp_incar
import lib.MachineLocal as MachineLocal
import lib.MachineSlurm as MachineSlurm
import lib.MachinePBS as MachinePBS
from lib.machine_exec import exec_hosts
from lib.machine_exec import exec_hosts_batch

template_name = 'template'
train_name = '00.train'
train_task_fmt = '%03d'
train_tmpl_path = os.path.join(template_name, train_name)
data_system_fmt = '%03d'
model_devi_name = '01.model_devi'
model_devi_task_fmt = data_system_fmt + '.%06d'
model_devi_conf_fmt = data_system_fmt + '.%04d'
fp_name = '02.fp'
fp_task_fmt = data_system_fmt + '.%06d'

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

def check_empty_iter(iter_index) :
    fp_path = os.path.join(make_iter_name(iter_index), fp_name)
    fp_data_sys = glob.glob(os.path.join(fp_path, "data.*"))
    return (len(fp_data_sys) == 0)

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

def make_train (iter_index, 
               jdata) :    
    # load json param
    train_param = jdata['train_param']
    if iter_index > 0 :
        stop_batch = jdata['res_stop_batch']
        start_lr = jdata['res_start_lr']
        decay_steps = jdata['res_decay_steps']
        decay_rate = jdata['res_decay_rate']
    numb_models = jdata['numb_models']
    init_data_sys_ = jdata['init_data_sys']    
    
    if iter_index > 0 and check_empty_iter(iter_index-1) :
        log_task('prev data is empty, copy prev model')
        copy_model(numb_models, iter_index-1, iter_index)
        return
    else :
        iter_name = make_iter_name(iter_index)
        work_path = os.path.join(iter_name, train_name)
        copy_flag = os.path.join(work_path, 'copied')
        if os.path.isfile(copy_flag) :
            os.remove(copy_flag)

    init_data_sys = []
    for ii in init_data_sys_ :
        init_data_sys.append(os.path.abspath(ii))
    if iter_index > 0 :
        for ii in range(iter_index) :
            fp_path = os.path.join(make_iter_name(ii), fp_name)
            fp_data_sys = glob.glob(os.path.join(fp_path, "data.*"))
            for jj in fp_data_sys :
                init_data_sys.append(os.path.abspath(jj))
    for ii in init_data_sys :
        if not os.path.isdir(ii) :
            raise RuntimeError ("data sys %s does not exists" % ii)
    # establish work path
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    create_path(work_path)
    # establish tasks
    jinput = jdata['default_training_param']
    jinput['systems'] = init_data_sys    
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        create_path(task_path)
        jinput['seed'] = random.randrange(sys.maxsize)
        with open(os.path.join(task_path, train_param), 'w') as outfile:
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

def run_train (iter_index,
               jdata, 
               exec_machine) :    
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = jdata['deepmd_path']
    train_nthreads = jdata['train_nthreads']
    train_param = jdata['train_param']
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
    command = os.path.join(deepmd_path, 'bin/dp_train') + ' ' + train_param
    if iter_index > 0:
        command += ' --init-model old/model.ckpt '
    command = cmd_append_log (command, 'train.log')
    # train models
    exec_hosts_batch(exec_machine, command, train_nthreads, all_task, None, verbose = True)

def post_train (iter_index,
                jdata) :
    # load json param
    numb_models = jdata['numb_models']
    deepmd_path = jdata['deepmd_path']
    # paths
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, train_name)
    # check if is copied
    copy_flag = os.path.join(work_path, 'copied')
    if os.path.isfile(copy_flag) :
        log_task('copied model, do not post train')
        return
    all_task = []
    for ii in range(numb_models) :
        task_path = os.path.join(work_path, train_task_fmt % ii)
        all_task.append(task_path)
    command = os.path.join(deepmd_path, 'bin/dp_frz')
    command = cmd_append_log(command, 'freeze.log')
    # frz models
    exec_hosts(MachineLocal, command, 1, all_task, None)
    # symlink models
    for ii in range(numb_models) :
        task_file = os.path.join(train_task_fmt % ii, 'frozen_model.pb')
        ofile = os.path.join(work_path, 'graph.%03d.pb' % ii)
        if os.path.isfile(ofile) :
            os.remove(ofile)
        os.symlink(task_file, ofile)    

def make_model_devi (iter_index, 
                     jdata) :
    model_devi_jobs = jdata['model_devi_jobs']
    ensemble = model_devi_jobs['ensemble']
    nsteps = model_devi_jobs['nsteps']
    trj_freq = model_devi_jobs['trj_freq']
    job_names = get_job_names (model_devi_jobs)
    assert (iter_index < len(job_names)) 
    cur_job_name = job_names[iter_index]    
    cur_job = model_devi_jobs[cur_job_name]
    if 'ensemble' in cur_job.keys() :
        ensemble = cur_job['ensemble']
    if 'nsteps' in cur_job.keys() :
        nsteps = cur_job['nsteps']
    if 'trj_freq' in cur_job.keys() :
        trj_freq = cur_job['trj_freq']

    conf_systems_glob = cur_job['systems']
    conf_systems = []
    for ss in conf_systems_glob :
        cur_systems = []
        for ii in ss :
            cur_systems += glob.glob(ii)
        cur_systems.sort()
        conf_systems.append (cur_systems)
    mass_map = cur_job['mass_map']
    temps = cur_job['temps']
    press = cur_job['press']

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

    all_task = []
    task_param = []
    conf_path = os.path.join(work_path, 'confs')
    create_path(conf_path)
    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        for cc in ss :            
            conf_name = make_model_devi_conf_name(sys_counter, conf_counter)
            poscar_name = conf_name + '.poscar'
            lmp_name = conf_name + '.lmp'
            os.symlink(cc, os.path.join(conf_path, poscar_name))
            all_task.append(conf_path)
            task_param.append(' ' + poscar_name + ' ' + lmp_name)
            conf_counter += 1
        sys_counter += 1
    exec_hosts(MachineLocal, os.path.join(os.getcwd(), "lib/ovito_file_convert.py"), 1, all_task, task_param)

    sys_counter = 0
    for ss in conf_systems:
        conf_counter = 0
        task_counter = 0
        for cc in ss :            
            for tt in temps:
                for pp in press:
                    task_name = make_model_devi_task_name(sys_counter, task_counter)
                    conf_name = make_model_devi_conf_name(sys_counter, conf_counter) + '.lmp'
                    task_path = os.path.join(work_path, task_name)
                    # print(task_path)
                    create_path(task_path)
                    create_path(os.path.join(task_path, 'traj'))
                    loc_conf_name = 'conf.lmp'
                    os.symlink(os.path.join(os.path.join('..','confs'), conf_name), 
                               os.path.join(task_path, loc_conf_name) )
                    file_c = make_lammps_input(ensemble,
                                               loc_conf_name,
                                               task_model_list,
                                               nsteps,
                                               trj_freq,
                                               mass_map,
                                               tt,
                                               pres = pp)
                    with open(os.path.join(task_path, 'input.lammps'), 'w') as fp :
                        fp.write(file_c)
                    # cvt_lammps_conf(cc, 'conf.lmp')
                    task_counter += 1
            conf_counter += 1
        sys_counter += 1

def run_model_devi (iter_index, 
                    jdata, 
                    exec_machine) :
    model_devi_np = jdata['model_devi_np']
    lmp_path = jdata['lmp_path']
    lmp_exec = os.path.join(lmp_path, 'lmp_mpi')

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = lmp_exec + " -i input.lammps"
    command = cmd_append_log(command, "model_devi.log")
    
    exec_hosts_batch(exec_machine, command, model_devi_np, all_task, verbose = True, mpi = True)

def post_model_devi (iter_index, 
                     jdata) :
    model_devi_trust = jdata['model_devi_trust']
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert(os.path.isdir(work_path))
    all_task = glob.glob(os.path.join(work_path, "task.*"))
    all_task.sort()
    command = "awk '{if ($2 > %e) {print $1,$2}}' model_devi.out > sel.out" % model_devi_trust

    exec_hosts (MachineLocal, command, 1, all_task)

def _make_fp_vasp_inner (modd_path,
                         work_path,
                         system_index,
                         fp_task_max,
                         fp_link_files,
                         fp_params):
    """
    modd_path           string          path of model devi
    work_path           string          path of fp
    system_index        [string]        index of systems
    fp_task_max         int             max number of tasks
    fp_link_files       [string]        linked files for fp, POTCAR for example
    fp_params           map             parameters for fp
    """
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    fp_tasks = []
    count_total = 0
    for ss in system_index :
        modd_system_glob = os.path.join(modd_path, 'task.' + ss + '.*')
        modd_system_task = glob.glob(modd_system_glob)
        modd_system_task.sort()
        cc = 0
        for tt in modd_system_task :
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sel_conf = np.loadtxt(os.path.join(tt, 'sel.out'))
                sel_conf = np.reshape(sel_conf, [-1,2])
                if sel_conf.shape[0] == 0:
                    continue
                sel_conf = sel_conf[:,0]
                sel_conf = sel_conf.astype(int)
                for ii in sel_conf :
                    conf_name = os.path.join(tt, "traj")
                    conf_name = os.path.join(conf_name, str(ii) + '.lammpstrj')
                    conf_name = os.path.abspath(conf_name)
                    fp_task_name = make_fp_task_name(int(ss), cc)
                    fp_task_path = os.path.join(work_path, fp_task_name)
                    create_path(fp_task_path)
                    fp_tasks.append(fp_task_path)
                    cwd = os.getcwd()
                    os.chdir(fp_task_path)
                    os.symlink(os.path.relpath(conf_name), 'conf.lmp')
                    incar = make_vasp_incar(ecut, ediff, npar, kpar, kspacing = kspacing, kgamma = True)
                    with open('INCAR', 'w') as fp:
                        fp.write(incar)
                    for pair in fp_link_files :
                        os.symlink(pair[0], pair[1])
                    os.chdir(cwd)
                    cc += 1
                    count_total += 1
                    if count_total > fp_task_max :
                        return fp_tasks
    return fp_tasks

def make_fp_vasp (iter_index, 
                  jdata) :
    fp_task_max = jdata['fp_task_max']
    fp_params = jdata['fp_params']
    fp_link_files = jdata['fp_link_files']

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    create_path(work_path)
    modd_path = os.path.join(iter_name, model_devi_name)
    modd_task = glob.glob(os.path.join(modd_path, "task.*"))
    system_index = []
    for ii in modd_task :        
        system_index.append(os.path.basename(ii).split('.')[1])
    system_index.sort()
    set_tmp = set(system_index)
    system_index = list(set_tmp)

    fp_tasks = _make_fp_vasp_inner(modd_path, work_path, system_index, fp_task_max, fp_link_files, fp_params)
    if len(fp_tasks) == 0 :
        return

    command = os.path.join(os.getcwd(), "lib/ovito_file_convert.py")
    command += " conf.lmp POSCAR"
    exec_hosts(MachineLocal, command, 1, fp_tasks, verbose = True)


def make_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        make_fp_vasp(iter_index, jdata) 
    else :
        raise RuntimeError ("unsupported fp style")

def run_fp_vasp (iter_index,
                 jdata,
                 exec_machine) :
    fp_command = jdata['fp_command']
    fp_np = jdata['fp_np']
#    fp_command = ("OMP_NUM_THREADS=1 mpirun -n %d " % fp_np) + fp_command
    fp_command = cmd_append_log(fp_command, "vasp.log")

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    fp_run_tasks = []
    for ii in fp_tasks :
        if os.path.isfile(os.path.join(ii, 'OUTCAR')) :
            with open(os.path.join(ii, 'OUTCAR'), 'r') as fp :
                content = fp.read()
                count = content.count('TOTAL-FORCE')
                if count != 1 :
                    fp_run_tasks.append(ii)
        else :
            fp_run_tasks.append(ii)

    exec_hosts_batch(exec_machine, fp_command, fp_np, fp_run_tasks, verbose = True, mpi = True)
        
def run_fp (iter_index,
            jdata,
            exec_machine) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        run_fp_vasp(iter_index, jdata, exec_machine) 
    else :
        raise RuntimeError ("unsupported fp style")    


def post_fp_vasp (iter_index,
                  jdata):
    model_devi_jobs = jdata['model_devi_jobs']
    job_names = get_job_names (model_devi_jobs)
    assert (iter_index < len(job_names)) 
    cur_job_name = job_names[iter_index]    
    cur_job = model_devi_jobs[cur_job_name]

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

    to_config = 'template/tools.vasp/cessp2force_lin.py'
    to_config = os.path.abspath(to_config)
    cmd_to_config = to_config + " OUTCAR "
    cmd_to_config = cmd_append_log(cmd_to_config, "to_config.log")
    exec_hosts(MachineLocal, cmd_to_config, 1, fp_tasks)

    convert_to_raw = 'template/tools.vasp/convert2raw.py'
    shuffle_raw = 'template/tools.raw/shuffle_raw.py'
    copy_raw = 'template/tools.raw/copy_raw.py'
    raw_to_set = 'template/tools.raw/raw_to_set.sh'
    convert_to_raw = os.path.abspath(convert_to_raw)
    shuffle_raw = os.path.abspath(shuffle_raw)
    copy_raw = os.path.abspath(copy_raw)
    raw_to_set = os.path.abspath(raw_to_set)
    cwd = os.getcwd()
    for ss in system_index :
        sys_config_data = glob.glob(os.path.join(work_path, "task.%s.*/test.configs"%ss))
        sys_data_path = os.path.join(work_path, 'data.%s/orig'%ss)
        create_path(sys_data_path)
        with open(os.path.join(sys_data_path, 'data.configs'), 'wb') as wfd:
            for f in sys_config_data :
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)
        os.chdir(sys_data_path)
        sp.check_call(convert_to_raw + ' data.configs', shell = True)
        os.chdir('..')
        sp.check_call(shuffle_raw + ' orig/ .', shell = True)
        if os.path.isfile('type.raw') :
            os.remove('type.raw')
        os.symlink('orig/type.raw', 'type.raw')
        sp.check_call(raw_to_set, shell = True)
        os.chdir(cwd)

def post_fp (iter_index,
             jdata) :
    fp_style = jdata['fp_style']

    if fp_style == "vasp" :
        post_fp_vasp(iter_index, jdata) 
    else :
        raise RuntimeError ("unsupported fp style")            
    
def run_iter (json_file, exec_machine) :
    fp = open (json_file, 'r')
    jdata = json.load (fp)
    numb_iter = jdata["numb_iter"]
    max_tasks = 10000
    numb_task = 9
    record = "record.dpgen"

    iter_rec = [0, -1]
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                iter_rec = [int(x) for x in line.split()]
        logging.info ("continue from iter %03d task %02d" % (iter_rec[0], iter_rec[1]))

    for ii in range (numb_iter) :
        for jj in range (numb_task) :
            if ii * max_tasks + jj <= iter_rec[0] * max_tasks + iter_rec[1] : 
                continue
            if   jj == 0 :
                log_iter ("make_train", ii, jj)
                make_train (ii, jdata) 
            elif jj == 1 :
                log_iter ("run_train", ii, jj)
                run_train  (ii, jdata, exec_machine)
            elif jj == 2 :
                log_iter ("post_train", ii, jj)
                post_train  (ii, jdata)
            elif jj == 3 :
                log_iter ("make_model_devi", ii, jj)
                make_model_devi  (ii, jdata)
            elif jj == 4 :
                log_iter ("run_model_devi", ii, jj)
                run_model_devi  (ii, jdata, exec_machine)
            elif jj == 5 :
                log_iter ("post_model_devi", ii, jj)
                post_model_devi  (ii, jdata)
            elif jj == 6 :
                log_iter ("make_fp", ii, jj)
                make_fp (ii, jdata)
            elif jj == 7 :
                log_iter ("run_fp", ii, jj)
                run_fp (ii, jdata, exec_machine)
            elif jj == 8 :
                log_iter ("post_fp", ii, jj)
                post_fp (ii, jdata)
            else :
                raise RuntimeError ("unknow task %d, something wrong" % jj)
            record_iter (record, ii, jj)


def _main () :
    parser = argparse.ArgumentParser()
    parser.add_argument("JSON", type=str, 
                        help="The json parameter")
    parser.add_argument("--machine", type=str, 
                        help="The machine settings")        
    args = parser.parse_args()

    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')
    # logging.basicConfig (filename="compute_string.log", filemode="a", level=logging.INFO, format='%(asctime)s %(message)s')

    machine_type = "local"    
    gmxrc = None
    vcores = None
    if args.machine != None :
        fp = open (args.machine, 'r')
        jdata = json.load (fp)
        machine_type = jdata["machine_type"]
        gmxrc = jdata["gmxrc"]
        vcores = jdata["virtual_cores"]

    if   machine_type == "local" :
        exec_machine = MachineLocal
    elif machine_type == "slurm" :
        exec_machine = MachineSlurm

    if vcores != None:
        exec_machine.has_virtual_cores(vcores)
    if gmxrc != None:
        exec_machine.add_source_file(gmxrc)

    logging.info ("start running")
    run_iter (args.JSON, exec_machine)
    logging.info ("finished!")

if __name__ == '__main__':
    # _main()
    fp = open ('param.json', 'r')
    jdata = json.load (fp)
    # post_train(0, jdata)
    logging.basicConfig (level=logging.INFO, format='%(asctime)s %(message)s')
    run_iter('param.json', MachinePBS)
    # make_model_devi(0, jdata, None)
    # run_model_devi(0, jdata, MachineLocal)
    # post_model_devi(0, jdata)
    # make_fp(0, jdata)
    # run_fp(0, jdata, MachineLocal)
#    post_fp(0, jdata)
#    make_train(1, jdata)
#    run_train(1, jdata, MachineLocal)

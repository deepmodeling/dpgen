#!/usr/bin/env python3

"""
init: crystal configuration
task:
        00.equi
        01.eos
        02.elastic
        03.vacancy
        04.interstitial
        05.surf
        06.phonon
"""


import sys
import os, re, argparse, filecmp, json, glob
import lib.vasp as vasp
import lib.lammps as lammps
import random
import logging
import warnings
import shutil
import time
import numpy as np
import subprocess as sp
from lib.utils import make_iter_name
from lib.utils import create_path
from lib.utils import copy_file_list
from lib.utils import replace
from lib.utils import cmd_append_log
from lib.utils import log_iter
from lib.utils import record_iter
from lib.utils import log_iter
from lib.pwscf import make_pwscf_input
import lib.MachineLocal as MachineLocal
import lib.MachineLocalGPU as MachineLocalGPU
import lib.MachineSlurm as MachineSlurm
import lib.MachinePBS as MachinePBS
from lib.machine_exec import exec_hosts
from lib.machine_exec import exec_hosts_batch
from lib.batch_exec import exec_batch
from lib.batch_exec import exec_batch_group
from lib.RemoteJob import SSHSession, JobStatus, SlurmJob, PBSJob, CloudMachineJob
import gen_00_equi,cmpt_00_equi
import gen_01_eos,cmpt_01_eos
import gen_02_elastic,cmpt_02_elastic
import gen_03_vacancy,cmpt_03_vacancy
import gen_04_interstitial,cmpt_04_interstitial
import gen_05_surf,cmpt_05_surf
import requests
from hashlib import sha1

def _verfy_ac(private_key, params):
    items= sorted(params.items())
    
    params_data = ""
    for key, value in items:
        params_data = params_data + str(key) + str(value)
    params_data = params_data + private_key
    sign = sha1()
    sign.update(params_data.encode())
    signature = sign.hexdigest()
    return signature

def _ucloud_remove_machine(machine, UHostId):
    ucloud_url = machine['url']
    ucloud_stop_param = {}
    ucloud_stop_param['Action'] = "StopUHostInstance"
    ucloud_stop_param['Region'] = machine['ucloud_param']['Region']
    ucloud_stop_param['UHostId'] = UHostId
    ucloud_stop_param['PublicKey'] = machine['ucloud_param']['PublicKey']
    ucloud_stop_param['Signature'] = _verfy_ac(machine['Private'], ucloud_stop_param)

    
    req = requests.get(ucloud_url, ucloud_stop_param)
    if req.json()['RetCode'] != 0 :
        raise RuntimeError ("failed to stop ucloud machine")

    terminate_fin = False
    try_time = 0 
    while not terminate_fin:
        ucloud_delete_param = {}
        ucloud_delete_param['Action'] = "TerminateUHostInstance"
        ucloud_delete_param['Region'] = machine['ucloud_param']['Region']
        ucloud_delete_param['UHostId'] = UHostId
        ucloud_delete_param['PublicKey'] = machine['ucloud_param']['PublicKey']
        ucloud_delete_param['Signature'] = _verfy_ac(machine['Private'], ucloud_delete_param)    
        req = requests.get(ucloud_url, ucloud_delete_param)
        if req.json()['RetCode'] == 0 :
            terminate_fin = True
        try_time = try_time + 1
        if try_time >= 200:
            raise RuntimeError ("failed to terminate ucloud machine")
        time.sleep(10)
    print("Machine ",UHostId,"has been successfully terminated!")   

def _ucloud_submit_jobs(machine,
                        resources,
                        command,
                        work_path,
                        tasks,
                        group_size,
                        forward_common_files,
                        forward_task_files,
                        backward_task_files) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    njob = len(task_chunks)
    continue_status = False
    if os.path.isfile("record.machine"):
        with open ("record.machine", "r") as fr:
            record_machine = json.load(fr)
            if record_machine["purpose"] == machine["purpose"] and record_machine["njob"] == njob:
                continue_status = True
                ucloud_machines = record_machine["ucloud_machines"]
                ucloud_hostids = record_machine["ucloud_hostids"]
        fr.close()
    ucloud_url = machine['url']
    if continue_status == False:
        assert machine['machine_type'] == 'ucloud'
        ucloud_start_param = machine['ucloud_param']
        ucloud_start_param['Action'] = "CreateUHostInstance"
        ucloud_start_param['Name'] = "train"
        ucloud_start_param['Signature'] = _verfy_ac(machine['Private'], ucloud_start_param)

        
        ucloud_machines = []
        ucloud_hostids = []
        for ii in range(njob) :
            req = requests.get(ucloud_url, ucloud_start_param)
            if req.json()['RetCode'] != 0 :
                print(json.dumps(req.json(),indent=2, sort_keys=True))
                raise RuntimeError ("failed to start ucloud machine")
            ucloud_machines.append(str(req.json()["IPs"][0]))
            ucloud_hostids.append(str(req.json()["UHostIds"][0]))

        new_record_machine = {}
        new_record_machine["purpose"] = machine["purpose"]
        new_record_machine["njob"] = njob
        new_record_machine["ucloud_machines"] = ucloud_machines
        new_record_machine["ucloud_hostids"] = ucloud_hostids
        with open ("record.machine", "w") as fw:
            json.dump(new_record_machine, fw)
        fw.close()

    machine_fin = [False for ii in ucloud_machines]
    total_machine_num = len(ucloud_machines)
    fin_machine_num = 0
    while not all(machine_fin):
        for idx,mac in enumerate(ucloud_machines):
            if not machine_fin[idx]:
                ucloud_check_param = {}
                ucloud_check_param['Action'] = "GetUHostInstanceVncInfo"
                ucloud_check_param['Region'] = machine['ucloud_param']['Region']
                ucloud_check_param['UHostId'] = ucloud_hostids[idx]
                ucloud_check_param['PublicKey'] = machine['ucloud_param']['PublicKey']
                ucloud_check_param['Signature'] = _verfy_ac(machine['Private'], ucloud_check_param)
                req = requests.get(ucloud_url, ucloud_check_param)
                print("the UHostId is", ucloud_hostids[idx])
                print(json.dumps(req.json(),indent=2, sort_keys=True))
                if req.json()['RetCode'] == 0 :
                    machine_fin[idx] = True
                    fin_machine_num = fin_machine_num + 1
        print("Current finish",fin_machine_num,"/", total_machine_num)

        
        ucloud_check_param1 = {}
        ucloud_check_param1['Action'] = "DescribeUHostInstance"
        ucloud_check_param1['Region'] = machine['ucloud_param']['Region']
        ucloud_check_param1["Limit"] = 100
        ucloud_check_param1['PublicKey'] = machine['ucloud_param']['PublicKey']
        ucloud_check_param1['Signature'] = _verfy_ac(machine['Private'], ucloud_check_param1)
        req1 = requests.get(ucloud_url, ucloud_check_param1).json()
        
        machine_all_fin = True
        for idx1 in range(int(req1["TotalCount"])):
            if req1["UHostSet"][idx1]["State"] != "Running":
                machine_all_fin = False
                break
        if machine_all_fin == True:
            machine_fin = [True for i in machine_fin]
        time.sleep(10)
    ssh_sess = []
    ssh_param = {}
    ssh_param['port'] = 22
    ssh_param['username'] = 'root'
    ssh_param['work_path'] = machine['work_path']
    for ii in ucloud_machines :
        ssh_param['hostname'] = ii
        ssh_sess.append(SSHSession(ssh_param))

    job_list = []
    for ii in range(njob) :
        chunk = task_chunks[ii]
        print("Current machine is", ucloud_machines[ii])
        rjob = CloudMachineJob(ssh_sess[ii], work_path)
        rjob.upload('.',  forward_common_files)
        rjob.upload(chunk, forward_task_files)
        rjob.submit(chunk, command, resources = resources)
        job_list.append(rjob)
    
    job_fin = [False for ii in job_list]
    while not all(job_fin) :
        for idx,rjob in enumerate(job_list) :
            if not job_fin[idx] :
                status = rjob.check_status()
                if status == JobStatus.terminated :
                    raise RuntimeError("find unsuccessfully terminated job on machine" % ucloud_machines[idx])
                elif status == JobStatus.finished :
                    rjob.download(task_chunks[idx], backward_task_files)
                    rjob.clean()
                    _ucloud_remove_machine(machine, ucloud_hostids[idx])
                    job_fin[idx] = True
        time.sleep(10)
    os.remove("record.machine")

def _group_slurm_jobs(ssh_sess,
                      resources,
                      command,
                      work_path,
                      tasks,
                      group_size,
                      forward_common_files,
                      forward_task_files,
                      backward_task_files,
                      remote_job = SlurmJob) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    job_list = []
    for chunk in task_chunks :
        rjob = remote_job(ssh_sess, work_path)
        cwd=[os.path.basename(ii) for ii in chunk]
        rjob.upload(cwd,forward_common_files)
        rjob.upload(chunk, forward_task_files)
        rjob.submit(chunk, command, resources = resources)
        job_list.append(rjob)

    job_fin = [False for ii in job_list]
    while not all(job_fin) :
        for idx,rjob in enumerate(job_list) :
            if not job_fin[idx] :
                status = rjob.check_status()
                if status == JobStatus.terminated :
                    raise RuntimeError("find unsuccessfully terminated job in %s" % rjob.get_job_root())
                elif status == JobStatus.finished :
                    rjob.download(task_chunks[idx], backward_task_files)
                    rjob.clean()
                    job_fin[idx] = True
        time.sleep(10)

def _group_local_jobs(ssh_sess,
                      resources,
                      command,
                      work_path,
                      tasks,
                      group_size,
                      forward_common_files,
                      forward_task_files,
                      backward_task_files) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    job_list = []
    for chunk in task_chunks :
        rjob = CloudMachineJob(ssh_sess, work_path)
        cwd=[os.path.basename(ii) for ii in chunk]
        rjob.upload(cwd,forward_common_files)
        rjob.upload(chunk, forward_task_files)
        rjob.submit(chunk, command, resources = resources)
        job_list.append(rjob)
        job_fin = False
        while not job_fin :
            status = rjob.check_status()
            if status == JobStatus.terminated :
                raise RuntimeError("find unsuccessfully terminated job in %s" % rjob.get_job_root())
            elif status == JobStatus.finished :
                rjob.download(chunk, backward_task_files)
                rjob.clean()
                job_fin = True
            time.sleep(10)

def _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files):

    print("group_size",group_size)
    if ssh_sess == None and machine_type == 'ucloud':
        print("The first situation!")
        _ucloud_submit_jobs(machine,
                            resources,
                            command,
                            work_path,
                            run_tasks,
                            group_size,
                            model_names,
                            forward_files,
                            backward_files)
    elif machine_type == 'slurm' :        
        print("The second situation!")
        _group_slurm_jobs(ssh_sess,
                           resources,
                           command,
                           work_path,
                           run_tasks,
                           group_size,
                           model_names,
                           forward_files,
                           backward_files)
    elif machine_type == 'pbs' :        
        _group_slurm_jobs(ssh_sess,
                           resources,
                           command,
                           work_path,
                           run_tasks,
                           group_size,
                           model_names,
                           forward_files,
                           backward_files,
                          remote_job = PBSJob)
    elif machine_type == 'local' :        
        _group_local_jobs(ssh_sess,
                           resources,
                           command,
                           work_path,
                           run_tasks,
                           group_size,
                           model_names,
                           forward_files,
                           backward_files)
    else :
        raise RuntimeError("unknow machine type")

def gen_equi(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cwd=os.getcwd()
    #vasp
    if task_type=="vasp":
        gen_00_equi.make_vasp(jdata, conf_dir)
    #deepmd
    elif task_type=="deepmd":
        gen_00_equi.make_deepmd_lammps (jdata, conf_dir)
    #meam
    elif task_type=="meam":
        gen_00_equi.make_meam_lammps (jdata, conf_dir)
    else :
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    os.chdir(cwd)
        
def run_equi(task_type,jdata,mdata,ssh_sess):
        #rmprint("This module has been run !")

    conf_dir=jdata['conf_dir']
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)

    conf_path = os.path.abspath(conf_dir)
    equi_path = re.sub('confs', '00.equi', conf_path)
    if task_type=="vasp":
        work_path=os.path.join(equi_path, 'vasp-k%.2f' % kspacing)
    elif task_type=="deepmd":
        work_path=os.path.join(equi_path, 'deepmd')
    elif task_type=="meam":
        work_path=os.path.join(equi_path, 'meam')
    assert(os.path.isdir(work_path))
    
    all_task = glob.glob(os.path.join(work_path,'.'))
    
    #vasp
    if task_type=="vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        machine_type = mdata['fp_machine']['machine_type']
        command = vasp_exec
        command = cmd_append_log(command, "log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'OUTCAR')
            if os.path.isfile(fres) :
                if not vasp.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR','CONTCAR']
        model_names=[]

    #lammps
    elif task_type=="deepmd" or task_type=="meam":
        lmp_exec = mdata['lmp_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        machine_type = mdata['model_devi_machine']['machine_type']
        command = lmp_exec + " -i lammps.in"
        command = cmd_append_log(command, "model_devi.log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                with open(fres, 'r') as fp :
                    lines = fp.read().split('\n')
                flag=False
                for jj in lines:
                    if ("Final energy per atoms" in jj) and (not 'print' in jj):
                        flag=True
                if not flag:
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)
        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['conf.lmp', 'lammps.in']
        backward_files = ['dump.relax','log.lammps','model_devi.out', 'model_devi.log']
        all_models = glob.glob(os.path.join(deepmd_model_dir, '*.pb'))
        model_names = [os.path.basename(ii) for ii in all_models]
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    
    _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files)

def cmpt_equi(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    stable=jdata['store_stable']
    #vasp
    if task_type=="vasp":
        vn, ve, vv = cmpt_00_equi.comput_vasp_nev(jdata, conf_dir, stable)
        print('conf_dir:\t EpA(eV)  VpA(A^3)')
        print("%s\t %8.4f  %7.3f " % (conf_dir, ve, vv))
    #deepmd
    elif task_type=="deepmd":
        ln, le, lv = cmpt_00_equi.comput_lmp_nev(conf_dir, 'deepmd', stable)
        print('conf_dir:\t EpA(eV)  VpA(A^3)')
        print("%s\t %8.4f  %7.3f " % (conf_dir, le, lv))
    #meam
    elif task_type=="meam":
        ln, le, lv = cmpt_00_equi.comput_lmp_nev(conf_dir, 'meam', stable)
        print('conf_dir:\t EpA(eV)  VpA(A^3)')
        print("%s\t %8.4f  %7.3f " % (conf_dir, le, lv))
    else :
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

def gen_eos(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    fix_shape=jdata['fix_shape']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_01_eos.make_vasp(jdata, conf_dir)  
    #deepmd             
    elif task_type == "deepmd" :
        if fix_shape :
            gen_01_eos.make_deepmd_lammps_fixv(jdata, conf_dir)
        else :
            gen_01_eos.make_deepmd_lammps(jdata, conf_dir)        
    #meam
    elif task_type == "meam" :
        if fix_shape :
            gen_01_eos.make_meam_lammps_fixv(jdata, conf_dir)
        else :
            raise RuntimeError("not implemented ", 'meam')            
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_eos(task_type,jdata,mdata,ssh_sess):
    conf_dir=jdata['conf_dir']
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    
    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', '01.eos', conf_path)
    if task_type=="vasp":
        work_path=os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    elif task_type=="deepmd":
        work_path=os.path.join(task_path, 'deepmd')
    elif task_type=="meam":
        work_path=os.path.join(task_path, 'meam')
    assert(os.path.isdir(work_path))
    print(work_path)
    
    all_task = glob.glob(os.path.join(work_path, "vol-*"))
    all_task.sort()
    
    #vasp
    if task_type=="vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        machine_type = mdata['fp_machine']['machine_type']
        command = vasp_exec
        command = cmd_append_log(command, "log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'OUTCAR')
            if os.path.isfile(fres) :
                if not vasp.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR']
        model_names=[]

    #lammps
    elif task_type=="deepmd" or task_type=="meam":
        lmp_exec = mdata['lmp_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        machine_type = mdata['model_devi_machine']['machine_type']
        command = lmp_exec + " -i lammps.in"
        command = cmd_append_log(command, "model_devi.log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                with open(fres, 'r') as fp :
                    lines = fp.read().split('\n')
                flag=False
                for jj in lines:
                    if ("Final energy per atoms" in jj) and (not 'print' in jj):
                        flag=True
                if not flag:
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['conf.lmp', 'lammps.in']
        backward_files = ['log.lammps','model_devi.out', 'model_devi.log']
        all_models = glob.glob(os.path.join(deepmd_model_dir, '*.pb'))
        model_names = [os.path.basename(ii) for ii in all_models]
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files)

def cmpt_eos(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    #vasp
    if task_type == "vasp":
        cmpt_01_eos.comput_vasp_eos(jdata, conf_dir)   
    #deepmd                
    elif task_type == "deepmd" :
        cmpt_01_eos.comput_lmp_eos(conf_dir, 'deepmd')
    #meam
    elif task_type == "meam" :
        cmpt_01_eos.comput_lmp_eos(conf_dir, 'meam')
    else :
        raise RuntimeError("unknow task ", task_type)

def gen_elastic(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_02_elastic.make_vasp(jdata, conf_dir)
    #deepmd
    elif task_type == "deepmd":
        gen_02_elastic.make_deepmd_lammps (jdata, conf_dir)
    #meam
    elif task_type == "meam":
        gen_02_elastic.make_meam_lammps (jdata, conf_dir)
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    os.chdir(cwd)

def run_elastic(task_type,jdata,mdata,ssh_sess):
    conf_dir=jdata['conf_dir']
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)
    
    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', '02.elastic', conf_path)
    if task_type == "vasp":
        work_path=os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    elif task_type == "deepmd":
        work_path=os.path.join(task_path, 'deepmd')
    elif task_type == "meam":
        work_path=os.path.join(task_path, 'meam')
    assert(os.path.isdir(work_path))
    print(work_path)
    
    all_task = glob.glob(os.path.join(work_path, "dfm-*"))
    all_task.sort()
    
    #vasp
    if task_type == "vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        machine_type = mdata['fp_machine']['machine_type']
        command = vasp_exec
        command = cmd_append_log(command, "log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'OUTCAR')
            if os.path.isfile(fres) :
                if not vasp.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)
        
        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['INCAR', 'POSCAR','POTCAR','KPOINTS']
        backward_files = ['OUTCAR','CONTCAR']
        model_names=[]

    #lammps
    elif task_type == "deepmd" or task_type == "meam":
        lmp_exec = mdata['lmp_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        machine_type = mdata['model_devi_machine']['machine_type']
        command = lmp_exec + " -i lammps.in"
        command = cmd_append_log(command, "model_devi.log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                with open(fres, 'r') as fp :
                    lines = fp.read().split('\n')
                flag=False
                for jj in lines:
                    if ('Final Stress' in jj) and (not 'print' in jj):
                        flag=True
                if not flag:
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['conf.lmp', 'lammps.in','strain.out']
        backward_files = ['log.lammps','model_devi.out', 'model_devi.log']
        all_models = glob.glob(os.path.join(deepmd_model_dir, '*.pb'))
        model_names = [os.path.basename(ii) for ii in all_models]
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files)
    
def cmpt_elastic(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    if task_type == "vasp":
        cmpt_02_elastic.cmpt_vasp(jdata, conf_dir)               
    elif task_type == "deepmd":
        cmpt_02_elastic.cmpt_deepmd_lammps(jdata, conf_dir, 'deepmd')
    elif task_type == "meam":
        cmpt_02_elastic.cmpt_deepmd_lammps(jdata, conf_dir, 'meam')
    else :
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    
def gen_vacancy(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    supercell=jdata['supercell']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_03_vacancy.make_vasp(jdata, conf_dir, supercell)
    #deepmd
    elif task_type == "deepmd":
        gen_03_vacancy.make_deepmd_lammps(jdata, conf_dir, supercell)
    #meam
    elif task_type == "meam":
        gen_03_vacancy.make_meam_lammps(jdata, conf_dir, supercell)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_vacancy(task_type,jdata,mdata,ssh_sess):
    conf_dir=jdata['conf_dir']
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', '03.vacancy', conf_path)
    if task_type == "vasp":
        work_path=os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    elif task_type == "deepmd":
        work_path=os.path.join(task_path, 'deepmd')
    elif task_type == "meam":
        work_path=os.path.join(task_path, 'meam')
    assert(os.path.isdir(work_path))
    
    all_task = glob.glob(os.path.join(work_path,'struct-*'))
    
    #vasp
    if task_type == "vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        machine_type = mdata['fp_machine']['machine_type']
        command = vasp_exec
        command = cmd_append_log(command, "log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'OUTCAR')
            if os.path.isfile(fres) :
                if not vasp.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)
        
        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR']
        model_names=[]

    #lammps
    elif task_type == "deepmd" or task_type == "meam":
        lmp_exec = mdata['lmp_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        machine_type = mdata['model_devi_machine']['machine_type']
        command = lmp_exec + " -i lammps.in"
        command = cmd_append_log(command, "model_devi.log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                with open(fres, 'r') as fp :
                    lines = fp.read().split('\n')
                flag=False
                for jj in lines:
                    if ("Final energy per atoms" in jj) and (not 'print' in jj):
                        flag=True
                if not flag:
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['conf.lmp', 'lammps.in']
        backward_files = ['log.lammps','model_devi.out', 'model_devi.log']
        all_models = glob.glob(os.path.join(deepmd_model_dir, '*.pb'))
        model_names = [os.path.basename(ii) for ii in all_models]
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files)

def cmpt_vacancy(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    supercell=jdata['supercell']
    #vasp
    if task_type == "vasp":
        cmpt_03_vacancy.cmpt_vasp(jdata, conf_dir, supercell)               
    #deepmd
    elif task_type == "deepmd":
        cmpt_03_vacancy.cmpt_deepmd_lammps(jdata, conf_dir, supercell, 'deepmd')
    #meam
    elif task_type == "meam":
        cmpt_03_vacancy.cmpt_deepmd_lammps(jdata, conf_dir, supercell,'meam')
    else :
        raise RuntimeError("unknow task ", task_type)

def gen_interstitial(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    supercell=jdata['supercell']
    insert_ele=jdata['insert_ele']
    reprod_opt=jdata['reprod-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_04_interstitial.make_vasp(jdata, conf_dir, supercell, insert_ele)
    #deepmd
    elif task_type == "deepmd":
        if not reprod_opt:
            gen_04_interstitial.make_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_type)
        else :
            gen_04_interstitial.make_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type)
    #meam
    elif task_type == "meam" :
        if not reprod_opt:
            gen_04_interstitial.make_meam_lammps(jdata, conf_dir, supercell, insert_ele, task_type)
        else :
            gen_04_interstitial.make_meam_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_interstitial(task_type,jdata,mdata,ssh_sess):
    conf_dir=jdata['conf_dir']
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', '04.interstitial', conf_path)
    if task_type == "vasp":
        work_path=os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    elif task_type == "deepmd":
        work_path=os.path.join(task_path, 'deepmd')
    elif task_type == "meam":
        work_path=os.path.join(task_path, 'meam')
    assert(os.path.isdir(work_path))
    
    all_task = glob.glob(os.path.join(work_path,'struct-*'))
    
    #vasp
    if task_type == "vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        machine_type = mdata['fp_machine']['machine_type']
        command = vasp_exec
        command = cmd_append_log(command, "log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'OUTCAR')
            if os.path.isfile(fres) :
                if not vasp.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)
        
        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR','XDATCAR']
        model_names=[]

    #lammps
    elif task_type == "deepmd" or task_type == "meam":
        lmp_exec = mdata['lmp_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        machine_type = mdata['model_devi_machine']['machine_type']
        command = lmp_exec + " -i lammps.in"
        command = cmd_append_log(command, "model_devi.log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                with open(fres, 'r') as fp :
                    lines = fp.read().split('\n')
                flag=False
                for jj in lines:
                    if ("Final energy per atoms" in jj) and (not 'print' in jj):
                        flag=True
                if not flag:
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['conf.lmp', 'lammps.in']
        backward_files = ['log.lammps','model_devi.out', 'model_devi.log']
        all_models = glob.glob(os.path.join(deepmd_model_dir, '*.pb'))
        model_names = [os.path.basename(ii) for ii in all_models]
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files)

def cmpt_interstitial(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    supercell=jdata['supercell']
    insert_ele=jdata['insert_ele']
    reprod_opt=jdata['reprod-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        cmpt_04_interstitial.cmpt_vasp(jdata, conf_dir, supercell, insert_ele)
    #deepmd
    elif task_type == "deepmd":
        if not reprod_opt:
            cmpt_04_interstitial.cmpt_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_type)
        else :
            cmpt_04_interstitial.cmpt_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type)
    #meam
    elif task_type == "meam" :
        if not reprod_opt:
            cmpt_04_interstitial.cmpt_meam_lammps(jdata, conf_dir, supercell, insert_ele, task_type)
        else :
            cmpt_04_interstitial.cmpt_meam_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def gen_surf(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    max_miller=jdata['max_miller']
    relax_box=jdata['relax_box']
    static_opt=jdata['static-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_05_surf.make_vasp(jdata, conf_dir, max_miller, static = static_opt, relax_box = relax_box)
    #deepmd
    elif task_type == "deepmd" :
        gen_05_surf.make_deepmd_lammps(jdata, conf_dir, max_miller, static = static_opt, relax_box = relax_box, task_name = 'deepmd')
    #meam
    elif task_type == "meam" :
        gen_05_surf.make_meam_lammps(jdata, conf_dir, max_miller, static = static_opt, relax_box = relax_box, task_name = 'meam')
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_surf(task_type,jdata,mdata,ssh_sess):
    conf_dir=jdata['conf_dir']
    fp_params = jdata['vasp_params']
    kspacing = fp_params['kspacing']
    deepmd_model_dir = jdata['deepmd_model_dir']
    deepmd_model_dir = os.path.abspath(deepmd_model_dir)

    conf_path = os.path.abspath(conf_dir)
    task_path = re.sub('confs', '05.surf', conf_path)
    if task_type == "vasp":
        work_path=os.path.join(task_path, 'vasp-k%.2f' % kspacing)
    elif task_type == "deepmd":
        work_path=os.path.join(task_path, 'deepmd')
    elif task_type == "meam":
        work_path=os.path.join(task_path, 'meam')
    assert(os.path.isdir(work_path))
    
    all_task = glob.glob(os.path.join(work_path,'struct-*'))
     
    #vasp
    if task_type == "vasp":
        vasp_exec=mdata['fp_command']
        group_size = mdata['fp_group_size']
        resources = mdata['fp_resources']
        machine=mdata['fp_machine']
        machine_type = mdata['fp_machine']['machine_type']
        command = vasp_exec
        command = cmd_append_log(command, "log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'OUTCAR')
            if os.path.isfile(fres) :
                if not vasp.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)
            
        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR']
        model_names=[]

    #lammps
    elif task_type == "deepmd" or task_type == "meam":
        lmp_exec = mdata['lmp_command']
        group_size = mdata['model_devi_group_size']
        resources = mdata['model_devi_resources']
        machine=mdata['model_devi_machine']
        machine_type = mdata['model_devi_machine']['machine_type']
        command = lmp_exec + " -i lammps.in"
        command = cmd_append_log(command, "model_devi.log")

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                with open(fres, 'r') as fp :
                    lines = fp.read().split('\n')
                flag=False
                for jj in lines:
                    if ("Final energy per atoms" in jj) and (not 'print' in jj):
                        flag=True
                if not flag:
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]
        forward_files = ['conf.lmp', 'lammps.in']
        backward_files = ['log.lammps','model_devi.out', 'model_devi.log']
        all_models = glob.glob(os.path.join(deepmd_model_dir, '*.pb'))
        model_names = [os.path.basename(ii) for ii in all_models]
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    _run(machine,
         machine_type,
         ssh_sess,
         resources,
         command,
         work_path,
         run_tasks,
         group_size,
         model_names,
         forward_files,
         backward_files)

def cmpt_surf(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    static_opt=jdata['static-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        cmpt_05_surf.cmpt_vasp(jdata, conf_dir, static = static_opt) 
    #deepmd        
    elif task_type == "deepmd" :
        cmpt_05_surf.cmpt_deepmd_lammps(jdata, conf_dir, 'deepmd', static = static_opt)
    #meam  
    elif task_type == "meam" :
        cmpt_05_surf.cmpt_deepmd_lammps(jdata, conf_dir, 'meam', static = static_opt)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_task (json_file, machine_file) :
    with open (json_file, 'r') as fp :
        jdata = json.load (fp)
    with open (machine_file, 'r') as fp:
        mdata = json.load (fp)

    record = "record.auto_test"
    
    train_machine = mdata['train_machine']    
    if ('machine_type' in train_machine) and  \
       (train_machine['machine_type'] == 'ucloud'):
        train_ssh_sess = None
    else :
        train_ssh_sess = SSHSession(train_machine)
    
    model_devi_machine = mdata['model_devi_machine']    
    if ('machine_type' in model_devi_machine) and  \
       (model_devi_machine['machine_type'] == 'ucloud'):
        model_devi_ssh_sess = None
    else :
        model_devi_ssh_sess = SSHSession(model_devi_machine)
    
    fp_machine = mdata['fp_machine']    
    if ('machine_type' in fp_machine) and  \
       (fp_machine['machine_type'] == 'ucloud'):
        fp_ssh_sess = None
    else :
        fp_ssh_sess = SSHSession(fp_machine)
    
    confs = jdata['conf_dir']
    ii = jdata['task_type']
    jj=jdata['task']
    task_list=['equi','eos','elastic','vacancy','interstitial','surf','all']
    #default task
    log_iter ("gen_equi", ii, "equi")
    gen_equi (ii, jdata, mdata) 
    log_iter ("run_equi", ii, "equi")
    run_equi  (ii, jdata, mdata,model_devi_ssh_sess)
    log_iter ("cmpt_equi", ii,"equi")
    cmpt_equi (ii, jdata, mdata)
    if  jj == "eos" or jj=="all":
        log_iter ("gen_eos", ii, "eos")
        gen_eos (ii, jdata, mdata) 
        log_iter ("run_eos", ii, "eos")
        run_eos  (ii, jdata, mdata,model_devi_ssh_sess)
        log_iter ("cmpt_eos", ii, "eos")
        cmpt_eos (ii, jdata, mdata)
    if jj=="elastic" or jj=="all":
        log_iter ("gen_elastic", ii, "elastic")
        gen_elastic (ii, jdata, mdata) 
        log_iter ("run_elastic", ii, "elastic")
        run_elastic  (ii, jdata, mdata,model_devi_ssh_sess)
        log_iter ("cmpt_elastic", ii, "elastic")
        cmpt_elastic (ii, jdata, mdata)
    if jj=="vacancy" or jj=="all":
        log_iter ("gen_vacancy", ii, "vacancy")
        gen_vacancy (ii, jdata, mdata) 
        log_iter ("run_vacancy", ii, "vacancy")
        run_vacancy  (ii, jdata, mdata,model_devi_ssh_sess)
        log_iter ("cmpt_vacancy", ii, "vacancy")
        cmpt_vacancy (ii, jdata, mdata)
    if jj=="interstitial" or jj=="all":
        log_iter ("gen_interstitial", ii, "interstitial")
        gen_interstitial (ii, jdata, mdata) 
        log_iter ("run_interstitial", ii, "interstitial")
        run_interstitial  (ii, jdata, mdata,model_devi_ssh_sess)
        log_iter ("cmpt_interstitial", ii, "interstitial")
        cmpt_interstitial (ii, jdata, mdata)
    if jj=="surf" or jj=="all":
        log_iter ("gen_surf", ii, "surf")
        gen_surf (ii, jdata, mdata) 
        log_iter ("run_surf", ii, "surf")
        run_surf  (ii, jdata, mdata,model_devi_ssh_sess)
        log_iter ("cmpt_surf", ii, "surf")
        cmpt_surf (ii, jdata, mdata)
    if jj not in task_list :
        raise RuntimeError ("unknow task %s, something wrong" % jj)
    record_iter (record, confs, ii, jj)

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

    logging.info ("start auto-testing")
    run_task (args.PARAM, args.MACHINE)
    logging.info ("finished!")

if __name__ == '__main__':
    _main()


import os,sys,glob
import numpy as np
import subprocess as sp
from dpgen.remote.RemoteJob import SlurmJob

import requests
from hashlib import sha1

def _verfy_ac(private_key, params):
    items= sorted(params.items())
    
    params_data = "";
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

def ucloud_submit_jobs(machine,
                       resources,
                       command,
                       work_path,
                       tasks,
                       group_size,
                       forward_common_files,
                       forward_task_files,
                       backward_task_files, 
                       forward_task_deference = True) :
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
        rjob.upload(chunk, forward_task_files, 
                    dereference = forward_task_deference)
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


def group_slurm_jobs(ssh_sess,
                     resources,
                     command,
                     work_path,
                     tasks,
                     group_size,
                     forward_common_files,
                     forward_task_files,
                     backward_task_files,
                     remote_job = SlurmJob, 
                     forward_task_deference = True) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    job_list = []
    for chunk in task_chunks :
        rjob = remote_job(ssh_sess, work_path)
        rjob.upload('.',  forward_common_files)
        rjob.upload(chunk, forward_task_files, 
                    dereference = forward_task_deference)
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

def group_local_jobs(ssh_sess,
                     resources,
                     command,
                     work_path,
                     tasks,
                     group_size,
                     forward_common_files,
                     forward_task_files,
                     backward_task_files, 
                     forward_task_deference = True) :
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    job_list = []
    for chunk in task_chunks :
        rjob = CloudMachineJob(ssh_sess, work_path)
        rjob.upload('.',  forward_common_files)
        rjob.upload(chunk, forward_task_files, 
                    dereference = forward_task_deference)
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

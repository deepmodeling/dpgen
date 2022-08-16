#!/usr/bin/env python
# coding: utf-8

import os,sys,glob,time
import numpy as np
import subprocess as sp
from monty.serialization import dumpfn,loadfn
from dpgen.remote.RemoteJob import SlurmJob, PBSJob, CloudMachineJob, JobStatus, awsMachineJob,SSHSession
from dpgen import dlog

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

def aws_submit_jobs(machine,
                    resources,
                    command,
                    work_path,
                    tasks,
                    group_size,
                    forward_common_files,
                    forward_task_files,
                    backward_task_files, 
                    forward_task_deference = True):
    import boto3
    task_chunks = [
        [os.path.basename(j) for j in tasks[i:i + group_size]] \
        for i in range(0, len(tasks), group_size)
    ]
    task_chunks = (str(task_chunks).translate((str.maketrans('','',' \'"[]'))).split(',')) 
    # flatten the task_chunks
    print('task_chunks=',task_chunks)
    njob = len(task_chunks)    
    print('njob=',njob)
    continue_status = False
    ecs=boto3.client('ecs')
    ec2=boto3.client('ec2')
    status_list=[]    
    containerInstanceArns=ecs.list_container_instances(cluster="tensorflow")    
    if  containerInstanceArns['containerInstanceArns']:
        containerInstances=ecs.describe_container_instances(cluster="tensorflow", \
                                                    containerInstances=containerInstanceArns['containerInstanceArns'])['containerInstances']
        status_list=[container['status'] for container in containerInstances]

    need_apply_num=group_size-len(status_list)
    print('need_apply_num=',need_apply_num)
    if need_apply_num>0:
        for ii in range(need_apply_num) :   #apply for machines,
            ec2.run_instances(**machine['run_instances'])      
    machine_fin = False    
    status_list=[]
    while not len(status_list)>=group_size:
        containerInstanceArns=ecs.list_container_instances(cluster="tensorflow")
        if  containerInstanceArns['containerInstanceArns']:
            containerInstances=ecs.describe_container_instances(cluster="tensorflow", \
                                                        containerInstances=containerInstanceArns['containerInstanceArns'])['containerInstances']
            status_list=[container['status'] for container in containerInstances]            
        if len(status_list)>=group_size:
            break
        else:
            time.sleep(20)
    print('current available containers status_list=',status_list)       
    print('remote_root=',machine['remote_root'])   
    rjob = awsMachineJob(machine['remote_root'],work_path)
    taskARNs=[]
    taskstatus=[]
    running_job_num=0
    rjob.upload('.',  forward_common_files)
    for ijob in range(njob) :   #uplaod && submit job
        containerInstanceArns=ecs.list_container_instances(cluster="tensorflow")
        containerInstances=ecs.describe_container_instances(cluster="tensorflow", \
                                                    containerInstances=containerInstanceArns['containerInstanceArns'])['containerInstances']
        status_list=[container['status'] for container in containerInstances]
        print('current available containers status_list=',status_list)      
        while  running_job_num>=group_size:
            taskstatus=[task['lastStatus'] for task in ecs.describe_tasks(cluster='tensorflow',tasks=taskARNs)['tasks']]
            running_job_num=len(list(filter(lambda str:(str=='PENDING' or str =='RUNNING'),taskstatus)))
            print('waiting for running job finished, taskstatus=',taskstatus,'running_job_num=',running_job_num)
            time.sleep(10)         
        chunk = str(task_chunks[ijob])
        print('current task chunk=',chunk)
        task_definition=command['task_definition']
        concrete_command=(command['concrete_command'] %(work_path,chunk))
        command_override=command['command_override']
        command_override['containerOverrides'][0]['command'][0]=concrete_command
        print('concrete_command=',concrete_command)
        rjob.upload(chunk, forward_task_files, 
                    dereference = forward_task_deference)
        taskres=ecs.run_task(cluster='tensorflow',\
                     taskDefinition=task_definition,overrides=command_override)
        while not taskres['tasks'][0]:
            print('task submit failed,taskres=',taskres,'trying to re-submit'+str(chunk),)
            time.sleep(10)
            taskres=ecs.run_task(cluster='tensorflow',\
                     taskDefinition=task_definition,overrides=command_override)

        taskARNs.append(taskres['tasks'][0]['taskArn'])
        taskstatus=[task['lastStatus'] for task in ecs.describe_tasks(cluster='tensorflow',tasks=taskARNs)['tasks']]
        running_job_num=len(list(filter(lambda str:(str=='PENDING' or str =='RUNNING'),taskstatus)))
        print('have submitted %s/%s,taskstatus=' %(work_path,chunk) ,taskstatus,'running_job_num=',running_job_num )       
    task_fin_flag=False    
    while not task_fin_flag:       
        taskstatus=[task['lastStatus'] for task in ecs.describe_tasks(cluster='tensorflow',tasks=taskARNs)['tasks']]
        task_fin_flag=all([status=='STOPPED' for status in taskstatus])
        if task_fin_flag:
            print('task finished,next step:copy files to local && taskstatus=',taskstatus)
        else:
            print('all tasks submitted,task running && taskstatus=',taskstatus)
            time.sleep(20)    
    for ii in range(njob):
        chunk = task_chunks[ii]
        print('downloading '+str(chunk),backward_task_files)
        rjob.download(chunk,backward_task_files)

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
    cwd=os.getcwd()
    _pmap=PMap(cwd)
    path_map=_pmap.load()
    dlog.debug("work_path: %s"% work_path)
    dlog.debug("curr_path: %s"% cwd)

    job_list = []
    task_chunks_=['+'.join(ii) for ii in task_chunks]
    for ii in task_chunks_:
        dlog.debug("task_chunk %s" % ii)

    #dlog.debug(path_map)
    for ii,chunk in enumerate(task_chunks) :
       
        # map chunk info. to uniq id    
        chunk_uni=task_chunks_[ii].encode('utf-8')
        chunk_sha1=sha1(chunk_uni).hexdigest() 

        if chunk_sha1 in path_map:
           job_uuid=path_map[chunk_sha1][1].split('/')[-1]
           dlog.debug("load uuid %s" % job_uuid)
        else:
           job_uuid=None

        rjob = remote_job(ssh_sess, work_path, job_uuid)
        dlog.debug('uuid %s'%job_uuid)
        rjob.upload('.',  forward_common_files)
        rjob.upload(chunk, forward_task_files, 
                   dereference = forward_task_deference)
        if job_uuid:
           rjob.submit(chunk, command, resources = resources,restart=True)
        else:
           rjob.submit(chunk, command, resources = resources)
        job_list.append(rjob)
        path_map[chunk_sha1]=[rjob.local_root,rjob.remote_root]        
    _pmap.dump(path_map)
    
    job_fin = [False for ii in job_list]
    lcount=[0]*len(job_list)
    count_fail = 0
    while not all(job_fin) :
        for idx,rjob in enumerate(job_list) :
            if not job_fin[idx] :
                try:
                  status = rjob.check_status()
                except Exception:
                  ssh_sess = SSHSession(ssh_sess.remote_profile)
                  for _idx,_rjob in enumerate(job_list):
                    job_list[_idx] = SlurmJob(ssh_sess, work_path, _rjob.job_uuid)
                  count_fail = count_fail +1
                  dlog.info("ssh_sess failed for %d times"%count_fail)
                  break
                if status == JobStatus.terminated :
                    lcount[idx]+=1
                    _job_uuid=rjob.remote_root.split('/')[-1]
                    dlog.info('Job at %s  terminated, submit again'% _job_uuid)
                    dlog.debug('try %s times for %s'% (lcount[idx], _job_uuid))
                    rjob.submit(task_chunks[idx], command, resources = resources,restart=True)
                    if lcount[idx]>3:
                       dlog.info('Too many errors for ! %s ' % _job_uuid)
                       rjob.download(task_chunks[idx], backward_task_files,back_error=True)
                       rjob.clean()
                       job_fin[idx] = True
                elif status == JobStatus.finished :
                    rjob.download(task_chunks[idx], backward_task_files)
                    rjob.clean()
                    job_fin[idx] = True
        time.sleep(10)
    dlog.debug('error count') 
    dlog.debug(lcount)
    # delete path map file when job finish
    _pmap.delete()

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

class PMap(object):
   '''
   Path map class to operate {read,write,delte} the pmap.json file
   '''

   def __init__(self,path,fname="pmap.json"):
       self.f_path_map=os.path.join(path,fname)

   def load(self):
      f_path_map=self.f_path_map
      if os.path.isfile(f_path_map):
         path_map=loadfn(f_path_map)
      else:
         path_map={}
      return path_map

   def dump(self,pmap,indent=4):
      f_path_map=self.f_path_map
      dumpfn(pmap,f_path_map,indent=indent)

   def delete(self):
      f_path_map=self.f_path_map
      try:
         os.remove(f_path_map)
      except Exception:
         pass

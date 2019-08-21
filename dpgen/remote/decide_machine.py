#!/usr/bin/env python
# coding: utf-8
from dpgen.remote.RemoteJob import SSHSession, JobStatus, SlurmJob, PBSJob, LSFJob, CloudMachineJob
import os
import json
import numpy as np


def decide_train_machine(mdata):
	if 'train' in mdata:
	    continue_flag = False
	    ## decide whether to use an existing machine
	    if 'record.machine' in os.listdir():
	        try:
	            with open('record.machine', 'r') as _infile:
	                profile = json.load(_infile)
	                if profile['purpose'] == 'train':
	                    mdata['train_machine'] = profile['machine']
	                    mdata['train_resources'] = profile['resources']
	                    mdata['deepmd_path'] = profile['deepmd_path']
	                    continue_flag = True
	        except:
	            pass
	    pd_flag = False
	    pd_count_list =[]
	    # pd for pending job in slurm
	    # if we need to launch new machine_idxines
	    if not continue_flag:

	        #assert isinstance(mdata['train']['machine'], list)
	        #assert isinstance(mdata['train']['resources'], list)
	        #assert len(mdata['train']['machine']) == len(mdata['train']['resources'])
	        # mdata['train'] is  a list 
	        for machine_idx in range(len(mdata['train'])):
	            temp_machine = mdata['train'][machine_idx]['machine']
	            temp_resources = mdata['train'][machine_idx]['resources']
	            #assert isinstance(temp_machine, dict), "unsupported type of train machine [%d]!" %machine_idx
	            #assert isinstance(temp_resources, dict), "unsupported type of train resources [%d]!"%machine_idx
	            #assert temp_machine['machine_type'] == 'slurm', "Currently only support for Slurm!"
	            temp_ssh_sess = SSHSession(temp_machine)
	            cwd = os.getcwd()
	            temp_rjob = SlurmJob(temp_ssh_sess, cwd)

	            ## By `squeue -u user -p partition | grep PD` 
	            command = temp_rjob._make_squeue(temp_machine, temp_resources)
	            stdin, stdout, stderr = temp_rjob.ssh.exec_command(command)
	            pd_response = stdout.read().decode('utf-8').split("\n")
	            pd_count = len(pd_response)

	            temp_rjob.clean()
	            ## If there is no need to waiting for allocation
	            if pd_count ==1:
	                mdata['train_machine'] = temp_machine   
	                mdata['train_resources'] = temp_resources
	                mdata['deepmd_path'] = mdata['train'][machine_idx]['deepmd_path']
	                ## No need to wait
	                pd_flag = True
	                break
	            else:
	                pd_count_list.append(pd_count)
	        if not pd_flag:
	        	## All machines need waiting, then compare waiting jobs
	        	## Select a machine which has fewest waiting jobs
	            min_machine_idx = np.argsort(pd_count_list)[0]
	            mdata['train_machine'] = mdata['train'][min_machine_idx]['machine']
	            mdata['train_resources'] = mdata['train'][min_machine_idx]['resources']
	            mdata['deepmd_path'] = mdata['train'][min_machine_idx]['deepmd_path']

            ## Record whihc machine is selected
	        with open("record.machine","w") as _outfile:
	            profile = {}
	            profile['purporse'] = 'train'
	            profile['machine'] = mdata['train_machine']
	            profile['resources'] = mdata['train_resources']
	            profile['deepmd_path'] = mdata['deepmd_path']
	            json.dump(profile, _outfile, indent = 4)
	return mdata

def decide_model_devi_machine(mdata):
	if 'model_devi' in mdata:
	    continue_flag = False
	    if 'record.machine' in os.listdir():
	        try:
	            with open('record.machine', 'r') as _infile:
	                profile = json.load(_infile)
	                if profile['purpose'] == 'model_devi':
	                    mdata['model_devi_machine'] = profile['machine']
	                    mdata['model_devi_resources'] = profile['resources']
	                    mdata['lmp_command'] = profile['command']
	                    mdata['model_devi_group_size'] = profile['group_size']
	                    continue_flag = True
	        except:
	            pass
	    pd_count_list =[]
	    pd_flag = False
	    if not continue_flag:

	        #assert isinstance(mdata['model_devi']['machine'], list)
	        #ssert isinstance(mdata['model_devi']['resources'], list)
	        #assert len(mdata['model_devi']['machine']) == len(mdata['model_devi']['resources'])
	    
	        for machine_idx in range(len(mdata['model_devi'])):
	            temp_machine = mdata['model_devi'][machine_idx]['machine']
	            temp_resources = mdata['model_devi'][machine_idx]['resources']
	            #assert isinstance(temp_machine, dict), "unsupported type of model_devi machine [%d]!" %machine_idx
	            #assert isinstance(temp_resources, dict), "unsupported type of model_devi resources [%d]!"%machine_idx
	            #assert temp_machine['machine_type'] == 'slurm', "Currently only support for Slurm!"
	            temp_ssh_sess = SSHSession(temp_machine)
	            cwd = os.getcwd()
	            temp_rjob = SlurmJob(temp_ssh_sess, cwd)
	            command = temp_rjob._make_squeue(temp_machine, temp_resources)
	            stdin, stdout, stderr = temp_rjob.ssh.exec_command(command)
	            pd_response = stdout.read().decode('utf-8').split("\n")
	            pd_count = len(pd_response)
	            temp_rjob.clean()
	            if pd_count ==0:
	                mdata['model_devi_machine'] = temp_machine   
	                mdata['model_devi_resources'] = temp_resources
	                mdata['lmp_command'] = mdata['model_devi'][machine_idx]['command']
	                mdata['model_devi_group_size'] =  mdata['model_devi'][machine_idx]['group_size']
	                pd_flag = True
	                break
	            else:
	                pd_count_list.append(pd_count)
	        if not pd_flag:
	            min_machine_idx = np.argsort(pd_count_list)[0]
	            mdata['model_devi_machine'] = mdata['model_devi'][min_machine_idx]['machine']
	            mdata['model_devi_resources'] = mdata['model_devi'][min_machine_idx]['resources']
	            mdata['lmp_command'] = mdata['model_devi'][min_machine_idx]['command']
	            mdata['model_devi_group_size'] =  mdata['model_devi'][min_machine_idx]['group_size']
	        with open("record.machine","w") as _outfile:
	            profile = {}
	            profile['purporse'] = 'model_devi'
	            profile['machine'] = mdata['model_devi_machine']
	            profile['resources'] = mdata['model_devi_resources']
	            profile['group_size'] = mdata['model_devi_group_size']
	            profile['command'] = mdata['lmp_command']

	            json.dump(profile, _outfile, indent = 4)
	return mdata
def decide_fp_machine(mdata):

	if 'fp' in mdata:
	    #ssert isinstance(mdata['fp']['machine'], list)
	    #assert isinstance(mdata['fp']['resources'], list)
	    #assert len(mdata['fp']['machine']) == len(mdata['fp']['resources'])
	    continue_flag = False
	    ## decide whether to use an existing machine
	    if 'record.machine' in os.listdir():
	        try:
	            with open('record.machine', 'r') as _infile:
	                profile = json.load(_infile)
	                if profile['purpose'] == 'fp':
	                    mdata['fp_machine'] = profile['machine']
	                    mdata['fp_resources'] = profile['resources']
	                    mdata['fp_command'] = profile['command']
	                    mdata['fp_group_size'] = profile['group_size']
	                    #mdata['deepmd_path'] = profile['deepmd_path']
	                    continue_flag = True
	        except:
	            pass
	    pd_count_list =[]
	    pd_flag = False
	    for machine_idx in range(len(mdata['fp'])):
	        temp_machine = mdata['fp'][machine_idx]['machine']
	        temp_resources = mdata['fp'][machine_idx]['resources']
	        #assert isinstance(temp_machine, dict), "unsupported type of fp machine [%d]!" %machine_idx
	        #assert isinstance(temp_resources, dict), "unsupported type of fp resources [%d]!"%machine_idx
	        #assert temp_machine['machine_type'] == 'slurm', "Currently only support for Slurm!"
	        temp_ssh_sess = SSHSession(temp_machine)
	        cwd = os.getcwd()
	        temp_rjob = SlurmJob(temp_ssh_sess, cwd)
	        command = temp_rjob._make_squeue(temp_machine, temp_resources)
	        stdin, stdout, stderr = temp_rjob.ssh.exec_command(command)
	        pd_response = stdout.read().decode('utf-8').split("\n")
	        pd_count = len(pd_response)
	        temp_rjob.clean()
	        if pd_count ==0:
	            mdata['fp_machine'] = temp_machine   
	            mdata['fp_resources'] = temp_resources
	            mdata['fp_command'] = mdata['fp'][machine_idx]['command']
	            mdata['fp_group_size'] =  mdata['fp'][machine_idx]['group_size']
	            pd_flag = True
	            break
	        else:
	            pd_count_list.append(pd_count)
	    if not pd_flag:
	        min_machine_idx = np.argsort(pd_count_list)[0]
	        mdata['fp_machine'] = mdata['fp'][min_machine_idx]['machine']
	        mdata['fp_resources'] = mdata['fp'][min_machine_idx]['resources']
	        mdata['fp_command'] = mdata['fp'][min_machine_idx]['command']
	        mdata['fp_group_size'] =  mdata['fp'][min_machine_idx]['group_size']

	    with open("record.machine","w") as _outfile:
	            profile = {}
	            profile['purporse'] = 'fp'
	            profile['machine'] = mdata['fp_machine']
	            profile['resources'] = mdata['fp_resources']
	            profile['group_size'] = mdata['fp_group_size']
	            profile['command'] = mdata['fp_command']
	            json.dump(profile, _outfile, indent = 4)
	return mdata


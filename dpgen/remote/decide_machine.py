#!/usr/bin/env python
# coding: utf-8

from dpgen.dispatcher.SSHContext import SSHSession
from dpgen.dispatcher.SSHContext import SSHContext
from dpgen.dispatcher.Slurm import Slurm
from dpgen import dlog
import os
import json
import numpy as np


def decide_train_machine(mdata):
	if 'train' in mdata:
	    continue_flag = False
	    if 'record.machine' in os.listdir():
	        try:
	            with open('record.machine', 'r') as _infile:
	                profile = json.load(_infile)
	                if profile['purpose'] == 'train':
	                    mdata['train_machine'] = profile['machine']
	                    mdata['train_resources'] = profile['resources']
	                    if 'deepmd_path' in profile:
	                    	mdata['deepmd_path'] = profile['deepmd_path']
	                    elif 'python_path' in profile:
	                    	mdata['python_path'] = profile['python_path']
	                    if "group_size" in profile:
	                    	mdata["train_group_size"] = profile["group_size"]
	                    if 'deepmd_version' in profile:
	                    	mdata["deepmd_version"] = profile['deepmd_version']
	                    continue_flag = True
	        except:
	            pass
	    if ("hostname" not in mdata["train"][0]["machine"]) or (len(mdata["train"]) == 1):
	    	mdata["train_machine"] = mdata["train"][0]["machine"]
	    	mdata["train_resources"] = mdata["train"][0]["resources"]
	    	if 'deepmd_path' in mdata["train"][0]:
	    		mdata["deepmd_path"] = mdata["train"][0]["deepmd_path"]
	    	elif 'python_path' in mdata["train"][0]:
	    		mdata["python_path"] = mdata["train"][0]["python_path"]
	    	if "group_size" in mdata["train"][0]:
	    		mdata["train_group_size"] = mdata["train"][0]["group_size"]
	    	if 'deepmd_version' in mdata["train"][0]:
	    		mdata["deepmd_version"] = mdata["train"][0]["deepmd_version"]
	    	if 'command' in mdata["train"][0]:
	    		mdata["train_command"] = mdata["train"][0]["command"]
	    	continue_flag = True
		
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
	            temp_ssh_sess = SSHSession(temp_machine)
	            cwd = os.getcwd()
	            temp_context = SSHContext(cwd, temp_ssh_sess)
	            temp_batch = Slurm(temp_context)
	            command = temp_batch._make_squeue(temp_machine, temp_resources)
	            ret, stdin, stdout, stderr = temp_batch.context.block_call(command)
	            pd_response = stdout.read().decode('utf-8').split("\n")
	            pd_count = len(pd_response)
	            temp_context.clean()
	            ## If there is no need to waiting for allocation
	            if pd_count ==1:
	                mdata['train_machine'] = temp_machine   
	                mdata['train_resources'] = temp_resources
	                if 'deepmd_path' in mdata['train'][machine_idx]:
	                	mdata['deepmd_path'] = mdata['train'][machine_idx]['deepmd_path']
	                elif 'python_path' in mdata['train'][machine_idx]:
	                	mdata['python_path'] = mdata['train'][machine_idx]['python_path']
	                if 'group_size' in mdata['train'][machine_idx]:
	                	mdata['train_group_size'] = mdata['train'][machine_idx]['group_size']
	                if 'deepmd_version' in mdata['train'][machine_idx]:
	                	mdata['deepmd_version'] = mdata['train'][machine_idx]['deepmd_version']

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
	            if 'deepmd_path' in mdata['train'][min_machine_idx]:
	            	mdata['deepmd_path'] = mdata['train'][min_machine_idx]['deepmd_path']
	            elif 'python_path' in mdata['train'][min_machine_idx]:
	            	mdata['python_path'] = mdata['train'][min_machine_idx]['python_path']
	            if "group_size" in mdata['train'][min_machine_idx]:
	            	mdata["train_group_size"] = mdata['train'][min_machine_idx]["group_size"]
	            if 'deepmd_version' in mdata['train'][min_machine_idx]:
	            	mdata['deepmd_version'] = mdata['train'][min_machine_idx]["deepmd_version"]

            ## Record which machine is selected
	        with open("record.machine","w") as _outfile:
	            profile = {}
	            profile['purpose'] = 'train'
	            profile['machine'] = mdata['train_machine']
	            profile['resources'] = mdata['train_resources']
	            if 'deepmd_path' in mdata:
	            	profile['deepmd_path'] = mdata['deepmd_path']
	            elif 'python_path' in mdata:
	            	profile['python_path'] = mdata['python_path']
	            if "train_group_size" in mdata:
	            	profile["group_size"] = mdata["train_group_size"]
	            if 'deepmd_version' in mdata:
	            	profile['deepmd_version'] = mdata['deepmd_version']

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
	    if ("hostname" not in mdata["model_devi"][0]["machine"]) or (len(mdata["model_devi"]) == 1):
	    	mdata["model_devi_machine"] = mdata["model_devi"][0]["machine"]
	    	mdata["model_devi_resources"] = mdata["model_devi"][0]["resources"]
	    	mdata["lmp_command"] = mdata["model_devi"][0]["command"]
	    	#if "group_size" in mdata["train"][0]:
	    	mdata["model_devi_group_size"] = mdata["model_devi"][0]["group_size"]
	    	continue_flag = True

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
	            temp_context = SSHContext(cwd, temp_ssh_sess)
	            temp_batch = Slurm(temp_context)
	            command = temp_batch._make_squeue(temp_machine, temp_resources)
	            ret, stdin, stdout, stderr = temp_batch.context.block_call(command)
	            pd_response = stdout.read().decode('utf-8').split("\n")
	            pd_count = len(pd_response)
	            temp_context.clean()
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
	            profile['purpose'] = 'model_devi'
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
	    if ("hostname" not in mdata["fp"][0]["machine"]) or (len(mdata["fp"]) == 1):
	    	mdata["fp_machine"] = mdata["fp"][0]["machine"]
	    	mdata["fp_resources"] = mdata["fp"][0]["resources"]
	    	mdata["fp_command"] = mdata["fp"][0]["command"]
	    	#if "group_size" in mdata["train"][0]:
	    	mdata["fp_group_size"] = mdata["fp"][0]["group_size"]
	    	continue_flag = True

	   
	    pd_count_list =[]
	    pd_flag = False
	    if not continue_flag:
		    for machine_idx in range(len(mdata['fp'])):
		        temp_machine = mdata['fp'][machine_idx]['machine']
		        temp_resources = mdata['fp'][machine_idx]['resources']
		        temp_ssh_sess = SSHSession(temp_machine)
		        cwd = os.getcwd()
		        temp_context = SSHContext(cwd, temp_ssh_sess)
		        temp_batch = Slurm(temp_context)
		        command = temp_batch._make_squeue(temp_machine, temp_resources)
		        ret, stdin, stdout, stderr = temp_batch.context.block_call(command)
		        pd_response = stdout.read().decode('utf-8').split("\n")
		        pd_count = len(pd_response)
		        temp_context.clean()
		        #dlog.info(temp_machine["username"] + " " + temp_machine["hostname"] +  " " + str(pd_count))
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
		            profile['purpose'] = 'fp'
		            profile['machine'] = mdata['fp_machine']
		            profile['resources'] = mdata['fp_resources']
		            profile['group_size'] = mdata['fp_group_size']
		            profile['command'] = mdata['fp_command']
		            json.dump(profile, _outfile, indent = 4)
	return mdata


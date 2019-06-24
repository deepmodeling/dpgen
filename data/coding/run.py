#!/usr/bin/env python3
# main program for preparing data
import os,json,shutil,re,glob,argparse
import numpy as np
import subprocess as sp
import tools.hcp as hcp
import tools.fcc as fcc
import tools.bcc as bcc
import tools.diamond as diamond
import tools.sc as sc
import os
import sys
import argparse
import glob
import json
import random
import logging
import warnings
import shutil
import time
import dpdata
import numpy as np
import subprocess as sp
from lib.vasp import write_incar_dict
from lib.vasp import make_vasp_incar_user_dict
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

import requests
from hashlib import sha1


### utils

import os, re, shutil, logging


task_format = "%02d"
log_iter_head = " task " + task_format + ": "

def make_iter_name (task_index) :
    return "task." + (task_format % task_index)

def create_path (path) :
    path += '/'
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)        
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
            counter += 1
    os.makedirs (path)

def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

def copy_file_list (file_list, from_path, to_path) :
    for jj in file_list : 
        if os.path.isfile(os.path.join(from_path, jj)) :
            shutil.copy (os.path.join(from_path, jj), to_path)
        elif os.path.isdir(os.path.join(from_path, jj)) :
            shutil.copytree (os.path.join(from_path, jj), os.path.join(to_path, jj))

def cmd_append_log (cmd,
                    log_file) :
    ret = cmd
    ret = ret + " 1> " + log_file
    ret = ret + " 2> " + log_file
    return ret

def log_iter (task, ii, ) :
    logging.info ((log_iter_head + "%s") % (ii, jj, task))

def repeat_to_length(string_to_expand, length):
    ret = ""
    for ii in range (length) : 
        ret += string_to_expand
    return ret

def log_task (message) :
    header = repeat_to_length (" ", len(log_iter_head % (0, 0)))
    logging.info (header + message)

def record_iter (record, ii, jj) :
    with open (record, "a") as frec :
        frec.write ("%d %d\n" % (ii, jj)) 


########
#### submit
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
        rjob.upload('.',  forward_common_files)
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
####
def create_path (path) :
    path += '/'
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)        
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
            counter += 1
    os.makedirs (path)
    return path

def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

"""
0, make unit cell
1, copy
2, place element
3, relax
4, perturb
"""
global_dirname_02 = '00.place_ele'
global_dirname_03 = '01.scale_pert'
global_dirname_04 = '02.md'

def out_dir_name(jdata) :
    cell_type = jdata['cell_type']
    elements = jdata['elements']
    super_cell = jdata['super_cell']    
    from_poscar = False
    if 'from_poscar' in jdata :
        from_poscar = jdata['from_poscar']
        from_poscar_path = jdata['from_poscar_path']

    if from_poscar:
        poscar_name = os.path.basename(from_poscar_path)
        cell_str = "%02d" % (super_cell[0])
        for ii in range(1,len(super_cell)) :
            cell_str = cell_str + ("x%02d" % super_cell[ii])
        return poscar_name + '.' + cell_str
    else :
        ele_str = ""
        for ii in elements:
            ele_str = ele_str + ii.lower()
        cell_str = "%02d" % (super_cell[0])
        for ii in range(1,len(super_cell)) :
            cell_str = cell_str + ("x%02d" % super_cell[ii])
        return ele_str + '.' + cell_type + '.' + cell_str

def class_cell_type(jdata) :
    ct = jdata['cell_type']
    if ct == "hcp" :
        cell_type = hcp
    elif ct == "fcc" :
        cell_type = fcc
    elif ct == "diamond" :
        cell_type = diamond
    elif ct == "sc" :
        cell_type = sc
    elif ct == "bcc" :
        cell_type = bcc
    else :
        raise RuntimeError("unknow cell type %s" % ct)
    return cell_type

def poscar_ele(poscar_in, poscar_out, eles, natoms) :
    ele_line = ""
    natom_line = ""
    for ii in eles :
        ele_line += str(ii) + " "
    for ii in natoms :
        natom_line += str(ii) + " "
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
        lines[5] = ele_line + "\n"
        lines[6] = natom_line + "\n"
    with open(poscar_out, 'w') as fout :
        fout.write("".join(lines))

def poscar_natoms(lines) :
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    return numb_atoms

def poscar_shuffle(poscar_in, poscar_out) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    numb_atoms = poscar_natoms(lines)
    idx = np.arange(8, 8+numb_atoms)
    np.random.shuffle(idx)
    out_lines = lines[0:8]
    for ii in range(numb_atoms) :
        out_lines.append(lines[idx[ii]])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(out_lines))

def poscar_scale_direct (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = poscar_natoms(lines)
    pscale = float(lines[1])
    pscale = pscale * scale
    lines[1] = str(pscale) + "\n"
    return lines

def poscar_scale_cartesian (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = poscar_natoms(lines)
    # scale box
    for ii in range(2,5) :
        boxl = lines[ii].split()
        boxv = [float(ii) for ii in boxl]
        boxv = np.array(boxv) * scale
        lines[ii] = "%.16e %.16e %.16e\n" % (boxv[0], boxv[1], boxv[2])
    # scale coord
    for ii in range(8, 8+numb_atoms) :
        cl = lines[ii].split()
        cv = [float(ii) for ii in cl]
        cv = np.array(cv) * scale
        lines[ii] = "%.16e %.16e %.16e\n" % (cv[0], cv[1], cv[2])
    return lines    

def poscar_scale (poscar_in, poscar_out, scale) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    if 'D' == lines[7][0] or 'd' == lines[7][0]: 
        lines = poscar_scale_direct(lines, scale)
    elif 'C' == lines[7][0] or 'c' == lines[7][0] :
        lines = poscar_scale_cartesian(lines, scale)
    else :
        raise RuntimeError("Unknow poscar style at line 7: %s" % lines[7])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(lines))

def make_unit_cell (jdata) :
    latt = jdata['latt']
    out_dir = jdata['out_dir']
    path_uc = os.path.join(out_dir, global_dirname_02)
    cell_type = class_cell_type(jdata)

    cwd = os.getcwd()    
    # for ii in scale :
    # path_work = create_path(os.path.join(path_uc, '%.3f' % ii))
    path_work = create_path(path_uc)    
    os.chdir(path_work)
    with open('POSCAR.unit', 'w') as fp:
        fp.write (cell_type.poscar_unit(latt))
    os.chdir(cwd)        

def make_super_cell (jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    path_uc = os.path.join(out_dir, global_dirname_02)
    path_sc = os.path.join(out_dir, global_dirname_02)
    assert(os.path.isdir(path_uc)), "path %s should exists" % path_uc
    assert(os.path.isdir(path_sc)), "path %s should exists" % path_sc

    # for ii in scale :
    from_path = path_uc
    from_file = os.path.join(from_path, 'POSCAR.unit')
    to_path = path_sc
    to_file = os.path.join(to_path, 'POSCAR')
    cmd = "./tools/poscar_copy.py -n %d %d %d " % (super_cell[0], super_cell[1], super_cell[2]) + \
          from_file + " " + \
          to_file
    sp.check_call(cmd, shell = True)

def make_super_cell_poscar(jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    path_sc = os.path.join(out_dir, global_dirname_02)    
    create_path(path_sc)
    from_poscar_path = jdata['from_poscar_path']
    assert(os.path.isfile(from_poscar_path)), "file %s should exists" % from_poscar_path
    
    from_file = os.path.join(path_sc, 'POSCAR.copied')
    shutil.copy2(from_poscar_path, from_file)
    to_file = os.path.join(path_sc, 'POSCAR')
    cmd = "./tools/poscar_copy.py -n %d %d %d " % (super_cell[0], super_cell[1], super_cell[2]) + \
          from_file + " " + \
          to_file
    sp.check_call(cmd, shell = True)    

    # make system dir (copy)
    lines = open(to_file, 'r').read().split('\n')
    natoms_str = lines[6]
    natoms_list = [int(ii) for ii in natoms_str.split()]
    print(natoms_list)
    comb_name = "sys-"
    for idx,ii in enumerate(natoms_list) :
        comb_name += "%04d" % ii
        if idx != len(natoms_list)-1 :
            comb_name += "-"
    path_work = os.path.join(path_sc, comb_name)
    create_path(path_work)
    cwd = os.getcwd()
    to_file = os.path.abspath(to_file)
    os.chdir(path_work)
    os.symlink(os.path.relpath(to_file), 'POSCAR')
    os.chdir(cwd)

def make_combines (dim, natoms) :
    if dim == 1 :
        return [[natoms]]
    else :
        res = []
        for ii in range(natoms+1) :
            rest = natoms - ii
            tmp_combines = make_combines(dim-1, rest)
            for jj in tmp_combines :
                jj.append(ii)
            if len(res) == 0 :
                res = tmp_combines
            else :
                res += tmp_combines
        return res

def place_element (jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    cell_type = class_cell_type(jdata)
    natoms = np.cumprod(super_cell)[-1] * cell_type.numb_atoms()
    elements = jdata['elements']
    path_sc = os.path.join(out_dir, global_dirname_02)
    path_pe = os.path.join(out_dir, global_dirname_02)    
    combines = np.array(make_combines(len(elements), natoms), dtype = int)
    
    assert(os.path.isdir(path_pe))
    cwd = os.getcwd()
    for ii in combines :
        if any(ii == 0) :
            continue
        comb_name = "sys-"
        for idx,jj in enumerate(ii) :            
            comb_name += "%04d" % jj
            if idx != len(ii)-1 :
                comb_name += "-"
        path_pos_in = path_sc
        path_work = os.path.join(path_pe, comb_name)
        create_path(path_work)
        pos_in = os.path.join(path_pos_in, 'POSCAR')
        pos_out = os.path.join(path_work, 'POSCAR')
        poscar_ele(pos_in, pos_out, elements, ii)
        poscar_shuffle(pos_out, pos_out)

def make_vasp_relax (jdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    encut = jdata['encut']
    kspacing = jdata['kspacing_relax']
    kgamma = jdata['kgamma']
    ismear = 1
    if 'ismear' in jdata :
        ismear = jdata['ismear']
    sigma = 0.2
    if 'sigma' in jdata :
        sigma = jdata['sigma']
    cwd = os.getcwd()
    vasp_dir = os.path.join(cwd, 'vasp.in')

    work_dir = os.path.join(out_dir, global_dirname_02)
    assert (os.path.isdir(work_dir))
    work_dir = os.path.abspath(work_dir)
    if os.path.isfile(os.path.join(work_dir, 'INCAR' )) :
        os.remove(os.path.join(work_dir, 'INCAR' ))
    if os.path.isfile(os.path.join(work_dir, 'POTCAR')) :
        os.remove(os.path.join(work_dir, 'POTCAR'))
    shutil.copy2(os.path.join(vasp_dir, 'INCAR.rlx' ), 
                 os.path.join(work_dir, 'INCAR'))
    out_potcar = os.path.join(work_dir, 'POTCAR')
    with open(out_potcar, 'w') as outfile:
        for fname in potcars:
            with open(fname) as infile:
                outfile.write(infile.read())
    
    os.chdir(work_dir)
    
    replace('INCAR', 'ENCUT=.*', 'ENCUT=%f' % encut)
    replace('INCAR', 'ISIF=.*', 'ISIF=3')
    replace('INCAR', 'KSPACING=.*', 'KSPACING=%f' % kspacing)
    if kgamma :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=T')
    else :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=F')
    replace('INCAR', 'ISMEAR=.*', 'ISMEAR=%d' % ismear)
    replace('INCAR', 'SIGMA=.*', 'SIGMA=%f' % sigma)
    
    sys_list = glob.glob('sys-*')
    for ss in sys_list:
        os.chdir(ss)
        ln_src = os.path.relpath(os.path.join(work_dir,'INCAR'))
        os.symlink(ln_src, 'INCAR')
        ln_src = os.path.relpath(os.path.join(work_dir,'POTCAR'))
        os.symlink(ln_src, 'POTCAR')
        os.chdir(work_dir)
    os.chdir(cwd)

def make_scale(jdata):
    out_dir = jdata['out_dir']
    scale = jdata['scale']    
    skip_relax = jdata['skip_relax']    

    cwd = os.getcwd()
    init_path = os.path.join(out_dir, global_dirname_02)
    init_path = os.path.abspath(init_path)
    work_path = os.path.join(out_dir, global_dirname_03)
    os.chdir(init_path)
    init_sys = glob.glob("sys-*")
    init_sys.sort()
    os.chdir(cwd)

    create_path(work_path)
    for ii in init_sys :
        for jj in scale :
            pos_src = os.path.join(os.path.join(init_path, ii), 'CONTCAR')
            if not os.path.isfile(pos_src):
                if skip_relax :
                    pos_src = os.path.join(os.path.join(init_path, ii), 'POSCAR')
                    assert(os.path.isfile(pos_src))
                else :
                    raise RuntimeError("not file %s, vasp relaxation should be run before scale poscar")
            scale_path = os.path.join(work_path, ii)
            scale_path = os.path.join(scale_path, "scale-%.3f" % jj)
            create_path(scale_path)
            os.chdir(scale_path) 
            poscar_scale(pos_src, 'POSCAR', jj)
            os.chdir(cwd)

def pert_scaled(jdata) :
    out_dir = jdata['out_dir']
    scale = jdata['scale']    
    pert_box = jdata['pert_box']
    pert_atom = jdata['pert_atom']
    pert_numb = jdata['pert_numb']
    from_poscar = False 
    if 'from_poscar' in jdata :
        from_poscar = jdata['from_poscar']
    
    cwd = os.getcwd()
    path_sp = os.path.join(out_dir, global_dirname_03)
    assert(os.path.isdir(path_sp))
    os.chdir(path_sp)
    sys_pe = glob.glob('sys-*')
    sys_pe.sort()
    os.chdir(cwd)    

    pert_cmd = cwd
    pert_cmd = os.path.join(pert_cmd, 'tools')
    pert_cmd = os.path.join(pert_cmd, 'create_random_disturb.py')
    pert_cmd += ' -etmax %f -ofmt vasp POSCAR %d %f > /dev/null' %(pert_box, pert_numb, pert_atom)
    for ii in sys_pe :
        for jj in scale :
            path_work = path_sp
            path_work = os.path.join(path_work, ii)
            path_work = os.path.join(path_work, 'scale-%.3f' % jj)
            assert(os.path.isdir(path_work))
            os.chdir(path_work)
            sp.check_call(pert_cmd, shell = True)
            for kk in range(pert_numb) :
                pos_in = 'POSCAR%d.vasp' % (kk+1)
                dir_out = '%06d' % (kk+1)
                create_path(dir_out)
                pos_out = os.path.join(dir_out, 'POSCAR')
                if not from_poscar:
                    poscar_shuffle(pos_in, pos_out)
                else :
                    shutil.copy2(pos_in, pos_out)
                os.remove(pos_in)
            kk = -1
            pos_in = 'POSCAR'
            dir_out = '%06d' % (kk+1)
            create_path(dir_out)
            pos_out = os.path.join(dir_out, 'POSCAR')
            if not from_poscar:
                poscar_shuffle(pos_in, pos_out)
            else :
                shutil.copy2(pos_in, pos_out)
            os.chdir(cwd)

def make_vasp_md(jdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    scale = jdata['scale']    
    encut = jdata['encut']
    kspacing = jdata['kspacing_md']
    kgamma = jdata['kgamma']
    pert_numb = jdata['pert_numb']
    md_temp = jdata['md_temp']
    md_nstep = jdata['md_nstep']
    ismear = 1
    if 'ismear' in jdata :
        ismear = jdata['ismear']
    sigma = 0.2
    if 'sigma' in jdata :
        sigma = jdata['sigma']

    cwd = os.getcwd()
    vasp_dir = os.path.join(cwd, 'vasp.in')
    vasp_dir = os.path.join(cwd, vasp_dir)
    path_ps = os.path.join(out_dir, global_dirname_03)
    path_ps = os.path.abspath(path_ps)
    assert(os.path.isdir(path_ps))
    os.chdir(path_ps)
    sys_ps = glob.glob('sys-*')
    sys_ps.sort()
    os.chdir(cwd) 
    path_md = os.path.join(out_dir, global_dirname_04)
    path_md = os.path.abspath(path_md)
    create_path(path_md)
    shutil.copy2(os.path.join(vasp_dir, 'INCAR.md'), 
                 os.path.join(path_md, 'INCAR'))
    out_potcar = os.path.join(path_md, 'POTCAR')
    with open(out_potcar, 'w') as outfile:
        for fname in potcars:
            with open(fname) as infile:
                outfile.write(infile.read())
    os.chdir(path_md)
    replace('INCAR', 'ENCUT=.*', 'ENCUT=%f' % encut)
    replace('INCAR', 'ISIF=.*', 'ISIF=2')
    replace('INCAR', 'KSPACING=.*', 'KSPACING=%f' % kspacing)
    if kgamma :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=T')
    else :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=F')    
    replace('INCAR', 'NSW=.*', 'NSW=%d' % md_nstep)
    replace('INCAR', 'TEBEG=.*', 'TEBEG=%d' % md_temp)
    replace('INCAR', 'TEEND=.*', 'TEEND=%d' % md_temp)
    replace('INCAR', 'ISMEAR=.*', 'ISMEAR=%d' % ismear)
    replace('INCAR', 'SIGMA=.*', 'SIGMA=%f' % sigma)
    os.chdir(cwd)    

    for ii in sys_ps :
        for jj in scale :
            for kk in range(pert_numb) :
                path_work = path_md
                path_work = os.path.join(path_work, ii)
                path_work = os.path.join(path_work, "scale-%.3f" % jj)
                path_work = os.path.join(path_work, "%06d" % kk)
                create_path(path_work)
                os.chdir(path_work)                
                path_pos = path_ps
                path_pos = os.path.join(path_pos, ii)
                path_pos = os.path.join(path_pos, "scale-%.3f" % jj)
                path_pos = os.path.join(path_pos, "%06d" % kk)
                init_pos = os.path.join(path_pos, 'POSCAR')
                shutil.copy2 (init_pos, 'POSCAR')
                file_incar = os.path.join(path_md, 'INCAR')
                file_potcar = os.path.join(path_md, 'POTCAR')
                os.symlink(os.path.relpath(file_incar), 'INCAR')
                os.symlink(os.path.relpath(file_potcar), 'POTCAR')
                os.chdir(cwd)                

def _vasp_check_fin (ii) :
    if os.path.isfile(os.path.join(ii, 'OUTCAR')) :
        with open(os.path.join(ii, 'OUTCAR'), 'r') as fp :
            content = fp.read()
            count = content.count('Elapse')
            if count != 1 :
                return False
    else :
        return False
    return True

def _qe_check_fin(ii) :
    if os.path.isfile(os.path.join(ii, 'output')) :
        with open(os.path.join(ii, 'output'), 'r') as fp :
            content = fp.read()
            count = content.count('JOB DONE')
            if count != 1 :
                return False
    else :
        return False
    return True
    
        
def run_fp_inner (iter_index,
                  jdata,
                  mdata,
                  ssh_sess,
                  forward_files,
                  backward_files,
                  check_fin,
                  log_file = "log") :
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']
    machine_type = mdata['fp_machine']['machine_type']
    # fp_command = ("OMP_NUM_THREADS=1 mpirun -n %d " % fp_np) + fp_command
    # cpu task in parallel
    # if (('numb_gpu' not in fp_resources) or (fp_resources['numb_gpu'] == 0)) and (machine_type == 'slurm'):
    #     fp_command = 'srun ' + fp_command
    # if (('numb_gpu' not in fp_resources) or (fp_resources['numb_gpu'] == 0)) and (machine_type == 'pbs'):
    #     fp_command = 'mpirun  ' + fp_command
    fp_command = cmd_append_log(fp_command, log_file)

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0 :
        return

    fp_run_tasks = []
    for ii in fp_tasks :
        if not check_fin(ii) :
            fp_run_tasks.append(ii)

    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]

    if ssh_sess == None and machine_type == 'ucloud':
        _ucloud_submit_jobs(mdata['fp_machine'],
                            mdata['fp_resources'],
                            fp_command,
                            work_path,
                            run_tasks,
                            fp_group_size,
                            [],
                            forward_files,
                            backward_files)
    elif machine_type == 'slurm' :        
        _group_slurm_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files)
    elif machine_type == 'pbs' :        
        _group_slurm_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files,
                          remote_job = PBSJob)
    elif machine_type == 'local' :        
        _group_local_jobs(ssh_sess,
                           fp_resources,
                           fp_command,
                           work_path,
                           run_tasks,
                           fp_group_size,
                           [],
                           forward_files,
                           backward_files)
    else :
        raise RuntimeError("unknown machine type")

#    exec_hosts_batch(exec_machine, fp_command, fp_np, fp_run_tasks, verbose = True, mpi = False, gpu=True)
    # exec_batch_group(fp_command,
    #                  fp_nn, fp_nppn, fp_ngpu,
    #                  fp_run_tasks,
    #                  group_size = fp_group_size,
    #                  time_limit = fp_tlimit,
    #                  modules = fp_modules,
    #                  sources = fp_sources)
        
def run_fp (iter_index,
            jdata,
            mdata,
            ssh_sess) :
    fp_pp_files = jdata['fp_pp_files']
    forward_files = ['POSCAR', 'INCAR'] + fp_pp_files
    backward_files = ['OUTCAR']
    run_fp_inner(iter_index, jdata, mdata, ssh_sess, forward_files, backward_files, _vasp_check_fin) 


def coll_vasp_md(jdata) :
    out_dir = jdata['out_dir']
    md_nstep = jdata['md_nstep']
    scale = jdata['scale']    
    pert_numb = jdata['pert_numb']
    deepgen_templ = jdata['deepgen_templ']
    coll_ndata = jdata['coll_ndata']
    raw_files = ['box.raw', 'coord.raw', 'energy.raw', 'force.raw', 'virial.raw']

    deepgen_templ = os.path.abspath(deepgen_templ)
    cmd_cvt = os.path.join(deepgen_templ, 'tools.vasp')
    cmd_cvt = os.path.join(cmd_cvt, 'cessp2force_lin.py')
    cmd_2raw = os.path.join(deepgen_templ, 'tools.vasp')
    cmd_2raw = os.path.join(cmd_2raw, 'convert2raw.py')
    cmd_shfl = os.path.join(deepgen_templ, 'tools.raw')
    cmd_shfl = os.path.join(cmd_shfl, 'shuffle_raw.py')
    cmd_2set = os.path.join(deepgen_templ, 'tools.raw')
    cmd_2set = os.path.join(cmd_2set, 'raw_to_set.sh')

    cwd = os.getcwd()
    path_md = os.path.join(out_dir, global_dirname_04)
    path_md = os.path.abspath(path_md)
    assert(os.path.isdir(path_md)), "md path should exists"
    os.chdir(path_md)
    sys_md = glob.glob('sys-*')
    sys_md.sort()

    for ii in sys_md :
        os.chdir(ii)
        # convert outcars
        valid_outcars = []
        for jj in scale :
            for kk in range(pert_numb) :
                path_work = os.path.join("scale-%.3f" % jj, "%06d" % kk)
                outcar = os.path.join(path_work, 'OUTCAR')
                if os.path.isfile(outcar) :
                    with open(outcar, 'r') as fin:
                        nforce = fin.read().count('TOTAL-FORCE')
                    if nforce == md_nstep :
                        valid_outcars.append(outcar)
        arg_cvt = " "
        if len(valid_outcars) == 0:
            raise RuntimeError("MD dir: %s: find no valid outcar in sys %s, "
                               "check if your vasp md simulation is correctly done" 
                               % (path_md, ii)) 
        for ii in valid_outcars :
            arg_cvt += re.sub('OUTCAR', '', ii) + " "
        tmp_cmd_cvt = cmd_cvt + arg_cvt
        tmp_cmd_cvt += ' 1> cvt.log 2> cvt.log '
        sp.check_call(tmp_cmd_cvt, shell = True)
        # create deepmd data
        if os.path.isdir('deepmd') :
            shutil.rmtree('deepmd')
        os.mkdir('deepmd')
        os.chdir('deepmd')
        os.mkdir('orig')
        os.mkdir('shuffled')
        os.chdir('orig')
        sp.check_call(cmd_2raw + ' ../../test.configs', shell = True)
        os.chdir('..')
        sp.check_call(cmd_shfl + ' orig shuffled ', shell = True)
        for ii in raw_files:
            sp.check_call('head -n %d shuffled/%s > %s' % (coll_ndata, ii, ii), shell = True)
        shutil.copy2(os.path.join('orig', 'type.raw'), 'type.raw')
        print(cmd_2set + (' %d '%coll_ndata))
        sp.check_call(cmd_2set + (' %d '%coll_ndata), shell = True)
        os.chdir(path_md)
    os.chdir(cwd)
 
def run_iter(PARAM, MACHINE):
	with open (PARAM, 'r') as fp :
        jdata = json.load (fp)
    with open (MACHINE, 'r') as fm:
    	mdata = json.load(fm)
    out_dir = out_dir_name(jdata)
    jdata['out_dir'] = out_dir
    from_poscar = False 
    if 'from_poscar' in jdata :
        from_poscar = jdata['from_poscar']
    print ("# working dir %s" % out_dir)
    
    fp_machine = mdata['fp_machine']    
    if ('machine_type' in fp_machine) and  \
       (fp_machine['machine_type'] == 'ucloud'):
        fp_ssh_sess = None
    else :
        fp_ssh_sess = SSHSession(fp_machine)

    max_stage = jdata["max_stage"]
    min_stage = 1
    record = "record.dpgen"
    if os.path.isfile (record) :
        with open (record) as frec :
            for line in frec : 
                min_stage = int(frec.readline()[:-1])
        logging.info ("continue from task %03d" % (min_stage+1))
    for stage in range(min_stage+1, max_stage+1)
	    if stage == 1 :
	        create_path(out_dir)
	        if from_poscar :
	            make_super_cell_poscar(jdata)
	        else :
	            make_unit_cell(jdata)
	            make_super_cell(jdata)
	            place_element(jdata)
	        make_vasp_relax(jdata)
	      	###run_vasp_relax()
	    elif stage == 2 :
	        make_scale(jdata)
	        pert_scaled(jdata)
	    elif stage == 3 :
	        make_vasp_md(jdata)
	    elif stage == 4 :
	        coll_vasp_md(jdata)
	    else :
	        raise RuntimeError("unknow stage %d" % stage)

def _main() :
    parser = argparse.ArgumentParser(
        description="gen init confs")
    parser.add_argument('PARAM', type=str, 
                        help="parameter file, json format")
    '''
    parser.add_argument('STAGE', type=int,
                        help="the stage of init, can be 1, 2, 3 or 4. "
                        "1: Setup vasp jobs for relaxation. "
                        "2: Collect vasp relaxed confs (if relax is not skiped). Perturb system. "
                        "3: Setup vasp jobs for MD of perturbed system. "
                        "4: Collect vasp md confs, make deepmd data. "
    )
    '''
    parser.add_argument('MACHINE', type=str, help = "machine filem, json format")
    args = parser.parse_args()
    run_iter(args.PARAM, args.MACHINE)
    logging.info ("finished!")
    
if __name__ == "__main__":
    _main()

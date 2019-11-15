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
import dpgen.auto_test.lib.util as util
import dpgen.auto_test.lib.vasp as vasp
import dpgen.auto_test.lib.lammps as lammps
import random
import logging
import warnings
import shutil
import time
import numpy as np
import subprocess as sp
from dpgen.auto_test.lib.utils import make_iter_name
from dpgen.auto_test.lib.utils import create_path
from dpgen.auto_test.lib.utils import copy_file_list
from dpgen.auto_test.lib.utils import replace
from dpgen.dispatcher.Dispatcher import make_dispatcher
from dpgen.auto_test.lib.utils import log_iter
from dpgen.auto_test.lib.utils import record_iter
from dpgen.auto_test.lib.utils import log_iter
from dpgen.auto_test.lib.pwscf import make_pwscf_input
from dpgen.auto_test.lib.siesta import make_siesta_input
from dpgen.remote.RemoteJob import SSHSession, JobStatus, SlurmJob, PBSJob, CloudMachineJob
from dpgen.remote.decide_machine import decide_fp_machine, decide_model_devi_machine
from dpgen.remote.group_jobs import *
from dpgen.auto_test import gen_00_equi,cmpt_00_equi
from dpgen.auto_test import gen_01_eos,cmpt_01_eos
from dpgen.auto_test import gen_02_elastic,cmpt_02_elastic
from dpgen.auto_test import gen_03_vacancy,cmpt_03_vacancy
from dpgen.auto_test import gen_04_interstitial,cmpt_04_interstitial
from dpgen.auto_test import gen_05_surf,cmpt_05_surf
#from dpgen.auto_test import gen_06_phonon,cmpt_06_phonon
from dpgen.auto_test import gen_confs
import requests
from hashlib import sha1

lammps_task_type=['deepmd','meam','eam']

def gen_equi(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cwd=os.getcwd()
    #vasp
    if task_type=="vasp":
        gen_00_equi.make_vasp(jdata, conf_dir)
    #lammps
    elif task_type in lammps_task_type:
        gen_00_equi.make_lammps (jdata, conf_dir,task_type)
    else :
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    os.chdir(cwd)

def run_equi(task_type,jdata,mdata):
        #rmprint("This module has been run !")

    work_path=util.make_work_path(jdata,'00.equi',False,False,False)
    all_task = glob.glob(os.path.join(work_path,'.'))

    #vasp
    if task_type=="vasp":
        mdata=decide_fp_machine(mdata)

        forward_files = ['INCAR', 'POTCAR']
        backward_files = ['OUTCAR','CONTCAR','OSZICAR']
        common_files=['POSCAR']

    #lammps
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)

        forward_files = ['conf.lmp', 'lammps.in']
        backward_files = ['dump.relax','log.lammps', 'model_devi.log']

        fp_params = jdata['lammps_params']
        model_dir = fp_params['model_dir']
        model_dir = os.path.abspath(model_dir)
        model_name =fp_params['model_name']
        if not model_name :
            models = glob.glob(os.path.join(model_dir, '*pb'))
            model_name = [os.path.basename(ii) for ii in models]
        else:
            models = [os.path.join(model_dir,ii) for ii in model_name]
        common_files = model_name

        if len(model_name)>1 and task_type == 'deepmd':
            backward_files = backward_files + ['model_devi.out']

    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    run_tasks = util.collect_task(all_task,task_type)

    machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)
    disp = make_dispatcher(machine)
    disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  common_files,
                  forward_files,
                  backward_files,
                  outlog='autotest.out',
                  errlog='autotest.err')


def cmpt_equi(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cmpt_shift=jdata['alloy_shift']
    #vasp
    if task_type=="vasp":
        n, e, v, s = cmpt_00_equi.comput_vasp_nev(jdata, conf_dir,cmpt_shift)
    #lammps
    elif task_type in lammps_task_type:
        n, e, v, s = cmpt_00_equi.comput_lmp_nev(conf_dir, task_type,cmpt_shift)
    else :
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    if cmpt_shift:
        print('conf_dir:\t EpA(eV)  VpA(A^3)  ener_shift(eV)')
        print("%s\t %8.4f  %7.3f %8.4f" % (conf_dir, e, v, s))
    else:
        print('conf_dir:\t EpA(eV)  VpA(A^3)')
        print("%s\t %8.4f  %7.3f " % (conf_dir, e, v))

def gen_eos(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    fix_shape = True
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_01_eos.make_vasp(jdata, conf_dir)
    #lammps
    elif task_type in lammps_task_type:
        if fix_shape :
            gen_01_eos.make_lammps_fixv(jdata, conf_dir,task_type)
        else :
            gen_01_eos.make_lammps(jdata, conf_dir,task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_eos(task_type,jdata,mdata):
    work_path=util.make_work_path(jdata,'01.eos',False,False,False)

    all_task = glob.glob(os.path.join(work_path, "vol-*"))
    all_task.sort()

    #vasp
    if task_type=="vasp":
        mdata=decide_fp_machine(mdata)

        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR','OSZICAR']
        common_files=['INCAR','POTCAR']

    #lammps
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)

        fp_params = jdata['lammps_params']
        model_dir = fp_params['model_dir']
        model_dir = os.path.abspath(model_dir)
        model_name =fp_params['model_name']
        if not model_name :
            models = glob.glob(os.path.join(model_dir, '*pb'))
            model_name = [os.path.basename(ii) for ii in models]
        else:
            models = [os.path.join(model_dir,ii) for ii in model_name]
        forward_files = ['conf.lmp', 'lammps.in']+model_name
        backward_files = ['log.lammps', 'model_devi.log']
        common_files=['lammps.in']+model_name

        if len(model_name)>1 and task_type == 'deepmd':
            backward_files = backward_files + ['model_devi.out']

    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    run_tasks = util.collect_task(all_task,task_type)

    machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)
    disp = make_dispatcher(machine)
    disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  common_files,
                  forward_files,
                  backward_files,
                  outlog='autotest.out',
                  errlog='autotest.err')

def cmpt_eos(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    #vasp
    if task_type == "vasp":
        cmpt_01_eos.comput_vasp_eos(jdata, conf_dir)
    #lammps
    elif task_type in lammps_task_type:
        cmpt_01_eos.comput_lmp_eos(jdata, conf_dir, task_type)
    else :
        raise RuntimeError("unknow task ", task_type)

def gen_elastic(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_02_elastic.make_vasp(jdata, conf_dir)
    #lammps
    elif task_type in lammps_task_type:
        gen_02_elastic.make_lammps (jdata, conf_dir,task_type)
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)
    os.chdir(cwd)

def run_elastic(task_type,jdata,mdata):
    work_path=util.make_work_path(jdata,'02.elastic',False,False,False)

    all_task = glob.glob(os.path.join(work_path, "dfm-*"))
    all_task.sort()

    #vasp
    if task_type == "vasp":
        mdata=decide_fp_machine(mdata)

        forward_files = ['INCAR', 'POSCAR','POTCAR','KPOINTS']
        backward_files = ['OUTCAR','CONTCAR','OSZICAR']
        common_files=['INCAR','POTCAR','KPOINTS']

    #lammps
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)

        fp_params = jdata['lammps_params']
        model_dir = fp_params['model_dir']
        model_dir = os.path.abspath(model_dir)
        model_name =fp_params['model_name']
        if not model_name :
            models = glob.glob(os.path.join(model_dir, '*pb'))
            model_name = [os.path.basename(ii) for ii in models]
        else:
            models = [os.path.join(model_dir,ii) for ii in model_name]
        forward_files = ['conf.lmp', 'lammps.in','strain.out']+model_name
        backward_files = ['log.lammps', 'model_devi.log']
        common_files=['lammps.in']+model_name

        if len(model_name)>1 and task_type == 'deepmd':
            backward_files = backward_files + ['model_devi.out']

    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    run_tasks = util.collect_task(all_task,task_type)
    machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)
    disp = make_dispatcher(machine)
    disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  common_files,
                  forward_files,
                  backward_files,
                  outlog='autotest.out',
                  errlog='autotest.err')

def cmpt_elastic(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    if task_type == "vasp":
        cmpt_02_elastic.cmpt_vasp(jdata, conf_dir)
    elif task_type in lammps_task_type:
        cmpt_02_elastic.cmpt_deepmd_lammps(jdata, conf_dir, task_type)
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
    elif task_type in lammps_task_type:
        gen_03_vacancy.make_lammps(jdata, conf_dir, task_type, supercell)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_vacancy(task_type,jdata,mdata):

    work_path=util.make_work_path(jdata,'03.vacancy',False,False,False)
    all_task = glob.glob(os.path.join(work_path,'struct-*'))

    #vasp
    if task_type == "vasp":
        mdata=decide_fp_machine(mdata)

        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR','OSZICAR']
        common_files=['INCAR','POTCAR']

    #lammps
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)

        fp_params = jdata['lammps_params']
        model_dir = fp_params['model_dir']
        model_dir = os.path.abspath(model_dir)
        model_name =fp_params['model_name']
        if not model_name :
            models = glob.glob(os.path.join(model_dir, '*pb'))
            model_name = [os.path.basename(ii) for ii in models]
        else:
            models = [os.path.join(model_dir,ii) for ii in model_name]
        common_files = model_name
        forward_files = ['conf.lmp', 'lammps.in']+model_name
        backward_files = ['log.lammps','model_devi.log']
        common_files=['lammps.in']+model_name

        if len(model_name)>1 and task_type == 'deepmd':
            backward_files = backward_files + ['model_devi.out']

    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    run_tasks = util.collect_task(all_task,task_type)
    machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)
    disp = make_dispatcher(machine)
    disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  common_files,
                  forward_files,
                  backward_files,
                  outlog='autotest.out',
                  errlog='autotest.err')

def cmpt_vacancy(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    supercell=jdata['supercell']
    #vasp
    if task_type == "vasp":
        cmpt_03_vacancy.cmpt_vasp(jdata, conf_dir, supercell)
    #lammps
    elif task_type in lammps_task_type:
        cmpt_03_vacancy.cmpt_deepmd_lammps(jdata, conf_dir, supercell, task_type)
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
    #lammps
    elif task_type in lammps_task_type:
        if not reprod_opt:
            gen_04_interstitial.make_lammps(jdata, conf_dir, supercell, insert_ele, task_type)
        else :
            gen_04_interstitial.make_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_interstitial(task_type,jdata,mdata):

    reprod_opt=jdata['reprod-opt']
    work_path=util.make_work_path(jdata,'04.interstitial',reprod_opt,False,False)
    all_task = glob.glob(os.path.join(work_path,'struct-*'))

    #vasp
    if task_type == "vasp":
        mdata=decide_fp_machine(mdata)

        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR','XDATCAR','OSZICAR']
        common_files=['INCAR']

    #lammps
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)

        if reprod_opt:
            all_frame=[]
            for ii in all_task:
                all_frame+=(glob.glob(os.path.join(ii,'frame.*')))
            work_path = all_task
            all_task = all_frame

        run_tasks_ = []
        for ii in all_task:
            fres = os.path.join(ii, 'log.lammps')
            if os.path.isfile(fres) :
                if not lammps.check_finished(fres):
                    run_tasks_.append(ii)
            else :
                run_tasks_.append(ii)

        fp_params = jdata['lammps_params']
        model_dir = fp_params['model_dir']
        model_dir = os.path.abspath(model_dir)
        model_name =fp_params['model_name']
        if not model_name :
            models = glob.glob(os.path.join(model_dir, '*pb'))
            model_name = [os.path.basename(ii) for ii in models]
        else:
            models = [os.path.join(model_dir,ii) for ii in model_name]
        forward_files = ['conf.lmp', 'lammps.in']+model_name
        backward_files = ['log.lammps', 'model_devi.log']
        common_files=['lammps.in']+model_name

        if len(model_name)>1 and task_type == 'deepmd':
            backward_files = backward_files + ['model_devi.out']

    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)
    disp = make_dispatcher(machine)
    if reprod_opt:
        for ii in work_path:
            run_tasks=[]
            for jj in run_tasks_:
                if ii in jj:
                    run_tasks.append(os.path.basename(jj))

            disp.run_jobs(resources,
                          command,
                          ii,
                          run_tasks,
                          group_size,
                          common_files,
                          forward_files,
                          backward_files,
                          outlog='autotest.out',
                          errlog='autotest.err')
    else:
        run_tasks = util.collect_task(all_task,task_type)
        disp.run_jobs(resources,
                      command,
                      work_path,
                      run_tasks,
                      group_size,
                      common_files,
                      forward_files,
                      backward_files,
                      outlog='autotest.log',
                      errlog='autotest.log')

def cmpt_interstitial(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    supercell=jdata['supercell']
    insert_ele=jdata['insert_ele']
    reprod_opt=jdata['reprod-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        cmpt_04_interstitial.cmpt_vasp(jdata, conf_dir, supercell, insert_ele)
    #lammps
    elif task_type in lammps_task_type:
        if not reprod_opt:
            cmpt_04_interstitial.cmpt_deepmd_lammps(jdata, conf_dir, supercell, insert_ele, task_type)
        else :
            task_name=task_type+'-reprod'
            cmpt_04_interstitial.cmpt_deepmd_reprod_traj(jdata, conf_dir, supercell, insert_ele, task_name)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def gen_surf(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    max_miller=jdata['max_miller']
    relax_box=jdata['relax_box']
    static=jdata['static-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_05_surf.make_vasp(jdata, conf_dir, max_miller, static = static, relax_box = relax_box)
    #lammps
    elif task_type in lammps_task_type :
        gen_05_surf.make_lammps(jdata, conf_dir, max_miller, static = static, relax_box = relax_box, task_type = task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_surf(task_type,jdata,mdata):
    static=jdata['static-opt']
    work_path=util.make_work_path(jdata,'05.surf',False,static,False)

    all_task = glob.glob(os.path.join(work_path,'struct-*'))

    #vasp
    if task_type == "vasp":
        mdata=decide_fp_machine(mdata)

        forward_files = ['INCAR', 'POSCAR','POTCAR']
        backward_files = ['OUTCAR','OSZICAR']
        common_files=['INCAR','POTCAR']

    #lammps
    elif task_type in lammps_task_type:
        mdata = decide_model_devi_machine(mdata)

        fp_params = jdata['lammps_params']
        model_dir = fp_params['model_dir']
        model_dir = os.path.abspath(model_dir)
        model_name =fp_params['model_name']
        if not model_name :
            models = glob.glob(os.path.join(model_dir, '*pb'))
            model_name = [os.path.basename(ii) for ii in models]
        else:
            models = [os.path.join(model_dir,ii) for ii in model_name]
        forward_files = ['conf.lmp', 'lammps.in']+model_name
        backward_files = ['log.lammps','model_devi.log']
        common_files=['lammps.in']+model_name

        if len(model_name)>1 and task_type == 'deepmd':
            backward_files = backward_files + ['model_devi.out']

    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

    run_tasks = util.collect_task(all_task,task_type)
    machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)
    disp = make_dispatcher(machine)
    disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  common_files,
                  forward_files,
                  backward_files,
                  outlog='autotest.out',
                  errlog='autotest.err')

def cmpt_surf(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    static_opt=jdata['static-opt']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        cmpt_05_surf.cmpt_vasp(jdata, conf_dir, static = static_opt)
    #lammps
    elif task_type in lammps_task_type :
        if static_opt:
            task_name =task_type+'-static'
        else:
            task_name =task_type
        cmpt_05_surf.cmpt_deepmd_lammps(jdata, conf_dir, task_name, static = static_opt)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def gen_phonon(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        gen_06_phonon.make_vasp(jdata, conf_dir)
    #lammps
    elif task_type in lammps_task_type:
        gen_06_phonon.make_lammps(jdata, conf_dir,  task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_phonon(task_type,jdata,mdata):
    user= ('user_incar' in jdata.keys())
    work_path=util.make_work_path(jdata,'06.phonon',False,False,user)

    all_task = glob.glob(os.path.join(work_path,'.'))

    #vasp
    if task_type == "vasp":
        mdata=decide_fp_machine(mdata)
        machine,machine_type,ssh_sess,resources,command,group_size=util.get_machine_info(mdata,task_type)

        run_tasks = util.collect_task(all_task,task_type)
        forward_files = ['INCAR', 'POTCAR','KPOINTS']
        backward_files = ['OUTCAR','OSZICAR','vasprun.xml']
        common_files=['POSCAR']

        disp = make_dispatcher(machine)
        disp.run_jobs(resources,
                  command,
                  work_path,
                  run_tasks,
                  group_size,
                  common_files,
                  forward_files,
                  backward_files,
                  outlog='autotest.out',
                  errlog='autotest.err')
    #lammps
    elif task_type in lammps_task_type:
        None
    else:
        raise RuntimeError ("unknow task %s, something wrong" % task_type)

def cmpt_phonon(task_type,jdata,mdata):
    conf_dir=jdata['conf_dir']
    cwd=os.getcwd()
    #vasp
    if task_type == "vasp":
        cmpt_06_phonon.cmpt_vasp(jdata, conf_dir)
    #lammps
    elif task_type in lammps_task_type :
        cmpt_06_phonon.cmpt_lammps(jdata,conf_dir, task_type)
    else :
        raise RuntimeError("unknow task ", task_type)
    os.chdir(cwd)

def run_task (json_file, machine_file) :
    with open (json_file, 'r') as fp :
        jdata = json.load (fp)
    with open (machine_file, 'r') as fp:
        mdata = json.load (fp)

    record = "record.auto_test"

    confs = jdata['conf_dir']
    ele_list=[key for key in jdata['potcar_map'].keys()]
    key_id = jdata['key_id']

    ii = jdata['task_type']
    jj=jdata['task']
    task_list=['equi','eos','elastic','vacancy','interstitial','surf','phonon','all']
    task_type_list=['vasp']+lammps_task_type
    if jj not in task_list :
        raise RuntimeError ("unknow task %s, something wrong" % jj)
    if ii not in task_type_list :
        raise RuntimeError ("unknow task type %s, something wrong" % ii)

    #gen_configuration
    if 'confs' in confs and (not os.path.exists(confs+'/POSCAR')) :
        print('generate %s' % (ele_list))
        if len(ele_list) == 1 :
                gen_confs.gen_element(ele_list[0],key_id)
        else :
                gen_confs.gen_alloy(ele_list,key_id)
    #default task
    log_iter ("gen_equi", ii, "equi")
    gen_equi (ii, jdata, mdata)
    log_iter ("run_equi", ii, "equi")
    run_equi  (ii, jdata, mdata)
    log_iter ("cmpt_equi", ii,"equi")
    cmpt_equi (ii, jdata, mdata)
    if  jj == "eos" or jj=="all":
        log_iter ("gen_eos", ii, "eos")
        gen_eos (ii, jdata, mdata)
        log_iter ("run_eos", ii, "eos")
        run_eos  (ii, jdata, mdata)
        log_iter ("cmpt_eos", ii, "eos")
        cmpt_eos (ii, jdata, mdata)
    if jj=="elastic" or jj=="all":
        log_iter ("gen_elastic", ii, "elastic")
        gen_elastic (ii, jdata, mdata)
        log_iter ("run_elastic", ii, "elastic")
        run_elastic  (ii, jdata, mdata)
        log_iter ("cmpt_elastic", ii, "elastic")
        cmpt_elastic (ii, jdata, mdata)
    if jj=="vacancy" or jj=="all":
        log_iter ("gen_vacancy", ii, "vacancy")
        gen_vacancy (ii, jdata, mdata)
        log_iter ("run_vacancy", ii, "vacancy")
        run_vacancy  (ii, jdata, mdata)
        log_iter ("cmpt_vacancy", ii, "vacancy")
        cmpt_vacancy (ii, jdata, mdata)
    if jj=="interstitial" or jj=="all":
        log_iter ("gen_interstitial", ii, "interstitial")
        gen_interstitial (ii, jdata, mdata)
        log_iter ("run_interstitial", ii, "interstitial")
        run_interstitial  (ii, jdata, mdata)
        log_iter ("cmpt_interstitial", ii, "interstitial")
        cmpt_interstitial (ii, jdata, mdata)
    if jj=="surf" or jj=="all":
        log_iter ("gen_surf", ii, "surf")
        gen_surf (ii, jdata, mdata)
        log_iter ("run_surf", ii, "surf")
        run_surf  (ii, jdata, mdata)
        log_iter ("cmpt_surf", ii, "surf")
        cmpt_surf (ii, jdata, mdata)
    '''
    if jj=="phonon":
        log_iter ("gen_phonon", ii, "phonon")
        gen_phonon (ii, jdata, mdata)
        log_iter ("run_phonon", ii, "phonon")
        run_phonon  (ii, jdata, mdata)
        log_iter ("cmpt_phonon", ii, "phonon")
        cmpt_phonon (ii, jdata, mdata)
    '''
    record_iter (record, confs, ii, jj)

def gen_test(args):
    logging.info ("start auto-testing")
    run_task (args.PARAM, args.MACHINE)
    logging.info ("finished!")


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

#/usr/bin/env python
# coding: utf-8
# Copyright (c) The Dpmodeling Team.

import os
import time
import json
from uuid import uuid4
from threading import Thread
from glob import glob
from dpgen import dlog
from dpgen import SHORT_CMD
from dpgen.database.entry import Entry
from dpgen.database.vasp import VaspInput
from dpdata import System,LabeledSystem
from monty.serialization import loadfn,dumpfn
import numpy as np
import traceback

OUTPUT=SHORT_CMD+'_db.json'
SUPPORTED_CACULATOR=['vasp','pwscf','gaussian']
ITERS_PAT="iter.*/02.fp/task*"
INIT_PAT="init/*/02.md/sys-*/scale-*/*"

def db_run(args):
    dlog.info ("collecting data")
    #print(args.ID_PREFIX)
    _main(args.PARAM)
    dlog.info ("finished")

def _main(param):
    with open(param, "r") as fp:
      jdata = json.load(fp)
    calculator = jdata["calculator"]
    path = jdata["path"]
    calulator = jdata["calculator"]
    output = jdata["output"]
    config_info_dict = jdata["config_info_dict"]
    id_prefix = jdata["id_prefix"]
    skip_init = False
    if "skip_init" in jdata:
      skip_init = jdata["skip_init"]
    ## The mapping from sys_info to sys_configs
    assert calculator.lower() in SUPPORTED_CACULATOR
    dlog.info('data collection from: %s'%path)
    if calculator == "vasp":
        parsing_vasp(path,config_info_dict,skip_init, output,id_prefix)
    elif calculator == 'gaussian': 
        parsing_gaussian(path,output)
    else:
        parsing_pwscf(path,output)

def parsing_vasp(path,config_info_dict, skip_init, output=OUTPUT,id_prefix=None):
    
    fp_iters=os.path.join(path,ITERS_PAT) 
    dlog.debug(fp_iters)
    f_fp_iters=glob(fp_iters)
    dlog.info("len iterations data: %s"%len(f_fp_iters))
    fp_init=os.path.join(path,INIT_PAT)
    dlog.debug(fp_init)
    f_fp_init=glob(fp_init)
    if skip_init:
      entries = _parsing_vasp(f_fp_iters,config_info_dict, id_prefix)
      dlog.info("len collected data: %s"%len(entries))
    else:
      dlog.info("len initialization data: %s"%len(f_fp_init))
      entries=_parsing_vasp(f_fp_init,config_info_dict, id_prefix,iters=False)
      entries.extend(_parsing_vasp(f_fp_iters,config_info_dict, id_prefix))
      dlog.info("len collected data: %s"%len(entries))
    #print(output)
    #print(entries)
    dumpfn(entries,output,indent=4) 

def _parsing_vasp(paths,config_info_dict, id_prefix,iters=True):
    entries=[]
    icount=0
    if iters:
      iter_record = []
      iter_record_new = []
      try:
        with open ("record.database", "r") as f_record:
          iter_record = [i.split()[0] for i in f_record.readlines()]
        iter_record.sort()
        dlog.info("iter_record")
        dlog.info(iter_record)
      except Exception:
        pass
    for path in paths:
      try:
        f_outcar = os.path.join(path,'OUTCAR')
        f_job = os.path.join(path,'job.json')
        tmp_iter = path.split('/')[-3]
        if (tmp_iter in iter_record) and (tmp_iter != iter_record[-1]):
          continue
        if tmp_iter not in iter_record_new:
          iter_record_new.append(tmp_iter)
        vi = VaspInput.from_directory(path)
        if os.path.isfile(f_job):
          attrib=loadfn(f_job)
        else:
          attrib={}

        if iters and attrib:
          # generator/Cu/iter.000031/02.fp/task.007.000000
          tmp_=path.split('/')[-1]
          #config_info=tmp_.split('.')[1]
          task_info=tmp_.split('.')[-1]
          tmp_iter = path.split('/')[-3]
          iter_info = tmp_iter.split('.')[-1]
          sys_info = path.split('/')[-4]
          config_info_int = int(tmp_.split('.')[1])
          for (key, value) in config_info_dict.items():
            if config_info_int in value:
              config_info = key
          attrib['config_info']=config_info
          attrib['task_info']=task_info
          attrib['iter_info']=iter_info
          attrib['sys_info']=sys_info
          with open(f_outcar , "r") as fin_outcar:
            infile_outcar = fin_outcar.readlines()
          for line in infile_outcar:
            if "running on" in line:
              attrib["core"] = int(line.split()[2])
            if "Elapse" in line:
              attrib["wall_time"] = float(line.split()[-1])
            if "executed on" in line:
              attrib["date"] = line.split()[-2]
              attrib["clocktime"] = line.split()[-1]
          dlog.info("Attrib")
          dlog.info(attrib)
        comp=vi['POSCAR'].structure.composition
        ls = LabeledSystem(f_outcar)
        lss=ls.to_list()
        for ls in lss:
          if id_prefix:
            eid=id_prefix+"_"+str(icount)
          else:
            eid = str(uuid4())
        entry=Entry(comp,'vasp',vi.as_dict(),ls.as_dict(),attribute=attrib,entry_id=eid)
        entries.append(entry)
        icount+=1
      except Exception:
        #dlog.info(str(Exception))
        dlog.info("failed for %s"%(path))
        #pass
    if iters:
      iter_record.sort()
      iter_record_new.sort()
      with open("record.database" , "w") as fw:
        for line in iter_record:
          fw.write(line + "\n")
        for line in iter_record_new:
          fw.write(line + "\n")  
    return entries

def parsing_pwscf(path,output=OUTPUT):
    pass

def parsing_gaussian(path,output=OUTPUT):
    pass


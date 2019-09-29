#/usr/bin/env python
# coding: utf-8
# Copyright (c) The Dpmodeling Team.

import os
import time
from uuid import uuid4
from threading import Thread
from glob import glob
from dpgen import dlog
from dpgen import SHORT_CMD
from dpgen.database.entry import Entry
from dpgen.database.vasp import VaspInput
from dpdata import System,LabeledSystem
from monty.serialization import loadfn,dumpfn

OUTPUT=SHORT_CMD+'_db.json'
SUPPORTED_CACULATOR=['vasp','pwscf','siesta','gaussian']
ITERS_PAT="iter.*/02.fp/task*"
INIT_PAT="init/*/02.md/sys-*/scale-*/*"

def db_run(args):
    dlog.info ("collecting data")
    print(args.ID_PREFIX)
    _main(args.PATH, args.CALCULATOR, args.OUTPUT,args.ID_PREFIX)
    dlog.info ("finished")

def _main(path,calculator,output,id_prefix):
    assert calculator.lower() in SUPPORTED_CACULATOR
    dlog.info('data collection from: %s'%path)
    if calculator == "vasp":
        parsing_vasp(path,output,id_prefix)
    elif calculator == 'gaussian': 
        parsing_gaussian(path,output)
    elif calculator == "siesta":
        parsing_siesta(path, output)
    else:
        parsing_pwscf(path,output)

def parsing_vasp(path,output=OUTPUT,id_prefix=None):
    
    fp_iters=os.path.join(path,ITERS_PAT) 
    dlog.debug(fp_iters)
    f_fp_iters=glob(fp_iters)
    dlog.info("len iterations data: %s"%len(f_fp_iters))
    fp_init=os.path.join(path,INIT_PAT)
    dlog.debug(fp_init)
    f_fp_init=glob(fp_init)
    dlog.info("len initialization data: %s"%len(f_fp_init))
    entries=_parsing_vasp(f_fp_init,id_prefix,iters=False)
    entries.extend(_parsing_vasp(f_fp_iters,id_prefix))
    dlog.info("len collected data: %s"%len(entries))
   
    dumpfn(entries,output,indent=4) 

def _parsing_vasp(paths,id_prefix,iters=True):
    entries=[]
    icount=0
    for path in paths:
        f_outcar = os.path.join(path,'OUTCAR')
        f_job = os.path.join(path,'job.json')
        
        try: 
           vi = VaspInput.from_directory(path)
           if os.path.isfile(f_job):
              attrib=loadfn(f_job)
           else:
              attrib={}

           if iters and attrib:
              tmp_=path.split('/')[-1]
              iter_info=tmp_.split('.')[1]
              task_info=tmp_.split('.')[-1]
              attrib['iter_info']=iter_info
              attrib['task_info']=task_info
           else:
              pass
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
        except:
           dlog.info("failed here : %s"%path)
    return entries

def parsing_pwscf(path,output=OUTPUT):
    pass
def parsing_siesta(path,output=OUTPUT):
    pass
def parsing_gaussian(path,output=OUTPUT):
    pass


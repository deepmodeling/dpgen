#!/usr/bin/env python3

import logging

from monty.serialization import loadfn

from dpgen import dlog
from dpgen.auto_test.common_equi import make_equi, run_equi, post_equi
from dpgen.auto_test.common_prop import make_property, run_property, post_property


#lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']

def run_task(step, param_file, machine_file=None):
    jdata=loadfn(param_file)
    confs = jdata['structures']
    inter_parameter = jdata['interaction']

    if machine_file:
       mdata=loadfn(machine_file)

    if step == 'make' and 'relaxation' in jdata:
        relax_param = jdata['relaxation']
        make_equi(confs, inter_parameter, relax_param)

    elif step == 'make' and 'properties' in jdata:
        property_list = jdata['properties']
        make_property(confs, inter_parameter, property_list)

    elif step == 'run' and 'relaxation' in jdata:
        if machine_file is None:
           print('Please supply the machine.json, exit now!')  
           return
        run_equi(confs, inter_parameter, mdata)

    elif step == 'run' and 'properties' in jdata:
        if machine_file is None:
           print('Please supply the machine.json, exit now!')  
           return
        property_list = jdata['properties']
        run_property(confs, inter_parameter, property_list, mdata)

    elif step == 'post' and 'relaxation' in jdata:
        post_equi(confs, inter_parameter)

    elif step == 'post' and 'properties' in jdata:
        property_list = jdata['properties']
        post_property(confs,inter_parameter, property_list)

    else:
        raise RuntimeError('unknown tasks')

def gen_test(args):
    dlog.info("start auto-testing")
    if args.debug:
        dlog.setLevel(logging.DEBUG)
    run_task(args.TASK, args.PARAM, args.MACHINE)
    dlog.info("finished!")




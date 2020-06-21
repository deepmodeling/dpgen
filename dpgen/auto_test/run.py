#!/usr/bin/env python3

import argparse, json
import logging
from dpgen.auto_test.common_prop import make_property
from dpgen.auto_test.common_equi import make_equi,run_equi,post_equi

lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']


def run_task(step, json_file, machine_file=''):
    with open(json_file, 'r') as fp:
        jdata = json.load(fp)

    confs = jdata['structures']
    inter_parameter = jdata['interaction']

    if step == 'make' and 'relaxation' in jdata:
        relax_param = jdata['relaxation']
        make_equi(confs, inter_parameter, relax_param)

    elif step == 'make' and 'properties' in jdata:
        property_list = jdata['properties']
        make_property(confs, inter_parameter, property_list)

    elif step == 'run' and 'relaxation' in jdata:
        with open(machine_file, 'r') as fp:
            mdata = json.load(fp)
        run_equi(confs, inter_parameter, mdata)

    elif step == 'run' and 'properties' in jdata:
        with open(machine_file, 'r') as fp:
            mdata = json.load(fp)
        property_list = jdata['properties']
        run_property(confs, inter_parameter, property_list, mdata)

    elif step == 'post' and 'relaxation' in jdata:
        post_equi(confs, inter_parameter)

    elif step == 'post' and 'properties' in jdata:
        property_list = jdata['properties']
        post_property(confs, property_list)

    else:
        raise RuntimeError('unknown tasks')


def gen_test(args):
    logging.info("start auto-testing")
    run_task(args.TASK, args.PARAM, args.MACHINE)
    logging.info("finished!")


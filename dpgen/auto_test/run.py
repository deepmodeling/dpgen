#!/usr/bin/env python3

import argparse, json
import logging
from dpgen.auto_test import common

lammps_task_type = ['deepmd', 'meam', 'eam_fs', 'eam_alloy']


def run_task(step, json_file, machine_file=''):
    with open(json_file, 'r') as fp:
        jdata = json.load(fp)

    confs = jdata['structures']
    inter_parameter = jdata['interaction']

    if step == 'make' and 'relaxation' in jdata:
        relax_param = jdata['relaxation']
        common.make_equi(confs, inter_parameter, relax_param)

    elif step == 'make' and 'properties' in jdata:
        property_list = jdata['properties']
        common.make_property(confs, inter_parameter, property_list)

    elif step == 'run' and 'relaxation' in jdata:
        with open(machine_file, 'r') as fp:
            mdata = json.load(fp)
        common.run_equi(confs, inter_parameter, mdata)

    elif step == 'run' and 'properties' in jdata:
        with open(machine_file, 'r') as fp:
            mdata = json.load(fp)
        property_list = jdata['properties']
        common.run_property(confs, inter_parameter, property_list, mdata)

    elif step == 'post' and 'relaxation' in jdata:
        common.post_equi(confs, inter_parameter)

    elif step == 'post' and 'properties' in jdata:
        property_list = jdata['properties']
        common.post_property(confs, property_list)

    else:
        raise RuntimeError('unknown tasks')


def gen_test(args):
    logging.info("start auto-testing")
    run_task(args.STEP, args.PARAM, args.MACHINE)
    logging.info("finished!")


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("STEP", type=str,
                        help="Make or run or post")
    parser.add_argument("PARAM", type=str,
                        help="The parameters of the generator")
    parser.add_argument("MACHINE", type=str,
                        help="The settings of the machine running the generator",
                        default=None)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
    # logging.basicConfig (filename="compute_string.log", filemode="a", level=logging.INFO, format='%(asctime)s %(message)s')
    logging.getLogger("paramiko").setLevel(logging.WARNING)

    logging.info("start auto-testing")
    run_task(args.STEP, args.PARAM, args.MACHINE)
    logging.info("finished!")


if __name__ == '__main__':
    _main()

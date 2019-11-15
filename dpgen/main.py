#!/usr/bin/env python
# coding: utf-8
# Copyright (c) DeepGenerator Development Team.


import argparse
import sys
import itertools
from dpgen.generator.run import gen_run
from dpgen.data.gen import gen_init_bulk
from dpgen.data.surf import gen_init_surf
from dpgen.data.reaction import gen_init_reaction
from dpgen.auto_test.run import gen_test
from dpgen.database.run import db_run
from dpgen.tools.run_report import run_report
from dpgen import info, __version__, __date__



"""
A master convenience script with many tools for driving dpgen.
"""

__author__ = "Han Wang"
__copyright__ = "Copyright 2019, The DP-GEN Project"
__maintainer__ = "Haidi Wang"
__email__ = ""


def main():
    info()
    print("Description\n------------")
    parser = argparse.ArgumentParser(description="""
    dpgen is a convenient script that uses DeepGenerator to prepare initial
    data, drive DeepMDkit and analyze results. This script works based on
    several sub-commands with their own options. To see the options for the
    sub-commands, type "dpgen sub-command -h".""")

    subparsers = parser.add_subparsers()

    # init surf model
    parser_init_surf = subparsers.add_parser(
        "init_surf", help="Generating initial data for surface systems.")
    parser_init_surf.add_argument('PARAM', type=str, 
                             help="parameter file, json/yaml format")
    parser_init_surf.add_argument('MACHINE', type=str,default=None,nargs="?",
                        help="machine file, json/yaml format")
    parser_init_surf.set_defaults(func=gen_init_surf)
    
    # init bulk model
    parser_init_bulk = subparsers.add_parser(
        "init_bulk", help="Generating initial data for bulk systems.")
    parser_init_bulk.add_argument('PARAM', type=str, 
                             help="parameter file, json/yaml format")
    parser_init_bulk.add_argument('MACHINE', type=str,default=None,nargs="?",
                        help="machine file, json/yaml format")
    parser_init_bulk.set_defaults(func=gen_init_bulk)
    # parser_init.add_argument("-p",'--parameter', type=str, dest='param',
    #                     help="parameter file, json/yaml format")
    # parser_init.add_argument("-s","--stage", type=int, dest='stage',
    #                     help="the stage of init, can be 1, 2, 3 or 4. "
    #                     "1: Setup vasp jobs for relaxation. "
    #                     "2: Collect vasp relaxed confs (if relax is not skiped). Perturb system. "
    #                     "3: Setup vasp jobs for MD of perturbed system. "
    #                     "4: Collect vasp md confs, make deepmd data. ")
    # parser_init.add_argument("directories", metavar="dir", default=".",
    #                             type=str, nargs="*",
    #                             help="directory to process (default to .)")
    # parser_init.set_defaults(func=gen_data)

    parser_init_reaction = subparsers.add_parser(
        "init_reaction", help="Generating initial data for reactive systems.")
    parser_init_reaction.add_argument('PARAM', type=str, 
                             help="parameter file, json/yaml format")
    parser_init_reaction.add_argument('MACHINE', type=str,default=None,nargs="?",
                        help="machine file, json/yaml format")
    parser_init_reaction.set_defaults(func=gen_init_reaction)

    # run 
    parser_run = subparsers.add_parser(
        "run",
        help="Main process of Deep Potential Generator.")
    parser_run.add_argument('PARAM', type=str,
                        help="parameter file, json/yaml format")
    parser_run.add_argument('MACHINE', type=str,
                        help="machine file, json/yaml format")
    parser_run.add_argument('-d','--debug', action='store_true',
                        help="log debug info")
    parser_run.set_defaults(func=gen_run)

    # run/report
    parser_rr = subparsers.add_parser(
        "run/report",
        help="Report the systems and the thermodynamic conditions of the labeled frames.")
    parser_rr.add_argument("JOB_DIR", type=str, 
                           help="the directory of the DP-GEN job,")
    parser_rr.add_argument('-s',"--stat-sys", action = 'store_true',
                           help="count the labeled frames for each system")
    parser_rr.add_argument('-i', "--stat-iter", action= 'store_true',
                            help="print the iteration candidate,failed,accurate count and fp calculation,success and fail count")
    parser_rr.add_argument('-t', "--stat-time", action= 'store_true',
                            help="print the iteration time, warning!! assume model_devi parallel cores == 1")
    parser_rr.add_argument('-p',"--param", type=str, default = 'param.json',
                           help="the json file provides DP-GEN paramters, should be located in JOB_DIR")
    parser_rr.add_argument('-v',"--verbose", action = 'store_true',
                           help="being loud")
    parser_rr.set_defaults(func=run_report)    

    # test 
    parser_test = subparsers.add_parser("test", help="Auto-test for Deep Potential.")
    parser_test.add_argument('PARAM', type=str,
                        help="parameter file, json/yaml format")
    parser_test.add_argument('MACHINE', type=str,
                        help="machine file, json/yaml format")
    parser_test.set_defaults(func=gen_test)    

    # db 
    parser_db = subparsers.add_parser(
        "db",
        help="Collecting data from DP-GEN.")

    parser_db.add_argument('PARAM', type=str,
                        help="parameter file, json format")

    parser_db.set_defaults(func=db_run)


    try:
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError:
        # argcomplete not present.
        pass

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()

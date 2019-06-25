#!/usr/bin/env python
# coding: utf-8
# Copyright (c) DeepGenerator Development Team.


import argparse
import sys
import itertools
from dpgen.generator.run import gen_run
from dpgen.data.gen import gen_init
from dpgen import info


"""
A master convenience script with many tools for driving dpgen.
"""

__author__ = "Han Wang"
__copyright__ = "Copyright 2019, The DP-GEN Project"
__version__ = "0.1.0"
__maintainer__ = "Haidi Wang"
__email__ = ""
__date__ = "2019.06.26"


def main():
    info()
    parser = argparse.ArgumentParser(description="""
    dpgen is a convenient script that uses DeepGenerator to prepare initial
    data, drive DeepMDkit and analyze results. This script works based on
    several sub-commands with their own options. To see the options for the
    sub-commands, type "dpgen sub-command -h".""",
                                     epilog="""
    Author: DeepGenTeam
    Version: {}
    Last updated: {}""".format(__version__, __date__))

    subparsers = parser.add_subparsers()
    
    # init model
    parser_init = subparsers.add_parser(
        "init", help="dpgen initial data preparation tools.")
    parser_init.add_argument('PARAM', type=str, 
                             help="parameter file, json format")
    parser_init.add_argument('STAGE', type=int,
                             help="the stage of init, can be 1, 2, 3 or 4. "
                             "1: Setup vasp jobs for relaxation. "
                             "2: Collect vasp relaxed confs (if relax is not skiped). Perturb system. "
                             "3: Setup vasp jobs for MD of perturbed system. "
                             "4: Collect vasp md confs, make deepmd data. "
    )
    parser_init.set_defaults(func=gen_init)
    # parser_init.add_argument("-p",'--parameter', type=str, dest='param',
    #                     help="parameter file, json format")
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

    # run model
    parser_run = subparsers.add_parser(
        "run",
        help="Runing DeepMD with generator model.")
    parser_run.add_argument('PARAM', type=str,
                        help="parameter file, json format")
    parser_run.add_argument('MACHINE', type=str,
                        help="machine file, json format")
    parser_run.set_defaults(func=gen_run)

    # # test model
    # parser_test = subparsers.add_parser("test", help="auto test for deep potential.")
    # parser_test.add_argument("-p",'--parameter', type=str, dest='param',
    #                     help="parameter file, json format")
    # parser_test.add_argument("-m",'--machine', type=str, dest='machine',
    #                     help="machine file, json format")
    # parser_test.set_defaults(func=gen_test)

    # # convert  model
    # parser_structure = subparsers.add_parser(
    #     "struct",
    #     help="structure conversion and analysis tools.")

    # parser_structure.add_argument(
    #     "-f", "--filenames", dest="filenames",
    #     metavar="filename", nargs="+",
    #     help="List of structure files.")

    # groups = parser_structure.add_mutually_exclusive_group(required=True)
    # groups.add_argument("-c", "--convert", dest="convert", action="store_true",
    #                     help="Convert from structure file 1 to structure "
    #                          "file 2. Format determined from filename. "
    #                          "Supported formats include POSCAR/CONTCAR, "
    #                          "CIF, lmp.")
    # groups.add_argument("-s", "--symmetry", dest="symmetry",
    #                     metavar="tolerance", type=float,
    #                    help="Determine the spacegroup using the "
    #                          "specified tolerance. 0.1 is usually a good "
    #                          "value for DFT calculations.")
    # parser_structure.set_defaults(func=gen_struct)


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

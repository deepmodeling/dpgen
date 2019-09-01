#!/usr/bin/env python
# coding: utf-8
# Copyright (c) DeepGenerator Development Team.


import argparse
import sys
import itertools
from dpgen.generator.run import gen_run
from dpgen.data.gen import gen_init_bulk
from dpgen.data.surf import gen_init_surf
from dpgen.auto_test.run import gen_test
from dpgen.database.run import db_run
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

    # init surf model
    parser_init_surf = subparsers.add_parser(
        "init_surf", help="Generating initial data for surface systems.")
    parser_init_surf.add_argument('PARAM', type=str, 
                             help="parameter file, json/yaml format")
    parser_init_surf.add_argument('MACHINE', type=str,
                        help="machine file, json/yaml format")
    parser_init_surf.set_defaults(func=gen_init_surf)
    
    # init bulk model
    parser_init_bulk = subparsers.add_parser(
        "init_bulk", help="Generating initial data for bulk systems.")
    parser_init_bulk.add_argument('PARAM', type=str, 
                             help="parameter file, json/yaml format")
    parser_init_bulk.add_argument('MACHINE', type=str,
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

    # run 
    parser_run = subparsers.add_parser(
        "run",
        help="Main process of Deep Generator.")
    parser_run.add_argument('PARAM', type=str,
                        help="parameter file, json/yaml format")
    parser_run.add_argument('MACHINE', type=str,
                        help="machine file, json/yaml format")
    parser_run.set_defaults(func=gen_run)

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
    parser_db.add_argument('PATH', type=str,
                        help="root path for dpgen modeling")
    parser_db.add_argument('CALCULATOR', type=str,
                        help="calculator used for labeling: vasp/pwscf/gaussian")
    parser_db.add_argument('OUTPUT', type=str,
                        help="output filename : file.json/file.yaml")
    parser_db.add_argument("ID_PREFIX", type=str, default=None,
                                 nargs="?",
                                 help="prefix of an  entry id")

    parser_db.set_defaults(func=db_run)

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

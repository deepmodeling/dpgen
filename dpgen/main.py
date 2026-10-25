#!/usr/bin/env python
# Copyright (c) DeepGenerator Development Team.


import argparse
import importlib
import sys
from functools import partial

from dpgen import info
from dpgen.skills import configure_parser, show_skills, workflow_help

"""
A master convenience script with many tools for driving dpgen.
"""

__author__ = "Han Wang"
__copyright__ = "Copyright 2019, The DP-GEN Project"
__maintainer__ = "Haidi Wang"
__email__ = ""


def _run_command(module, name, args):
    """Load only the requested workflow when it is actually executed."""
    return getattr(importlib.import_module(module), name)(args)


def main_parser() -> argparse.ArgumentParser:
    """Returns parser for `dpgen` command.

    Returns
    -------
    argparse.ArgumentParser
        parser for `dpgen` command
    """
    parser = argparse.ArgumentParser(
        description="""
    dpgen is a convenient script that uses DeepGenerator to prepare initial
    data, drive DeepMDkit and analyze results. This script works based on
    several sub-commands with their own options. To see the options for the
    sub-commands, type "dpgen sub-command -h"."""
    )

    parser.epilog = "Discover bundled workflow guides: dpgen skills list."
    subparsers = parser.add_subparsers()
    configure_parser(
        subparsers.add_parser(
            "skills", help="Read bundled offline workflow instructions."
        )
    )

    # init surf model
    parser_init_surf = subparsers.add_parser(
        "init_surf", help="Generating initial data for surface systems."
    )
    parser_init_surf.add_argument(
        "PARAM", type=str, help="parameter file, json/yaml format"
    )
    parser_init_surf.add_argument(
        "MACHINE",
        type=str,
        default=None,
        nargs="?",
        help="machine file, json/yaml format",
    )
    parser_init_surf.set_defaults(
        func=partial(_run_command, "dpgen.data.surf", "gen_init_surf")
    )

    # init bulk model
    parser_init_bulk = subparsers.add_parser(
        "init_bulk", help="Generating initial data for bulk systems."
    )
    parser_init_bulk.add_argument(
        "PARAM", type=str, help="parameter file, json/yaml format"
    )
    parser_init_bulk.add_argument(
        "MACHINE",
        type=str,
        default=None,
        nargs="?",
        help="machine file, json/yaml format",
    )
    parser_init_bulk.set_defaults(
        func=partial(_run_command, "dpgen.data.gen", "gen_init_bulk")
    )

    parser_auto_gen_param = subparsers.add_parser(
        "auto_gen_param", help="auto gen param.json"
    )
    # parser_auto_gen_param.add_argument('meltpoint', type=float, help="melt point")
    parser_auto_gen_param.add_argument(
        "PARAM", type=str, help="parameter file, json/yaml format"
    )
    parser_auto_gen_param.set_defaults(
        func=partial(_run_command, "dpgen.tools.auto_gen_param", "auto_gen_param")
    )

    # parser_init_reaction
    parser_init_reaction = subparsers.add_parser(
        "init_reaction", help="Generating initial data for reactive systems."
    )
    parser_init_reaction.add_argument(
        "PARAM", type=str, help="parameter file, json/yaml format"
    )
    parser_init_reaction.add_argument(
        "MACHINE",
        type=str,
        default=None,
        nargs="?",
        help="machine file, json/yaml format",
    )
    parser_init_reaction.set_defaults(
        func=partial(_run_command, "dpgen.data.reaction", "gen_init_reaction")
    )

    # run
    parser_run = subparsers.add_parser(
        "run",
        help="Main process of Deep Potential Generator.",
        epilog=workflow_help("dpgen-run"),
    )
    parser_run.add_argument("PARAM", type=str, help="parameter file, json/yaml format")
    parser_run.add_argument("MACHINE", type=str, help="machine file, json/yaml format")
    parser_run.add_argument("-d", "--debug", action="store_true", help="log debug info")
    parser_run.set_defaults(
        func=partial(_run_command, "dpgen.generator.run", "gen_run")
    )

    # run/report
    parser_rr = subparsers.add_parser(
        "run/report",
        help="Report the systems and the thermodynamic conditions of the labeled frames.",
    )
    parser_rr.add_argument("JOB_DIR", type=str, help="the directory of the DP-GEN job,")
    parser_rr.add_argument(
        "-s",
        "--stat-sys",
        action="store_true",
        help="count the labeled frames for each system",
    )
    parser_rr.add_argument(
        "-i",
        "--stat-iter",
        action="store_true",
        help="print the iteration candidate,failed,accurate count and fp calculation,success and fail count",
    )
    parser_rr.add_argument(
        "-t",
        "--stat-time",
        action="store_true",
        help="print the iteration time, warning!! assume model_devi parallel cores == 1",
    )
    parser_rr.add_argument(
        "-p",
        "--param",
        type=str,
        default="param.json",
        help="the json file provides DP-GEN paramters, should be located in JOB_DIR",
    )
    parser_rr.add_argument("-v", "--verbose", action="store_true", help="being loud")
    parser_rr.set_defaults(
        func=partial(_run_command, "dpgen.tools.run_report", "run_report")
    )

    # collect
    parser_coll = subparsers.add_parser("collect", help="Collect data.")
    parser_coll.add_argument(
        "JOB_DIR", type=str, help="the directory of the DP-GEN job"
    )
    parser_coll.add_argument("OUTPUT", type=str, help="the output directory of data")
    parser_coll.add_argument(
        "-p",
        "--parameter",
        type=str,
        default="param.json",
        help="the json file provides DP-GEN paramters, should be located in JOB_DIR",
    )
    parser_coll.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print number of data in each system",
    )
    parser_coll.add_argument(
        "-m",
        "--merge",
        action="store_true",
        help="merge the systems with the same chemical formula",
    )
    parser_coll.add_argument(
        "-s", "--shuffle", action="store_true", help="shuffle the data systems"
    )
    parser_coll.set_defaults(
        func=partial(_run_command, "dpgen.collect.collect", "gen_collect")
    )

    # simplify
    parser_run = subparsers.add_parser(
        "simplify", help="Simplify data.", epilog=workflow_help("dpgen-simplify")
    )
    parser_run.add_argument("PARAM", type=str, help="parameter file, json/yaml format")
    parser_run.add_argument("MACHINE", type=str, help="machine file, json/yaml format")
    parser_run.add_argument("-d", "--debug", action="store_true", help="log debug info")
    parser_run.set_defaults(
        func=partial(_run_command, "dpgen.simplify.simplify", "gen_simplify")
    )

    # test
    parser_test = subparsers.add_parser(
        "autotest", help="Auto-test for Deep Potential."
    )
    parser_test.add_argument("TASK", type=str, help="task can be make, run or post")
    parser_test.add_argument("PARAM", type=str, help="parameter file, json/yaml format")
    parser_test.add_argument(
        "MACHINE",
        type=str,
        default=None,
        nargs="?",
        help="machine file, json/yaml format",
    )
    parser_test.add_argument(
        "-d", "--debug", action="store_true", help="log debug info"
    )
    parser_test.set_defaults(
        func=partial(_run_command, "dpgen.auto_test.run", "gen_test")
    )

    # db
    parser_db = subparsers.add_parser("db", help="Collecting data from DP-GEN.")

    parser_db.add_argument("PARAM", type=str, help="parameter file, json format")

    parser_db.set_defaults(func=partial(_run_command, "dpgen.database.run", "db_run"))

    # gui
    parser_gui = subparsers.add_parser(
        "gui",
        help="Serve DP-GUI.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser_gui.add_argument(
        "-p",
        "--port",
        type=int,
        default=6042,
        help="The port to serve DP-GUI on.",
    )
    parser_gui.add_argument(
        "--bind_all",
        action="store_true",
        help=(
            "Serve on all public interfaces. This will expose your DP-GUI instance "
            "to the network on both IPv4 and IPv6 (where available)."
        ),
    )
    parser_gui.set_defaults(func=partial(_run_command, "dpgen.gui", "start_dpgui"))
    return parser


def main():
    parser = main_parser()
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
    if args.func is not show_skills:
        info()
        print("Description\n------------")
    args.func(args)


if __name__ == "__main__":
    main()

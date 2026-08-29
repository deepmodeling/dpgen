#!/usr/bin/env python3

import argparse
import warnings

from dpgen.collect.collect import collect_data as collect_current_data


def collect_data(target_folder, param_file, output, verbose=True):
    """Delegate the legacy helper to the maintained collection implementation."""
    warnings.warn(
        "dpgen.tools.collect_data is deprecated; use `dpgen collect` or "
        "`dpgen.collect.collect.collect_data` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return collect_current_data(
        target_folder,
        param_file,
        output,
        verbose=verbose,
        shuffle=True,
        merge=False,
        include_init_data=False,
        iter_output_prefix="system.",
        discover_existing_iters=True,
    )


def _main():
    parser = argparse.ArgumentParser(
        description="Deprecated wrapper for `dpgen collect`"
    )
    parser.add_argument("JOB_DIR", type=str, help="the directory of the DP-GEN job")
    parser.add_argument("OUTPUT", type=str, help="the output directory of data")
    parser.add_argument(
        "-p",
        "--parameter",
        type=str,
        default="param.json",
        help="the json file provides DP-GEN paramters, should be located in JOB_DIR",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print number of data in each system",
    )
    args = parser.parse_args()

    collect_data(args.JOB_DIR, args.parameter, args.OUTPUT, args.verbose)


if __name__ == "__main__":
    _main()

# SPDX-License-Identifier: LGPL-3.0-or-later
"""DP-GUI entrypoint."""
import argparse


def start_dpgui(args: argparse.Namespace):
    """Host DP-GUI server.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments from argparse.

    Raises
    ------
    ModuleNotFoundError
        The dpgui package is not installed
    """
    try:
        from dpgui import (
            start_dpgui,
        )
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "To use DP-GUI, please install the dpgui package:\npip install dpgui"
        ) from e
    start_dpgui(port=args.port, bind_all=args.bind_all)

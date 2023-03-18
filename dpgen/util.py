#!/usr/bin/env python
# coding: utf-8
import json
import os
from pathlib import Path
from typing import List, Union

import dpdata
import h5py
from dargs import Argument

from dpgen import dlog

"""
some common utilities for generator, auto_test and data
"""

# constants define
MaxLength = 70


def sepline(ch="-", sp="-", screen=False):
    r"""
    seperate the output by '-'
    """
    if screen:
        print(ch.center(MaxLength, sp))
    else:
        dlog.info(ch.center(MaxLength, sp))


def box_center(ch="", fill=" ", sp="|"):
    r"""
    put the string at the center of |  |
    """
    strs = ch.center(Len, fill)
    dlog.info(sp + strs[1 : len(strs) - 1 :] + sp)


def expand_sys_str(root_dir: Union[str, Path]) -> List[str]:
    """Recursively iterate over directories taking those that contain `type.raw` file.

    If root_dir is a file but not a directory, it will be assumed as an HDF5 file.

    Parameters
    ----------
    root_dir : Union[str, Path]
        starting directory

    Returns
    -------
    List[str]
        list of string pointing to system directories
    """
    root_dir = Path(root_dir)
    if root_dir.is_dir():
        matches = [str(d) for d in root_dir.rglob("*") if (d / "type.raw").is_file()]
        if (root_dir / "type.raw").is_file():
            matches.append(str(root_dir))
    elif root_dir.is_file():
        # HDF5 file
        with h5py.File(root_dir, "r") as f:
            # list of keys in the h5 file
            f_keys = ["/"]
            f.visit(lambda x: f_keys.append("/" + x))
        matches = [
            "%s#%s" % (root_dir, d)
            for d in f_keys
            if str(Path(d) / "type.raw") in f_keys
        ]
    else:
        raise OSError(f"{root_dir} does not exist.")
    return matches


def normalize(arginfo: Argument, data: dict, strict_check: bool = True) -> dict:
    """Normalize and check input data.

    Parameters
    ----------
    arginfo : dargs.Argument
        argument information
    data : dict
        input data
    strict_check : bool, default=True
        strict check data or not

    Returns
    -------
    dict
        normalized data
    """
    data = arginfo.normalize_value(data, trim_pattern="_*")
    arginfo.check_value(data, strict=strict_check)
    return data


def convert_training_data_to_hdf5(input_files: List[str], h5_file: str):
    """Convert training data to HDF5 format and update the input files.

    Parameters
    ----------
    input_files : list of str
        DeePMD-kit input file names
    h5_file : str
        HDF5 file name
    """
    systems = []
    h5_dir = Path(h5_file).parent.absolute()
    cwd = Path.cwd().absolute()
    for ii in input_files:
        ii = Path(ii)
        dd = ii.parent.absolute()
        with open(ii, "r+") as f:
            jinput = json.load(f)
            if "training_data" in jinput["training"]:
                # v2.0
                p_sys = jinput["training"]["training_data"]["systems"]
            else:
                # v1.x
                p_sys = jinput["training"]["systems"]
            for ii, pp in enumerate(p_sys):
                if "#" in pp:
                    # HDF5 file
                    p1, p2 = pp.split("#")
                    ff = os.path.normpath(str((dd / p1).absolute().relative_to(cwd)))
                    pp = ff + "#" + p2
                    new_pp = os.path.normpath(os.path.relpath(ff, h5_dir)) + "/" + p2
                else:
                    pp = os.path.normpath(str((dd / pp).absolute().relative_to(cwd)))
                    new_pp = os.path.normpath(os.path.relpath(pp, h5_dir))
                p_sys[ii] = (
                    os.path.normpath(os.path.relpath(h5_file, dd)) + "#/" + str(new_pp)
                )
                systems.append(pp)
            f.seek(0)
            json.dump(jinput, f, indent=4)
    systems = list(set(systems))

    dlog.info("Combining %d training systems to %s...", len(systems), h5_file)

    with h5py.File(h5_file, "w") as f:
        for ii in systems:
            if "#" in ii:
                p1, p2 = ii.split("#")
                p1 = os.path.normpath(os.path.relpath(p1, h5_dir))
                group = f.create_group(str(p1) + "/" + p2)
                s = dpdata.LabeledSystem(ii, fmt="deepmd/hdf5")
                s.to("deepmd/hdf5", group)
            else:
                pp = os.path.normpath(os.path.relpath(ii, h5_dir))
                group = f.create_group(str(pp))
                s = dpdata.LabeledSystem(ii, fmt="deepmd/npy")
                s.to("deepmd/hdf5", group)

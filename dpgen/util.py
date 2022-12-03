#!/usr/bin/env python
# coding: utf-8
from typing import Union, List
from pathlib import Path

import h5py
from dargs import Argument

from dpgen import dlog

"""
some common utilities for generator, auto_test and data
"""

# constants define
MaxLength=70

def sepline(ch='-',sp='-',screen=False):
    r'''
    seperate the output by '-'
    '''
    if screen:
       print(ch.center(MaxLength,sp))
    else:
       dlog.info(ch.center(MaxLength,sp))

def box_center(ch='',fill=' ',sp="|"):
    r'''
    put the string at the center of |  |
    '''
    strs=ch.center(Len,fill)
    dlog.info(sp+strs[1:len(strs)-1:]+sp)


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
        with h5py.File(root_dir, 'r') as f:
            # list of keys in the h5 file
            f_keys = ["/"]
            f.visit(lambda x: f_keys.append("/" + x))
        matches = ["%s#%s"%(root_dir, d) for d in f_keys if str(Path(d) / "type.raw") in f_keys]
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

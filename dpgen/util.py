#!/usr/bin/env python
# coding: utf-8
from typing import Union, List
from pathlib import Path

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
    matches = [str(d) for d in root_dir.rglob("*") if (d / "type.raw").is_file()]
    if (root_dir / "type.raw").is_file():
        matches.append(str(root_dir))
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

#!/usr/bin/env python
# coding: utf-8


def convert_mdata(mdata, task_types=["train", "model_devi", "fp"]):
    """
    Convert mdata for DP-GEN main process.
    New convension is like mdata["fp"]["machine"],
    DP-GEN needs mdata["fp_machine"]

    Notice that we deprecate the function which can automatically select one most avalaible machine,
    since this function was only used by Angus, and only supports for Slurm.
    In the future this can be implemented.

    Parameters
    ----------
    mdata : dict
        Machine parameters to be converted.
    task_types : list of string
        Type of tasks, default is ["train", "model_devi", "fp"]

    Returns
    -------
    dict
        mdata converted
    """
    for task_type in task_types:
        if task_type in mdata:
            if isinstance(mdata[task_type], dict):
                task_data = mdata[task_type]
            elif isinstance(mdata[task_type], (list, tuple)):
                task_data = mdata[task_type][0]
            else:
                raise TypeError("mdata/%s should be dict or list!" % task_type)
            for key, item in task_data.items():
                if "comments" not in key:
                    mdata[task_type + "_" + key] = item
            group_size = task_data["resources"].get("group_size", 1)
            if group_size == 1:
                group_size = task_data.get("group_size", 1)
            mdata[task_type + "_" + "group_size"] = group_size
    return mdata

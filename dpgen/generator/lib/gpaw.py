"""Move all the GPAW related functions to here.
- functions for `arginfo.py`
- functions for `run.py`
"""


import os
import glob
from pathlib import Path

from dargs import Argument
from dpgen.generator.lib.utils import make_iter_name
from dpgen.util import set_directory
import dpdata
import numpy as np

from ..run import fp_name
# fp_name = "02.fp"


### ANCHOR functions for `arginfo.py`
def fp_style_gpaw_args() -> list[Argument]:
    args = [Argument("fp_gpaw_runfile",
                     str,
                     optional=True,
                     default="gpaw_singlepoint.py",
                     doc="Input file to run GPAW.",
                     )
            ]
    return args


### ANCHOR functions for `run.py`
def make_fp_gpaw(iter_index, jdata):
    """Make input file for customized FP style.

    Parameters
    ----------
    iter_index : int
        iter index
    jdata : dict
        Run parameters.
    """
    ## create symbolic link of the gpaw input file in the task directory
    work_path = os.path.join(make_iter_name(iter_index), fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, "task.*"))
    gpaw_runfile = jdata["fp_gpaw_runfile"]
    gpaw_runfile_source = Path(gpaw_runfile).resolve()
    assert os.path.exists(
        gpaw_runfile_source
    ), f"Can not find gpaw runfile {gpaw_runfile_source}"
    for ii in fp_tasks:
        with set_directory(Path(ii)):
            # create file `gpaw_runfile` in the current directory and symlink it to the source file
            Path(gpaw_runfile).symlink_to(gpaw_runfile_source)


def post_fp_gpaw(iter_index, jdata):
    """Post fp for custom fp. Collect data from user-defined `output_fn`.

    Parameters
    ----------
    iter_index : int
        The index of the current iteration.
    jdata : dict
        The parameter data.
    """
    model_devi_jobs = jdata["model_devi_jobs"]
    assert iter_index < len(model_devi_jobs)

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, fp_name)
    fp_tasks = glob.glob(os.path.join(work_path, "task.*"))
    fp_tasks.sort()
    if len(fp_tasks) == 0:
        return

    system_index = []
    for ii in fp_tasks:
        system_index.append(os.path.basename(ii).split(".")[1])
    system_index.sort()
    set_tmp = set(system_index)
    system_index = list(set_tmp)
    system_index.sort()

    output_fn = "CONF_ASE.traj"
    output_fmt = "ase/traj"

    for ss in system_index:
        sys_output = glob.glob(os.path.join(work_path, f"task.{ss}.*"))
        sys_output.sort()
        all_sys = dpdata.MultiSystems(type_map=jdata["type_map"])
        for oo in sys_output:
            if os.path.exists(os.path.join(oo, output_fn)):
                sys = dpdata.LabeledSystem(os.path.join(oo, output_fn), fmt=output_fmt)
                all_sys.append(sys)
        sys_data_path = os.path.join(work_path, f"data.{ss}")
        all_sys.to_deepmd_raw(sys_data_path)
        all_sys.to_deepmd_npy(sys_data_path, set_size=len(sys_output), prec=np.float64)

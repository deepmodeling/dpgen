
import os
import shutil
import glob
from packaging.version import Version
from dpgen.generator.lib.utils import symlink_user_forward_files
from dpgen.dispatcher.Dispatcher import make_submission

global_dirname_02 = "00.place_ele"
global_dirname_03 = "01.scale_pert"
global_dirname_04 = "02.md"


##### ANCHOR: Stage 1 - Geometry Optimization/ relaxation
def make_gpaw_relax(jdata, mdata):
    out_dir = jdata["out_dir"]
    cwd = os.getcwd()
    work_dir = os.path.join(out_dir, global_dirname_02)
    assert os.path.isdir(work_dir)
    work_dir = os.path.abspath(work_dir)

    gpaw_file = os.path.join(work_dir, "gpaw_runfile.py")
    shutil.copy2(jdata["relax_incar"], gpaw_file)

    ### Generate symlinks for GPAW input files
    os.chdir(work_dir)
    sys_list = glob.glob("sys-*")
    for ss in sys_list:
        os.chdir(ss)
        ln_src = os.path.relpath(gpaw_file)
        try:
            os.symlink(ln_src, "gpaw_runfile.py")
        except FileExistsError:
            pass
        os.chdir(work_dir)

    os.chdir(cwd)
    symlink_user_forward_files(
        mdata=mdata,
        task_type="fp",
        work_path=os.path.join(os.path.basename(out_dir), global_dirname_02),
        task_format={"fp": "sys-*"},
    )
    return


def run_gpaw_relax(jdata, mdata):
    fp_command = mdata["fp_command"].strip() + " gpaw_runfile.py"
    fp_group_size = mdata["fp_group_size"]
    # machine_type = mdata['fp_machine']['machine_type']
    work_dir = os.path.join(jdata["out_dir"], global_dirname_02)

    forward_files = ["POSCAR", "gpaw_runfile.py"]
    user_forward_files = mdata.get("fp" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files = ["conf_ase.traj", "calc.txt", "run.log"]
    backward_files += mdata.get("fp" + "_user_backward_files", [])
    forward_common_files = []

    relax_tasks = glob.glob(os.path.join(work_dir, "sys-*"))
    relax_tasks.sort()
    # dlog.info("work_dir",work_dir)
    # dlog.info("relax_tasks",relax_tasks)
    if len(relax_tasks) == 0:
        return

    relax_run_tasks = relax_tasks
    run_tasks = [os.path.basename(ii) for ii in relax_run_tasks]

    api_version = mdata.get("api_version", "1.0")
    if Version(api_version) < Version("1.0"):
        raise RuntimeError(
            f"API version {api_version} has been removed. Please upgrade to 1.0."
        )

    elif Version(api_version) >= Version("1.0"):
        submission = make_submission(
            mdata["fp_machine"],
            mdata["fp_resources"],
            commands=[fp_command],
            work_path=work_dir,
            run_tasks=run_tasks,
            group_size=fp_group_size,
            forward_common_files=forward_common_files,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog="fp.log",
            errlog="fp.log",
        )
        submission.run_submission()
    return


##### ANCHOR: Stage 2 - scale and perturb
def make_scale_gpaw(jdata):

    return


##### ANCHOR: Stage 3 - run AIMD
def make_gpaw_md(jdata, mdata):

    print("not implemented yet")
    return

def run_gpaw_md(jdata, mdata):

    print("not implemented yet")
    return

##### ANCHOR: Stage 4 - collect data
def coll_gpaw_md(jdata):

    print("not implemented yet")
    return
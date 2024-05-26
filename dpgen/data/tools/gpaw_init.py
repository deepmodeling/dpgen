import os
import shutil
import glob

import dpdata
from packaging.version import Version
from dpgen.generator.lib.utils import symlink_user_forward_files
from dpgen.dispatcher.Dispatcher import make_submission

### use from...import... may cause circular import
import dpgen.data.gen.create_path as create_path

from ase.io import Trajectory
from ase.io.vasp import write_vasp

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

    gpaw_file = os.path.join(work_dir, "gpaw_relax.py")
    shutil.copy2(jdata["relax_incar"], gpaw_file)

    ### Generate symlinks for GPAW input files
    os.chdir(work_dir)
    sys_list = glob.glob("sys-*")
    for ss in sys_list:
        os.chdir(ss)
        ln_src = os.path.relpath(gpaw_file)
        try:
            os.symlink(ln_src, "gpaw_relax.py")
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
    check_gpaw_input(jdata["relax_incar"])
    fp_command = mdata["fp_command"].strip() + " gpaw_relax.py"
    fp_group_size = mdata["fp_group_size"]
    # machine_type = mdata['fp_machine']['machine_type']
    work_dir = os.path.join(jdata["out_dir"], global_dirname_02)

    forward_files = ["POSCAR", "gpaw_relax.py"]
    user_forward_files = mdata.get("fp" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files = ["CONF_ASE.traj", "calc.txt", "fp.log"]
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

    ### Convert `CONF_ASE.traj` to `CONTCAR` to be used in the next step
    for ii in relax_tasks:
        os.chdir(ii)
        if os.path.isfile("CONF_ASE.traj"):
            traj = Trajectory('CONF_ASE.traj')
            write_vasp('CONTCAR', traj[-1])
        os.chdir(work_dir)
    return


##### ANCHOR: Stage 2 - scale and perturb
# Use the same `make_scale(jdata)` function as VASP


##### ANCHOR: Stage 3 - run AIMD
def make_gpaw_md(jdata, mdata):
    out_dir = jdata["out_dir"]
    scale = jdata["scale"]
    pert_numb = jdata["pert_numb"]

    cwd = os.getcwd()
    path_ps = os.path.join(out_dir, global_dirname_03)
    path_ps = os.path.abspath(path_ps)
    assert os.path.isdir(path_ps)
    os.chdir(path_ps)
    sys_ps = glob.glob("sys-*")
    sys_ps.sort()
    os.chdir(cwd)
    path_md = os.path.join(out_dir, global_dirname_04)
    path_md = os.path.abspath(path_md)
    create_path(path_md)

    gpaw_file = os.path.join(path_md, "gpaw_aimd.py")
    shutil.copy2(jdata["md_incar"], os.path.join(path_md, gpaw_file))

    for ii in sys_ps:
        for jj in scale:
            for kk in range(pert_numb + 1):
                path_work = path_md
                path_work = os.path.join(path_work, ii)
                path_work = os.path.join(path_work, f"scale-{jj:.3f}")
                path_work = os.path.join(path_work, "%06d" % kk)
                create_path(path_work)
                os.chdir(path_work)
                path_pos = path_ps
                path_pos = os.path.join(path_pos, ii)
                path_pos = os.path.join(path_pos, f"scale-{jj:.3f}")
                path_pos = os.path.join(path_pos, "%06d" % kk)
                init_pos = os.path.join(path_pos, "POSCAR")
                shutil.copy2(init_pos, "POSCAR")
                try:
                    os.symlink(os.path.relpath(gpaw_file), "gpaw_aimd.py")
                except FileExistsError:
                    pass
                os.chdir(cwd)

    symlink_user_forward_files(
        mdata=mdata,
        task_type="fp",
        work_path=os.path.join(os.path.basename(out_dir), global_dirname_04),
        task_format={"fp": "sys-*/scale*/00*"},
    )
    return


def run_gpaw_md(jdata, mdata):
    check_gpaw_input(jdata["md_incar"])
    fp_command = mdata["fp_command"].strip() + " gpaw_aimd.py"
    fp_group_size = mdata["fp_group_size"]
    # machine_type = mdata['fp_machine']['machine_type']
    work_dir = os.path.join(jdata["out_dir"], global_dirname_04)

    forward_files = ["POSCAR", "gpaw_aimd.py"]
    user_forward_files = mdata.get("fp" + "_user_forward_files", [])
    forward_files += [os.path.basename(file) for file in user_forward_files]
    backward_files = ["CONF_ASE.traj", "calc.txt", "fp.log"]
    backward_files += mdata.get("fp" + "_user_backward_files", [])
    forward_common_files = []

    path_md = work_dir
    path_md = os.path.abspath(path_md)
    cwd = os.getcwd()
    assert os.path.isdir(path_md), "md path should exists"
    md_tasks = glob.glob(os.path.join(work_dir, "sys-*/scale*/00*"))
    md_tasks.sort()

    if len(md_tasks) == 0:
        return

    md_run_tasks = md_tasks
    run_tasks = [ii.replace(work_dir + "/", "") for ii in md_run_tasks]
    # dlog.info("md_work_dir", work_dir)
    # dlog.info("run_tasks",run_tasks)
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

##### ANCHOR: Stage 4 - collect data


def coll_gpaw_md(jdata):
    out_dir = jdata["out_dir"]
    md_nstep = jdata["md_nstep"]
    scale = jdata["scale"]
    pert_numb = jdata["pert_numb"]
    coll_ndata = jdata["coll_ndata"]

    cwd = os.getcwd()
    path_md = os.path.join(out_dir, global_dirname_04)
    path_md = os.path.abspath(path_md)
    assert os.path.isdir(path_md), "md path should exists"
    os.chdir(path_md)
    sys_md = glob.glob("sys-*")
    sys_md.sort()

    for ii in sys_md:
        os.chdir(ii)
        # convert outcars
        valid_outcars = []
        for jj in scale:
            for kk in range(pert_numb):
                path_work = os.path.join(f"scale-{jj:.3f}", "%06d" % kk)
                outcar = os.path.join(path_work, "OUTCAR")
                # dlog.info("OUTCAR",outcar)
                if os.path.isfile(outcar):
                    # dlog.info("*"*40)
                    with open(outcar) as fin:
                        nforce = fin.read().count("TOTAL-FORCE")
                    # dlog.info("nforce is", nforce)
                    # dlog.info("md_nstep", md_nstep)
                    if nforce == md_nstep:
                        valid_outcars.append(outcar)
                    elif md_nstep == 0 and nforce == 1:
                        valid_outcars.append(outcar)
                    else:
                        dlog.info(
                            f"WARNING : in directory {os.getcwd()} nforce in OUTCAR is not equal to settings in INCAR"
                        )
        arg_cvt = " "
        if len(valid_outcars) == 0:
            raise RuntimeError(
                f"MD dir: {path_md}: find no valid outcar in sys {ii}, "
                "check if your vasp md simulation is correctly done"
            )

        flag = True
        if ("type_map" in jdata) and isinstance(jdata["type_map"], list):
            type_map = jdata["type_map"]
        else:
            type_map = None
        for oo in valid_outcars:
            if flag:
                _sys = dpdata.LabeledSystem(oo, type_map=type_map)
                if len(_sys) > 0:
                    all_sys = _sys
                    flag = False
                else:
                    pass
            else:
                _sys = dpdata.LabeledSystem(oo, type_map=type_map)
                if len(_sys) > 0:
                    all_sys.append(_sys)
        # create deepmd data
        if all_sys.get_nframes() >= coll_ndata:
            all_sys = all_sys.sub_system(np.arange(coll_ndata))
        all_sys.to_deepmd_raw("deepmd")
        all_sys.to_deepmd_npy("deepmd", set_size=all_sys.get_nframes())
        os.chdir(path_md)
    os.chdir(cwd)

    print(" not implemented yet")

    return


##### ANCHOR: Support functions
def check_gpaw_input(input_file: str) -> None:
    """
    Check the input files for the GPAW calculation
    """
    with open(input_file, "r") as f:
        text = f.read()

    if "calc.txt" not in text:
        raise ValueError(
            f"The GPAW calculator in file {input_file} did not contain field: txt='calc.txt'. It should be set for backward files."
        )

    if "CONF_ASE.traj" not in text:
        raise ValueError(
            f"The GPAW input file {input_file} did not output the trajectory file 'CONF_ASE.traj'. It should be set for backward files."
        )
    return

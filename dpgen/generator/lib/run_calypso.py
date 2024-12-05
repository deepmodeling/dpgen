"""calypso as model devi engine:
1. gen_structures
2. analysis
3. model devi.
"""

import glob
import os
import random
import shutil
import sys
from itertools import combinations
from pathlib import Path

import dpdata
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp
from packaging.version import Version

from dpgen import dlog
from dpgen.dispatcher.Dispatcher import make_submission
from dpgen.generator.lib.parse_calypso import _parse_calypso_input
from dpgen.generator.lib.utils import create_path, make_iter_name

train_name = "00.train"
model_devi_name = "01.model_devi"
fp_name = "02.fp"
calypso_run_opt_name = "gen_stru_analy"
calypso_model_devi_name = "model_devi_results"


def gen_structures(
    iter_index, jdata, mdata, caly_run_path, current_idx, length_of_caly_runopt_list
):
    # run calypso
    # vsc means generate elemental, binary and ternary at the same time
    vsc = jdata.get("vsc", False)  # take CALYPSO as confs generator

    model_devi_group_size = mdata["model_devi_group_size"]
    model_devi_resources = mdata["model_devi_resources"]
    api_version = mdata.get("api_version", "1.0")

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert os.path.isdir(work_path)

    calypso_run_opt_path = caly_run_path
    calypso_model_devi_path = os.path.join(work_path, calypso_model_devi_name)

    calypso_path = mdata.get("model_devi_calypso_path")
    # calypso_input_path = jdata.get('calypso_input_path')

    all_models = glob.glob(os.path.join(calypso_run_opt_path, "graph*pb"))
    model_names = [os.path.basename(ii) for ii in all_models]

    deepmdkit_python = mdata.get("model_devi_deepmdkit_python")
    command = (
        f"{deepmdkit_python} calypso_run_opt.py  1>> model_devi.log 2>> model_devi.log"
    )
    # command = "%s calypso_run_opt.py %s 1>> model_devi.log 2>> model_devi.log" % (deepmdkit_python,os.path.abspath(calypso_run_opt_path))
    # command += "  ||  %s check_outcar.py %s " % (deepmdkit_python,os.path.abspath(calypso_run_opt_path))
    command += f"  ||  {deepmdkit_python} check_outcar.py  "
    commands = [command]

    cwd = os.getcwd()
    os.chdir(calypso_run_opt_path)

    forward_files = ["POSCAR", "calypso_run_opt.py", "check_outcar.py", "input.dat"]
    backward_files = ["OUTCAR", "CONTCAR", "traj.traj", "model_devi.log"]

    run_calypso = calypso_path + "/calypso.x | tee log"
    if not vsc:
        Lpickup = _parse_calypso_input("PickUp", ".")
        PickUpStep = _parse_calypso_input("PickStep", ".")
        if os.path.exists(f"tag_pickup_{str(PickUpStep)}"):
            dlog.info(f"caution! tag_pickup_{str(PickUpStep)} exists!")
            Lpickup = "F"
        if Lpickup == "T":
            ftag = open(f"tag_pickup_{str(PickUpStep)}", "w")
            ftag.close()
            os.remove("step")
            fstep = open("step", "w")
            fstep.write("%12s" % str(PickUpStep))  # noqa: UP031
            fstep.close()
        else:
            PickUpStep = 1
            try:
                os.mkdir("opt")
            except Exception:
                pass

        popsize = int(_parse_calypso_input("PopSize", "."))
        maxstep = int(_parse_calypso_input("MaxStep", "."))

        for ii in range(int(PickUpStep) - 1, maxstep + 1):
            dlog.info(f"CALYPSO step {ii}")
            if ii == maxstep:
                os.system(f"{run_calypso}")
                break
            # run calypso

            os.system(f"{run_calypso}")

            for pop in range(ii * int(popsize), (ii + 1) * int(popsize)):
                try:
                    os.mkdir("task.%03d" % pop)  # noqa: UP031
                except Exception:
                    shutil.rmtree("task.%03d" % pop)  # noqa: UP031
                    os.mkdir("task.%03d" % pop)  # noqa: UP031
                shutil.copyfile(
                    "calypso_run_opt.py",
                    os.path.join("task.%03d" % pop, "calypso_run_opt.py"),  # noqa: UP031
                )
                shutil.copyfile(
                    "check_outcar.py",
                    os.path.join("task.%03d" % pop, "check_outcar.py"),  # noqa: UP031
                )
                shutil.copyfile(
                    f"POSCAR_{str(pop - ii * int(popsize) + 1)}",
                    os.path.join("task.%03d" % (pop), "POSCAR"),  # noqa: UP031
                )
                shutil.copyfile(
                    "input.dat",
                    os.path.join("task.%03d" % pop, "input.dat"),  # noqa: UP031
                )
            # for iii in range(1,popsize+1):
            #    shutil.copyfile('POSCAR_%s'%str(iii),os.path.join('task.%03d'%(iii-1),'POSCAR'))

            all_task = glob.glob("task.*")
            all_task.sort()

            run_tasks_ = all_task

            run_tasks = [os.path.basename(ii) for ii in run_tasks_]

            if Version(api_version) < Version("1.0"):
                raise RuntimeError(
                    f"API version {api_version} has been removed. Please upgrade to 1.0."
                )
            elif Version(api_version) >= Version("1.0"):
                os.chdir(cwd)
                submission = make_submission(
                    mdata["model_devi_machine"],
                    mdata["model_devi_resources"],
                    commands=commands,
                    work_path=calypso_run_opt_path,
                    run_tasks=run_tasks,
                    group_size=model_devi_group_size,
                    forward_common_files=model_names,
                    forward_files=forward_files,
                    backward_files=backward_files,
                    outlog="model_devi.log",
                    errlog="model_devi.log",
                )
                submission.run_submission()
                os.chdir(calypso_run_opt_path)

            sstep = os.path.join("opt", str(ii))
            os.mkdir(sstep)
            if not os.path.exists("traj"):
                os.mkdir("traj")

            for jjj in range(ii * int(popsize), (ii + 1) * int(popsize)):
                # to opt directory
                shutil.copyfile(
                    f"POSCAR_{str(jjj + 1 - ii * int(popsize))}",
                    os.path.join(sstep, f"POSCAR_{str(jjj + 1 - ii * int(popsize))}"),
                )
                shutil.copyfile(
                    os.path.join("task.%03d" % (jjj), "OUTCAR"),  # noqa: UP031
                    os.path.join(sstep, f"OUTCAR_{str(jjj + 1 - ii * int(popsize))}"),
                )
                shutil.copyfile(
                    os.path.join("task.%03d" % (jjj), "CONTCAR"),  # noqa: UP031
                    os.path.join(sstep, f"CONTCAR_{str(jjj + 1 - ii * int(popsize))}"),
                )
                # to run calypso directory
                shutil.copyfile(
                    os.path.join("task.%03d" % (jjj), "OUTCAR"),  # noqa: UP031
                    f"OUTCAR_{str(jjj + 1 - ii * int(popsize))}",
                )
                shutil.copyfile(
                    os.path.join("task.%03d" % (jjj), "CONTCAR"),  # noqa: UP031
                    f"CONTCAR_{str(jjj + 1 - ii * int(popsize))}",
                )
                # to traj
                shutil.copyfile(
                    os.path.join("task.%03d" % (jjj), "traj.traj"),  # noqa: UP031
                    os.path.join("traj", f"{str(jjj + 1)}.traj"),
                )

            tlist = glob.glob("task.*")
            for t in tlist:
                shutil.rmtree(t)

    else:
        # --------------------------------------------------------------
        type_map = jdata["type_map"]
        how_many_spec = len(type_map)
        if how_many_spec == 1:
            dlog.info("vsc mode can not work in one-element situation")
            sys.exit()

        comp_temp = list(map(list, list(combinations(type_map, 1))))
        for hms in range(2, how_many_spec + 1):
            # comp_temp = [['Mg'],['Al'],['Cu'],['Mg','Al'],['Mg','Cu'],['Al','Cu'],['Mg','Al','Cu']]
            comp_temp.extend(list(map(list, list(combinations(type_map, hms)))))

        component = []
        for comp_temp_ in comp_temp:
            component.append(
                "".join(comp_temp_)
            )  # component = ['Mg','Al','Cu','MgAl','MgCu','AlCu','MgAlCu']

        dlog.info(component)
        # calypso_input_path = jdata.get('calypso_input_path')

        pwd = os.getcwd()
        if len(glob.glob(f"input.dat.{component[0]}.*")) != 0:
            os.system("for i in input.dat.*;do mv $i ${i%.*};done")
        for idx, com in enumerate(component):
            if not os.path.exists(com):
                os.mkdir(com)
            # shutil.copyfile(os.path.join(calypso_input_path,'input.dat.%s'%com),os.path.join(com,'input.dat'))
            shutil.copyfile(f"input.dat.{com}", os.path.join(com, "input.dat"))
            os.chdir(com)
            os.system(run_calypso)
            os.chdir(pwd)

        shutil.copyfile(f"input.dat.{component[-1]}", "input.dat")

        name_list = Path(".").glob("*/POSCAR_*")
        for idx, name in enumerate(name_list):
            shutil.copyfile(name, "POSCAR_%s" % (idx + 1))
            try:
                os.mkdir("task.%04d" % (idx + 1))  # noqa: UP031
            except Exception:
                shutil.rmtree("task.%04d" % (idx + 1))  # noqa: UP031
                os.mkdir("task.%04d" % (idx + 1))  # noqa: UP031
            shutil.copyfile(
                "calypso_run_opt.py",
                os.path.join("task.%04d" % (idx + 1), "calypso_run_opt.py"),  # noqa: UP031
            )
            shutil.copyfile(
                "check_outcar.py",
                os.path.join("task.%04d" % (idx + 1), "check_outcar.py"),  # noqa: UP031
            )
            shutil.copyfile(
                f"POSCAR_{str(idx + 1)}",
                os.path.join("task.%04d" % (idx + 1), "POSCAR"),  # noqa: UP031
            )
            shutil.copyfile(
                "input.dat",
                os.path.join("task.%04d" % (idx + 1), "input.dat"),  # noqa: UP031
            )

        # sys.exit()

        all_task = glob.glob("task.*")
        all_task.sort()

        run_tasks_ = all_task

        run_tasks = [os.path.basename(ii) for ii in run_tasks_]

        if Version(api_version) < Version("1.0"):
            raise RuntimeError(
                f"API version {api_version} has been removed. Please upgrade to 1.0."
            )
        elif Version(api_version) >= Version("1.0"):
            os.chdir(cwd)
            submission = make_submission(
                mdata["model_devi_machine"],
                mdata["model_devi_resources"],
                commands=commands,
                work_path=calypso_run_opt_path,
                run_tasks=run_tasks,
                group_size=model_devi_group_size,
                forward_common_files=model_names,
                forward_files=forward_files,
                backward_files=backward_files,
                outlog="model_devi.log",
                errlog="model_devi.log",
            )
            submission.run_submission()
            os.chdir(calypso_run_opt_path)

        os.mkdir("opt")
        if not os.path.exists("traj"):
            os.mkdir("traj")
        for jjj in range(len(all_task)):
            # to opt directory
            shutil.copyfile(
                f"POSCAR_{str(jjj + 1)}",
                os.path.join("opt", f"POSCAR_{str(jjj + 1)}"),
            )
            shutil.copyfile(
                os.path.join("task.%04d" % (jjj + 1), "OUTCAR"),  # noqa: UP031
                os.path.join("opt", f"OUTCAR_{str(jjj + 1)}"),
            )
            shutil.copyfile(
                os.path.join("task.%04d" % (jjj + 1), "CONTCAR"),  # noqa: UP031
                os.path.join("opt", f"CONTCAR_{str(jjj + 1)}"),
            )
            # to run calypso directory
            shutil.copyfile(
                os.path.join("task.%04d" % (jjj + 1), "OUTCAR"),  # noqa: UP031
                f"OUTCAR_{str(jjj + 1)}",
            )
            shutil.copyfile(
                os.path.join("task.%04d" % (jjj + 1), "CONTCAR"),  # noqa: UP031
                f"CONTCAR_{str(jjj + 1)}",
            )
            # to traj
            shutil.copyfile(
                os.path.join("task.%04d" % (jjj + 1), "traj.traj"),  # noqa: UP031
                os.path.join("traj", f"{str(jjj + 1)}.traj"),
            )

        tlist = glob.glob("task.*")
        for t in tlist:
            shutil.rmtree(t)
        # --------------------------------------------------------------

    if current_idx < length_of_caly_runopt_list - 1:
        tobewrite = f"1 {str(current_idx + 1)}\n"
    elif current_idx == length_of_caly_runopt_list - 1:
        tobewrite = "2\n"

    os.chdir(cwd)
    os.chdir(work_path)
    f = open("record.calypso", "a+")
    f.write(tobewrite)
    f.close()
    os.chdir(cwd)


def gen_main(iter_index, jdata, mdata, caly_run_opt_list, gen_idx):
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)

    current_gen_path = os.path.join(
        work_path,
        "%s.%03d" % (calypso_run_opt_name, int(gen_idx)),  # noqa: UP031
    )
    if current_gen_path not in caly_run_opt_list:
        dlog.info(
            f"current gen path {current_gen_path} not in caly run opt list {caly_run_opt_list}"
        )
        sys.exit()

    indice = caly_run_opt_list.index(current_gen_path)
    for iidx, temp_path in enumerate(caly_run_opt_list):
        if iidx >= indice:
            gen_structures(
                iter_index, jdata, mdata, temp_path, iidx, len(caly_run_opt_list)
            )


def analysis(iter_index, jdata, calypso_model_devi_path):
    # Analysis

    ms = dpdata.MultiSystems(type_map=jdata["type_map"])

    cwd = os.getcwd()
    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)

    deepmd_data_path = os.path.join(work_path, "confs", "deepmd")
    traj_pos_path = os.path.join(work_path, "confs", "traj_confs")
    create_path(deepmd_data_path)
    create_path(traj_pos_path)

    # trajs to be model devi
    # traj_path = os.path.join(calypso_run_opt_path,'traj')
    # traj_list = glob.glob(traj_path+'/*.traj')
    # 'gen_struc_analy.000/traj/*.traj' 'gen_struc_analy.001/traj/*.traj' 'gen_struc_analy.002/traj/*.traj'
    traj_list = glob.glob(f"{work_path}/*/traj/*.traj")

    # read poscar from traj file in confs/traj/*.traj
    record_traj_num = 0
    for traj_name in traj_list:
        traj_num = os.path.basename(traj_name).split(".")[0]
        press_num = traj_name.split("/")[-3].split(".")[-1]
        trajs_origin = Trajectory(traj_name)

        record_traj_num += len(trajs_origin)
        if len(trajs_origin) >= 20:
            trajs = [trajs_origin[iii] for iii in [4, 9, -10, -5, -1]]
        elif 5 <= len(trajs_origin) < 20:
            trajs = [
                trajs_origin[random.randint(1, len(trajs_origin) - 1)]
                for iii in range(4)
            ]
            trajs.append(trajs[-1])
        elif 3 <= len(trajs_origin) < 5:
            trajs = [trajs_origin[round((len(trajs_origin) - 1) / 2)]]
            trajs.append(trajs[-1])
        elif len(trajs_origin) == 2:
            trajs = [trajs_origin[0], trajs_origin[-1]]
        elif len(trajs_origin) == 1:
            trajs = [trajs_origin[0]]
        else:
            pass

        for idx, traj in enumerate(trajs):
            write_vasp(
                os.path.join(
                    traj_pos_path,
                    "%d.%03d.%03d.poscar" % (int(press_num), int(traj_num), int(idx)),  # noqa: UP031
                ),
                traj,
            )

    traj_pos_list = glob.glob(traj_pos_path + "/*.poscar")

    for npos in traj_pos_list:
        try:
            ms.append(dpdata.System(npos, type_map=jdata["type_map"]))
        except Exception as e:
            dlog.info(npos, "failed : ", e)

    if len(ms) == 0:
        dlog.info("too little confs, ")
        raise RuntimeError(
            "no confs found in Analysis part and this should not happen!"
        )

    if os.path.exists(deepmd_data_path):
        shutil.rmtree(deepmd_data_path)
    ms.to_deepmd_raw(deepmd_data_path)
    ms.to_deepmd_npy(deepmd_data_path)

    split_lists = glob.glob(os.path.join(deepmd_data_path, "*"))
    for i, split_list in enumerate(split_lists):
        strus_path = os.path.join(calypso_model_devi_path, "%03d.structures" % i)  # noqa: UP031
        if not os.path.exists(strus_path):
            shutil.copytree(split_list, strus_path)
        else:
            shutil.rmtree(strus_path)
            shutil.copytree(split_list, strus_path)

    os.chdir(cwd)
    os.chdir(work_path)
    f = open("record.calypso", "a+")
    f.write("3\n")
    f.close()
    os.chdir(cwd)


def run_calypso_model_devi(iter_index, jdata, mdata):
    dlog.info("start running CALYPSO")

    iter_name = make_iter_name(iter_index)
    work_path = os.path.join(iter_name, model_devi_name)
    assert os.path.isdir(work_path)

    calypso_model_devi_path = os.path.join(work_path, calypso_model_devi_name)

    _caly_run_opt_list = glob.glob(
        os.path.join(work_path, f"{str(calypso_run_opt_name)}.*")
    )
    caly_run_opt_list = _caly_run_opt_list.copy()
    # check if gen_struc_analy.000.bk000 in caly_run_opt_list
    for temp_value in _caly_run_opt_list:
        if "bk" in temp_value:
            caly_run_opt_list.remove(temp_value)
    caly_run_opt_list.sort()

    cwd = os.getcwd()
    record_calypso_path = os.path.join(work_path, "record.calypso")
    while True:
        if not os.path.exists(record_calypso_path):
            f = open(record_calypso_path, "w")
            f.write("1 0\n")
            lines = ["1 0\n"]
            f.close()
        else:
            f = open(record_calypso_path)
            lines = f.readlines()
            f.close()

        if lines[-1].strip().strip("\n").split()[0] == "1":
            # Gen Structures
            gen_index = lines[-1].strip().strip("\n").split()[1]
            gen_main(iter_index, jdata, mdata, caly_run_opt_list, gen_index)

        elif lines[-1].strip().strip("\n") == "2":
            # Analysis & to deepmd/raw
            analysis(iter_index, jdata, calypso_model_devi_path)

        elif lines[-1].strip().strip("\n") == "3":
            # Model Devi
            _calypso_run_opt_path = os.path.abspath(caly_run_opt_list[0])
            all_models = glob.glob(os.path.join(_calypso_run_opt_path, "graph*pb"))
            cwd = os.getcwd()
            os.chdir(calypso_model_devi_path)
            args = " ".join(
                [
                    "calypso_run_model_devi.py",
                    "--all_models",
                    " ".join(all_models),
                    "--type_map",
                    " ".join(jdata.get("type_map")),
                ]
            )
            deepmdkit_python = mdata.get("model_devi_deepmdkit_python")
            os.system(f"{deepmdkit_python} {args} ")
            # Modd(iter_index,calypso_model_devi_path,all_models,jdata)
            os.chdir(cwd)

        elif lines[-1].strip().strip("\n") == "4":
            dlog.info("Model Devi is done.")
            # return
            break

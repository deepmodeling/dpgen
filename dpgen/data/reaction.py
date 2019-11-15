"""
input: trajectory
00: build dataset (mddatasetbuilder)
01: fp (gaussian)
02: convert to deepmd data
output: data
"""

import argparse
import glob
import json
import os

import dpdata
from dpgen import dlog
from dpgen.dispatcher.Dispatcher import make_dispatcher
from dpgen.generator.run import create_path, make_fp_task_name
from dpgen.util import sepline

build_path = "00.build"
fp_path = "01.fp"
data_path = "02.data"

trj_path = "lammpstrj"
dataset_name = "dpgen_init"


def link_trj(jdata):
    """link lammpstrj"""
    create_path(build_path)
    task_path = os.path.join(build_path, "task.000")
    create_path(task_path)

    os.symlink(os.path.abspath(jdata["lammpstrj"]), os.path.abspath(
        os.path.join(task_path, trj_path)))


def run_build_dataset(jdata, mdata, dispatcher, log_file="log"):
    work_path = build_path
    build_command = "{cmd} -n {dataset_name} -a {type_map} -d {lammpstrj} -c {cutoff} -i {interval} -s {dataset_size} -k \"{qmkeywords}\" --nprocjob {nprocjob} --nproc {nproc}".format(
        cmd=mdata["build_command"],
        type_map=" ".join(jdata["type_map"]),
        lammpstrj=trj_path,
        cutoff=jdata["cutoff"],
        interval=jdata["interval"],
        dataset_size=jdata["dataset_size"],
        qmkeywords=jdata["qmkeywords"],
        nprocjob=mdata["fp_resources"]["task_per_node"],
        nproc=mdata["build_resources"]["task_per_node"],
        dataset_name=dataset_name
    )
    run_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    run_tasks.sort()
    run_tasks = [os.path.basename(ii) for ii in run_tasks]

    dispatcher.run_jobs(mdata['build_resources'],
                        [build_command],
                        work_path,
                        run_tasks,
                        1,
                        [],
                        [trj_path],
                        [f"{dataset_name}_gjf"],
                        outlog=log_file,
                        errlog=log_file)


def link_fp_input():
    all_input_file = glob.glob(os.path.join(
        build_path, "task.*", f"{dataset_name}_gjf", "*", "*.gjf"))
    work_path = fp_path
    create_path(work_path)

    for ii, fin in enumerate(all_input_file):
        dst_path = os.path.join(work_path, make_fp_task_name(0, ii))
        create_path(dst_path)
        os.symlink(os.path.abspath(fin), os.path.abspath(
            os.path.join(dst_path, "input")))


def run_fp(jdata,
           mdata,
           dispatcher,
           log_file="log",
           forward_common_files=[]):
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    work_path = fp_path

    fp_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    fp_tasks.sort()
    if len(fp_tasks) == 0:
        return

    fp_run_tasks = fp_tasks

    run_tasks = [os.path.basename(ii) for ii in fp_run_tasks]

    dispatcher.run_jobs(mdata['fp_resources'],
                        [fp_command],
                        work_path,
                        run_tasks,
                        fp_group_size,
                        [],
                        ["input"],
                        ["output"],
                        outlog=log_file,
                        errlog=log_file)


def convert_data(jdata):
    s = dpdata.MultiSystems(*[dpdata.LabeledSystem(x, fmt="gaussian/log")
                              for x in glob.glob(os.path.join(fp_path, "*", "output"))],
                            type_map=jdata["type_map"])
    s.to_deepmd_npy(data_path)


def gen_init_reaction(args):
    try:
        import ruamel
        from monty.serialization import loadfn, dumpfn
        warnings.simplefilter(
            'ignore', ruamel.yaml.error.MantissaNoDotYAML1_1Warning)
        jdata = loadfn(args.PARAM)
        if args.MACHINE is not None:
            mdata = loadfn(args.MACHINE)
    except:
        with open(args.PARAM, 'r') as fp:
            jdata = json.load(fp)
        if args.MACHINE is not None:
            with open(args.MACHINE, "r") as fp:
                mdata = json.load(fp)

    record = "record.reaction"
    iter_rec = -1
    numb_task = 5
    if os.path.isfile(record):
        with open(record) as frec:
            for line in frec:
                iter_rec = int(line.strip())
        dlog.info("continue from task %02d" % iter_rec)
    for ii in range(numb_task):
        sepline(ii, '-')
        if ii <= iter_rec:
            continue
        elif ii == 0:
            dispatcher = make_dispatcher(mdata["build_machine"])
            link_trj(jdata)
        elif ii == 1:
            run_build_dataset(jdata, mdata, dispatcher)
        elif ii == 2:
            link_fp_input()
        elif ii == 3:
            dispatcher = make_dispatcher(mdata["fp_machine"])
            run_fp(jdata, mdata, dispatcher)
        elif ii == 4:
            convert_data(jdata)
        with open(record, "a") as frec:
            frec.write(ii)

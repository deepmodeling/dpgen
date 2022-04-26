"""
input: trajectory
00: ReaxFF MD (lammps)
01: build dataset (mddatasetbuilder)
02: fp (gaussian)
03: convert to deepmd data
output: data
"""

import argparse
import glob
import json
import os
import random

import dpdata
from dpgen import dlog
from dpgen.dispatcher.Dispatcher import make_dispatcher
from dpgen.generator.run import create_path, make_fp_task_name
from dpgen.util import sepline

reaxff_path = "00.reaxff"
build_path = "01.build"
fp_path = "02.fp"
data_path = "03.data"

trj_path = "lammpstrj"
ff_path = "ffield.reax"
data_init_path = "data.init"
control_path = "lmp_control"
lmp_path = "in.lmp"
dataset_name = "dpgen_init"


def link_reaxff(jdata):
    create_path(reaxff_path)
    task_path = os.path.join(reaxff_path, "task.000")
    create_path(task_path)

    rdata = jdata['reaxff']
    os.symlink(os.path.abspath(rdata["data"]), os.path.abspath(
        os.path.join(task_path, data_init_path)))
    os.symlink(os.path.abspath(rdata["ff"]), os.path.abspath(
        os.path.join(task_path, ff_path)))
    os.symlink(os.path.abspath(rdata["control"]), os.path.abspath(
        os.path.join(task_path, control_path)))
    with open(os.path.join(task_path, lmp_path), 'w') as f:
        f.write(make_lmp(jdata))


def make_lmp(jdata):
    rdata = jdata['reaxff']
    lmp_string = """units real
atom_style charge
read_data data.init
pair_style reax/c lmp_control
pair_coeff * * ffield.reax {type_map}
velocity all create {temp} {rand}
fix 1 all nvt temp {temp} {temp} {tau_t}
fix 2 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c
dump 1 all custom {dump_freq} lammpstrj id type x y z 
timestep {dt}
run	{nstep}
""".format(
        type_map=" ".join(jdata['type_map']),
        temp=rdata['temp'],
        rand=random.randrange(1000000-1)+1,
        tau_t=rdata['tau_t'],
        dump_freq=rdata['dump_freq'],
        dt=rdata['dt'],
        nstep=rdata['nstep']
    )
    return lmp_string


def run_reaxff(jdata, mdata, dispatcher, log_file="reaxff_log"):
    work_path = reaxff_path
    reaxff_command = "{} -in {}".format(mdata["reaxff_command"], lmp_path)
    run_tasks = glob.glob(os.path.join(work_path, 'task.*'))
    run_tasks.sort()
    run_tasks = [os.path.basename(ii) for ii in run_tasks]

    dispatcher.run_jobs(mdata['reaxff_resources'],
                        [reaxff_command],
                        work_path,
                        run_tasks,
                        1,
                        [],
                        [ff_path, data_init_path, control_path, lmp_path],
                        [trj_path],
                        outlog=log_file,
                        errlog=log_file)


def link_trj(jdata):
    """link lammpstrj"""
    create_path(build_path)
    task_path = os.path.join(build_path, "task.000")
    create_path(task_path)

    os.symlink(os.path.abspath(os.path.join(reaxff_path, "task.000", trj_path)), os.path.abspath(
        os.path.join(task_path, trj_path)))


def run_build_dataset(jdata, mdata, dispatcher, log_file="build_log"):
    work_path = build_path
    build_command = "{cmd} -n {dataset_name} -a {type_map} -d {lammpstrj} -c {cutoff} -s {dataset_size} -k \"{qmkeywords}\" --nprocjob {nprocjob} --nproc {nproc}".format(
        cmd=mdata["build_command"],
        type_map=" ".join(jdata["type_map"]),
        lammpstrj=trj_path,
        cutoff=jdata["cutoff"],
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
                        [f"dataset_{dataset_name}_gjf"],
                        outlog=log_file,
                        errlog=log_file)


def link_fp_input():
    all_input_file = glob.glob(os.path.join(
        build_path, "task.*", f"dataset_{dataset_name}_gjf", "*", "*.gjf"))
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
           log_file="output",
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
                        [log_file],
                        outlog=log_file,
                        errlog=log_file)


def convert_data(jdata):
    s = dpdata.MultiSystems(*[dpdata.LabeledSystem(x, fmt="gaussian/log")
                              for x in glob.glob(os.path.join(fp_path, "*", "output"))],
                            type_map=jdata["type_map"])
    s.to_deepmd_npy(data_path)
    dlog.info("Initial data is avaiable in %s" % os.path.abspath(data_path))


def gen_init_reaction(args):
    try:
        import ruamel
        from monty.serialization import loadfn, dumpfn
        warnings.simplefilter(
            'ignore', ruamel.yaml.error.MantissaNoDotYAML1_1Warning)
        jdata = loadfn(args.PARAM)
        if args.MACHINE is not None:
            mdata = loadfn(args.MACHINE)
    except Exception:
        with open(args.PARAM, 'r') as fp:
            jdata = json.load(fp)
        if args.MACHINE is not None:
            with open(args.MACHINE, "r") as fp:
                mdata = json.load(fp)

    record = "record.reaction"
    iter_rec = -1
    numb_task = 7
    if os.path.isfile(record):
        with open(record) as frec:
            for line in frec:
                iter_rec = int(line.strip())
        dlog.info("continue from task %02d" % iter_rec)
    for ii in range(numb_task):
        sepline(str(ii), '-')
        if ii <= iter_rec:
            continue
        elif ii == 0:
            link_reaxff(jdata)
        elif ii == 1:
            dispatcher = make_dispatcher(mdata["reaxff_machine"])
            run_reaxff(jdata, mdata, dispatcher)
        elif ii == 2:
            link_trj(jdata)
        elif ii == 3:
            dispatcher = make_dispatcher(mdata["build_machine"])
            run_build_dataset(jdata, mdata, dispatcher)
        elif ii == 4:
            link_fp_input()
        elif ii == 5:
            dispatcher = make_dispatcher(mdata["fp_machine"])
            run_fp(jdata, mdata, dispatcher)
        elif ii == 6:
            convert_data(jdata)
        with open(record, "a") as frec:
            frec.write(str(ii)+'\n')

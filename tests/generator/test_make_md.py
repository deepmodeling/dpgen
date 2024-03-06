import copy
import glob
import json
import os
import re
import shutil
import sys
import unittest

import dpdata
import numpy as np

from dpgen.generator.run import _read_model_devi_file, parse_cur_job_sys_revmat

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "generator"
import tempfile

from .comp_sys import test_atom_names, test_atom_types, test_cell, test_coord
from .context import (
    find_only_one_key,
    machine_file,
    machine_file_v1,
    make_model_devi,
    my_file_cmp,
    param_amber_file,
    param_file,
    param_pimd_file,
    parse_cur_job,
    parse_cur_job_revmat,
    revise_by_keys,
    revise_lmp_input_dump,
    revise_lmp_input_model,
    revise_lmp_input_plm,
    run_model_devi,
)


def _make_fake_models(idx, numb_models):
    train_dir = os.path.join("iter.%06d" % idx, "00.train")
    os.makedirs(train_dir, exist_ok=True)
    pwd = os.getcwd()
    os.chdir(train_dir)
    for ii in range(numb_models):
        os.makedirs("%03d" % ii, exist_ok=True)
        with open(os.path.join("%03d" % ii, "forzen_model.pb"), "w") as fp:
            fp.write(str(ii))
        if not os.path.isfile("graph.%03d.pb" % ii):
            os.symlink(
                os.path.join("%03d" % ii, "forzen_model.pb"), "graph.%03d.pb" % ii
            )
    os.chdir(pwd)


def _check_confs(testCase, idx, jdata):
    md_dir = os.path.join("iter.%06d" % idx, "01.model_devi")
    tasks = glob.glob(os.path.join(md_dir, "task.*"))
    tasks.sort()
    cur_job = jdata["model_devi_jobs"][idx]
    sys_idx = cur_job["sys_idx"]
    sys_configs = jdata["sys_configs"]
    poscars = []
    for ii in sys_idx:
        sys_poscars = []
        for ss in sys_configs[ii]:
            tmp_poscars = sorted(glob.glob(ss))
            sys_poscars += tmp_poscars
        poscars.append(sys_poscars)
    for ii in tasks:
        conf_file = os.path.join(ii, "conf.lmp")
        l_conf_file = os.path.basename(os.readlink(conf_file))
        poscar_file = poscars[int(l_conf_file.split(".")[0])][
            int(l_conf_file.split(".")[1])
        ]
        sys_0 = dpdata.System(conf_file, type_map=jdata["type_map"])
        sys_1 = dpdata.System(poscar_file, type_map=jdata["type_map"])
        test_atom_names(testCase, sys_0, sys_1)
        test_atom_types(testCase, sys_0, sys_1)
        test_cell(testCase, sys_0, sys_1)
        test_coord(testCase, sys_0, sys_1)


def _check_pb(testCase, idx):
    md_dir = os.path.join("iter.%06d" % idx, "01.model_devi")
    tr_dir = os.path.join("iter.%06d" % idx, "00.train")
    md_pb = glob.glob(os.path.join(md_dir, "grapb*pb"))
    tr_pb = glob.glob(os.path.join(tr_dir, "grapb*pb"))
    md_pb.sort()
    tr_pb.sort()
    for ii, jj in zip(md_pb, tr_pb):
        my_file_cmp(testCase, ii, jj)


def _check_traj_dir(testCase, idx):
    md_dir = os.path.join("iter.%06d" % idx, "01.model_devi")
    tasks = glob.glob(os.path.join(md_dir, "task.*"))
    tasks.sort()
    for ii in tasks:
        testCase.assertTrue(os.path.isdir(os.path.join(ii, "traj")))


def _get_lammps_pt(lmp_input):
    with open(lmp_input) as fp:
        for ii in fp:
            if "variable" in ii and "TEMP" in ii and "TEMP_NBEADS" not in ii:
                lt = float(ii.split()[3])
            if "variable" in ii and "PRES" in ii:
                lp = float(ii.split()[3])
    return lt, lp


def _check_pt(testCase, idx, jdata):
    md_dir = os.path.join("iter.%06d" % idx, "01.model_devi")
    tasks = glob.glob(os.path.join(md_dir, "task.*"))
    tasks.sort()
    cur_job = jdata["model_devi_jobs"][idx]
    ensemble, nsteps, trj_freq, temps, press, pka_e, dt, nbeads = parse_cur_job(cur_job)
    testCase.assertTrue(ensemble, "npt")
    # get poscars
    sys_idx = cur_job["sys_idx"]
    sys_configs = jdata["sys_configs"]
    poscars = []
    for ii in sys_idx:
        sys_poscars = []
        for ss in sys_configs[ii]:
            tmp_poscars = glob.glob(ss)
            sys_poscars += tmp_poscars
        sys_poscars.sort()
        poscars.append(sys_poscars)
    for sidx, ii in enumerate(poscars):
        count = 0
        for ss in ii:
            for tt in temps:
                for pp in press:
                    task_dir = os.path.join(md_dir, "task.%03d.%06d" % (sidx, count))
                    lt, lp = _get_lammps_pt(os.path.join(task_dir, "input.lammps"))
                    testCase.assertAlmostEqual(lt, tt)
                    testCase.assertAlmostEqual(lp, pp)
                    count += 1


class TestMakeModelDevi(unittest.TestCase):
    def tearDown(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        if os.path.isdir("test_model_devi_pimd"):
            shutil.rmtree("test_model_devi_pimd")

    def test_make_model_devi(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_file) as fp:
            jdata = json.load(fp)
        with open(machine_file) as fp:
            mdata = json.load(fp)
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        _check_confs(self, 0, jdata)
        _check_traj_dir(self, 0)
        _check_pt(self, 0, jdata)
        # shutil.rmtree('iter.000000')

    def test_make_model_devi_pimd(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_pimd_file) as fp:
            jdata = json.load(fp)
        with open(machine_file) as fp:
            mdata = json.load(fp)
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        _check_confs(self, 0, jdata)
        _check_traj_dir(self, 0)
        _check_pt(self, 0, jdata)

    def test_make_model_devi_nopbc_npt(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_file) as fp:
            jdata = json.load(fp)
            jdata["model_devi_nopbc"] = True
        with open(machine_file) as fp:
            mdata = json.load(fp)
        _make_fake_models(0, jdata["numb_models"])
        cwd = os.getcwd()
        with self.assertRaises(RuntimeError):
            make_model_devi(0, jdata, mdata)
        os.chdir(cwd)

    def test_make_model_devi_nopbc_nvt(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_file) as fp:
            jdata = json.load(fp)
            jdata["model_devi_nopbc"] = True
            jdata["model_devi_jobs"][0]["ensemble"] = "nvt"
        with open(machine_file) as fp:
            mdata = json.load(fp)
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        # _check_confs(self, 0, jdata)
        _check_traj_dir(self, 0)
        _check_pt(self, 0, jdata)
        # shutil.rmtree('iter.000000')

    def test_run_model_devi(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_file) as fp:
            jdata = json.load(fp)
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, {})
        with tempfile.TemporaryDirectory() as remote_root:
            run_model_devi(
                0,
                jdata,
                {
                    "api_version": "1.0",
                    "model_devi_command": (
                        "test -f input.lammps"
                        "&& touch model_devi.out log.lammps traj/0.lammpstrj"
                        "&& echo lmp"
                    ),
                    "model_devi_machine": {
                        "batch_type": "shell",
                        "local_root": "./",
                        "remote_root": remote_root,
                        "context_type": "local",
                    },
                    "model_devi_resources": {
                        "group_size": 1,
                    },
                    "model_devi_group_size": 1,
                },
            )

    def test_run_model_devi_pimd(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_pimd_file) as fp:
            jdata = json.load(fp)
        with open(machine_file_v1) as fp:
            mdata = json.load(fp)
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        with tempfile.TemporaryDirectory() as remote_root:
            run_model_devi(
                0,
                jdata,
                {
                    "api_version": "1.0",
                    "model_devi_command": (
                        "touch model_devi1.out model_devi2.out model_devi3.out model_devi4.out log.lammps.0 log.lammps.1 log.lammps.2 log.lammps.3"
                        "&& echo lmp"
                    ),
                    "model_devi_machine": {
                        "batch_type": "shell",
                        "local_root": "./",
                        "remote_root": remote_root,
                        "context_type": "local",
                    },
                    "model_devi_group_size": 1,
                    "model_devi_resources": {
                        "group_size": 1,
                        "cpu_per_node": 4,
                    },
                },
            )

    def test_read_model_devi_file_pimd(self):
        path = "test_model_devi_pimd"
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
        path = os.path.join(path, "iter.000000/01.model_devi/task.000.000000")
        os.makedirs(os.path.join(path, "traj"), exist_ok=True)
        for i in range(4):
            for j in range(0, 5, 2):
                with open(os.path.join(path, f"traj/{j}.lammpstrj{i+1}"), "a") as fp:
                    fp.write(f"{i} {j}\n")
        model_devi_array = np.zeros([3, 7])
        model_devi_array[:, 0] = np.array([0, 2, 4])
        model_devi_total_array = np.zeros([12, 7])
        total_steps = np.array([0, 2, 4, 5, 7, 9, 10, 12, 14, 15, 17, 19])
        model_devi_total_array[:, 0] = total_steps
        for i in range(4):
            model_devi_array[:, 4] = 0.1 * (i + 1) * np.arange(1, 4)
            model_devi_total_array[i * 3 : (i + 1) * 3, 4] = model_devi_array[:, 4]
            np.savetxt(
                os.path.join(path, f"model_devi{i+1}.out"),
                model_devi_array,
                fmt="%.12e",
            )
        _read_model_devi_file(path)
        model_devi_out = np.loadtxt(os.path.join(path, "model_devi.out"))
        traj_files = glob.glob(os.path.join(path, "traj/*.lammpstrj"))
        traj_files = sorted(
            traj_files, key=lambda x: int(re.search(r"(\d+).lammpstrj", x).group(1))
        )
        for idx, traj in enumerate(traj_files):
            traj_content = np.loadtxt(traj)
            ibead = idx // 3
            istep = (idx % 3) * 2
            np.testing.assert_array_almost_equal(traj_content, np.array([ibead, istep]))
        np.testing.assert_array_almost_equal(model_devi_out, model_devi_total_array)
        for istep in total_steps:
            self.assertTrue(
                os.path.isfile(os.path.join(path, f"traj/{istep}.lammpstrj"))
            )


class TestMakeModelDeviRevMat(unittest.TestCase):
    def tearDown(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")

    def test_make_model_devi(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        jdata = {
            "type_map": ["Mg", "Al"],
            "mass_map": [24, 27],
            "init_data_prefix": "data",
            "init_data_sys": ["deepmd"],
            "init_batch_size": [16],
            "sys_configs_prefix": os.getcwd(),
            "sys_configs": [
                ["data/al.fcc.02x02x02/01.scale_pert/sys-0032/scale*/000001/POSCAR"],
                ["data/al.fcc.02x02x02/01.scale_pert/sys-0032/scale*/000000/POSCAR"],
            ],
            "numb_models": 4,
            "shuffle_poscar": False,
            "model_devi_f_trust_lo": 0.050,
            "model_devi_f_trust_hi": 0.150,
            "model_devi_plumed": True,
            "model_devi_jobs": [
                {
                    "sys_idx": [0, 1],
                    "traj_freq": 10,
                    "template": {"lmp": "lmp/input.lammps", "plm": "lmp/input.plumed"},
                    "rev_mat": {
                        "lmp": {"V_NSTEPS": [1000], "V_TEMP": [50, 100]},
                        "plm": {"V_DIST0": [3, 4]},
                    },
                    "sys_rev_mat": {
                        "0": {"lmp": {"V_PRES": [1, 10]}, "plm": {"V_DIST1": [5, 6]}},
                        "1": {
                            "lmp": {"V_PRES": [1, 5, 10]},
                            "plm": {"V_DIST1": [5, 6, 7]},
                        },
                    },
                }
            ],
        }
        mdata = {"deepmd_version": "1"}
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        _check_confs(self, 0, jdata)
        _check_traj_dir(self, 0)
        # check the first task
        md_dir = os.path.join("iter.%06d" % 0, "01.model_devi")
        tasks = glob.glob(os.path.join(md_dir, "task.*"))
        tasks.sort()
        # each system contains 2 frames
        self.assertEqual(
            len(tasks),
            (
                len(jdata["model_devi_jobs"][0]["rev_mat"]["lmp"]["V_NSTEPS"])
                * len(jdata["model_devi_jobs"][0]["rev_mat"]["lmp"]["V_TEMP"])
                * len(jdata["model_devi_jobs"][0]["rev_mat"]["plm"]["V_DIST0"])
                * (
                    len(
                        jdata["model_devi_jobs"][0]["sys_rev_mat"]["0"]["lmp"]["V_PRES"]
                    )
                    * len(
                        jdata["model_devi_jobs"][0]["sys_rev_mat"]["0"]["plm"][
                            "V_DIST1"
                        ]
                    )
                    + len(
                        jdata["model_devi_jobs"][0]["sys_rev_mat"]["1"]["lmp"]["V_PRES"]
                    )
                    * len(
                        jdata["model_devi_jobs"][0]["sys_rev_mat"]["1"]["plm"][
                            "V_DIST1"
                        ]
                    )
                )
                * 2
            ),
        )

        cur_job = jdata["model_devi_jobs"][0]
        rev_keys = ["V_NSTEPS", "V_TEMP", "V_PRES", "V_DIST0", "V_DIST1"]
        rev_matrix = []
        # 2 systems with each 2 frames
        for i0 in cur_job["rev_mat"]["lmp"]["V_NSTEPS"]:
            for i1 in cur_job["rev_mat"]["lmp"]["V_TEMP"]:
                for i3 in cur_job["rev_mat"]["plm"]["V_DIST0"]:
                    for i2 in cur_job["sys_rev_mat"]["0"]["lmp"]["V_PRES"]:
                        for i4 in cur_job["sys_rev_mat"]["0"]["plm"]["V_DIST1"]:
                            rev_matrix.append([i0, i1, i2, i3, i4])
        for i0 in cur_job["rev_mat"]["lmp"]["V_NSTEPS"]:
            for i1 in cur_job["rev_mat"]["lmp"]["V_TEMP"]:
                for i3 in cur_job["rev_mat"]["plm"]["V_DIST0"]:
                    for i2 in cur_job["sys_rev_mat"]["0"]["lmp"]["V_PRES"]:
                        for i4 in cur_job["sys_rev_mat"]["0"]["plm"]["V_DIST1"]:
                            rev_matrix.append([i0, i1, i2, i3, i4])
        for i0 in cur_job["rev_mat"]["lmp"]["V_NSTEPS"]:
            for i1 in cur_job["rev_mat"]["lmp"]["V_TEMP"]:
                for i3 in cur_job["rev_mat"]["plm"]["V_DIST0"]:
                    for i2 in cur_job["sys_rev_mat"]["1"]["lmp"]["V_PRES"]:
                        for i4 in cur_job["sys_rev_mat"]["1"]["plm"]["V_DIST1"]:
                            rev_matrix.append([i0, i1, i2, i3, i4])
        for i0 in cur_job["rev_mat"]["lmp"]["V_NSTEPS"]:
            for i1 in cur_job["rev_mat"]["lmp"]["V_TEMP"]:
                for i3 in cur_job["rev_mat"]["plm"]["V_DIST0"]:
                    for i2 in cur_job["sys_rev_mat"]["1"]["lmp"]["V_PRES"]:
                        for i4 in cur_job["sys_rev_mat"]["1"]["plm"]["V_DIST1"]:
                            rev_matrix.append([i0, i1, i2, i3, i4])
        numb_rev = len(rev_matrix)
        for ii in range(len(tasks)):
            with open(os.path.join(tasks[ii], "job.json")) as fp:
                rev_values = rev_matrix[ii % numb_rev]
                job_recd = json.load(fp)
                for kk in job_recd.keys():
                    kidx = rev_keys.index(kk)
                    self.assertEqual(rev_values[kidx], job_recd[kk])

        cwd_ = os.getcwd()
        os.chdir(tasks[0])
        with open("input.lammps") as fp:
            lines = fp.readlines()
            for ii in lines:
                if "variable" in ii and "TEMP" in ii:
                    self.assertEqual("variable TEMP equal 50", " ".join(ii.split()))
                if "variable" in ii and "PRES" in ii:
                    self.assertEqual("variable PRES equal 1", " ".join(ii.split()))
                if "variable" in ii and "NSTEPS" in ii:
                    self.assertEqual("variable NSTEPS equal 1000", " ".join(ii.split()))
        with open("input.plumed") as fp:
            lines = fp.readlines()
            for ii in lines:
                if "RESTRAINT" in ii:
                    self.assertEqual(
                        "RESTRAINT ARG=d1,d2 AT=3,5 KAPPA=150.0,150.0 LABEL=restraint",
                        " ".join(ii.split()),
                    )
        os.chdir(cwd_)

    def test_make_model_devi_null(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        jdata = {
            "type_map": ["Mg", "Al"],
            "mass_map": [24, 27],
            "init_data_prefix": "data",
            "init_data_sys": ["deepmd"],
            "init_batch_size": [16],
            "sys_configs_prefix": os.getcwd(),
            "sys_configs": [
                ["data/al.fcc.02x02x02/01.scale_pert/sys-0032/scale*/000001/POSCAR"],
                ["data/al.fcc.02x02x02/01.scale_pert/sys-0032/scale*/000000/POSCAR"],
            ],
            "numb_models": 4,
            "shuffle_poscar": False,
            "model_devi_f_trust_lo": 0.050,
            "model_devi_f_trust_hi": 0.150,
            "model_devi_plumed": True,
            "model_devi_jobs": [
                {
                    "sys_idx": [0, 1],
                    "traj_freq": 10,
                    "template": {"lmp": "lmp/input.lammps", "plm": "lmp/input.plumed"},
                }
            ],
        }
        mdata = {"deepmd_version": "1"}
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        _check_confs(self, 0, jdata)
        _check_traj_dir(self, 0)
        # check the first task
        md_dir = os.path.join("iter.%06d" % 0, "01.model_devi")
        tasks = glob.glob(os.path.join(md_dir, "task.*"))
        # 4 accounts for 2 systems each with 2 frames
        self.assertEqual(len(tasks), (4))
        tasks.sort()
        cwd_ = os.getcwd()
        os.chdir(tasks[0])
        with open("input.lammps") as fp:
            lines = fp.readlines()
            for ii in lines:
                if "variable" in ii and "TEMP" in ii:
                    self.assertEqual("variable TEMP equal V_TEMP", " ".join(ii.split()))
                if "variable" in ii and "PRES" in ii:
                    self.assertEqual("variable PRES equal V_PRES", " ".join(ii.split()))
                if "variable" in ii and "NSTEPS" in ii:
                    self.assertEqual(
                        "variable NSTEPS equal V_NSTEPS", " ".join(ii.split())
                    )
        with open("input.plumed") as fp:
            lines = fp.readlines()
            for ii in lines:
                if "RESTRAINT" in ii:
                    self.assertEqual(
                        "RESTRAINT ARG=d1,d2 AT=V_DIST0,V_DIST1 KAPPA=150.0,150.0 LABEL=restraint",
                        " ".join(ii.split()),
                    )
        os.chdir(cwd_)


class TestParseCurJobRevMat(unittest.TestCase):
    def setUp(self):
        self.cur_job = {
            "sys_idx": [0, 1],
            "template": {"lmp": "lmp/input.lammps", "plm": "lmp/input.plumed"},
            "rev_mat": {
                "lmp": {"V_NSTEPS": [1000], "V_TEMP": [50, 100], "V_PRES": [1, 10]},
                "plm": {"V_DIST0": [3, 4], "V_DIST1": [5, 6]},
            },
        }
        self.ref_matrix = []
        for i0 in self.cur_job["rev_mat"]["lmp"]["V_NSTEPS"]:
            for i1 in self.cur_job["rev_mat"]["lmp"]["V_TEMP"]:
                for i2 in self.cur_job["rev_mat"]["lmp"]["V_PRES"]:
                    for i3 in self.cur_job["rev_mat"]["plm"]["V_DIST0"]:
                        for i4 in self.cur_job["rev_mat"]["plm"]["V_DIST1"]:
                            self.ref_matrix.append([i0, i1, i2, i3, i4])
        self.ref_keys = ["V_NSTEPS", "V_TEMP", "V_PRES", "V_DIST0", "V_DIST1"]
        self.ref_nlmp = 3

    def test_parse_cur_job(self):
        rk, rm, nl = parse_cur_job_revmat(self.cur_job, use_plm=True)
        self.assertEqual(rk, self.ref_keys)
        self.assertEqual(nl, self.ref_nlmp)
        self.assertEqual(rm, self.ref_matrix)


class TestParseCurJobSysRevMat(unittest.TestCase):
    def setUp(self):
        self.cur_job = {
            "sys_idx": [0, 1],
            "template": {"lmp": "lmp/input.lammps", "plm": "lmp/input.plumed"},
            "rev_mat": {
                "lmp": {"V_NSTEPS": [1000], "V_TEMP": [50, 100]},
                "plm": {"V_DIST0": [3, 4]},
            },
            "sys_rev_mat": {
                "0": {"lmp": {"V_PRES": [1, 10]}, "plm": {"V_DIST1": [5, 6]}},
                "1": {"lmp": {"V_PRES": [1, 10, 20]}, "plm": {"V_DIST1": [5, 6, 7]}},
            },
        }
        self.sys_ref_matrix = [[], []]
        for i0 in self.cur_job["sys_rev_mat"]["0"]["lmp"]["V_PRES"]:
            for i1 in self.cur_job["sys_rev_mat"]["0"]["plm"]["V_DIST1"]:
                self.sys_ref_matrix[0].append([i0, i1])
        for i0 in self.cur_job["sys_rev_mat"]["1"]["lmp"]["V_PRES"]:
            for i1 in self.cur_job["sys_rev_mat"]["1"]["plm"]["V_DIST1"]:
                self.sys_ref_matrix[1].append([i0, i1])
        self.sys_ref_keys = ["V_PRES", "V_DIST1"]
        self.sys_ref_nlmp_0 = 1
        self.sys_ref_nlmp_1 = 1

    def test_parse_cur_job(self):
        rk0, rm0, nl0 = parse_cur_job_sys_revmat(self.cur_job, 0, use_plm=True)
        rk1, rm1, nl1 = parse_cur_job_sys_revmat(self.cur_job, 1, use_plm=True)
        self.assertEqual(rk0, self.sys_ref_keys)
        self.assertEqual(nl0, self.sys_ref_nlmp_0)
        self.assertEqual(rm0, self.sys_ref_matrix[0])
        self.assertEqual(rk1, self.sys_ref_keys)
        self.assertEqual(nl1, self.sys_ref_nlmp_1)
        self.assertEqual(rm1, self.sys_ref_matrix[1])


class MakeModelDeviByReviseMatrix(unittest.TestCase):
    def test_find_only_one_key_1(self):
        lines = ["aaa bbb ccc\n", "bbb ccc\n", "ccc bbb ccc\n"]
        idx = find_only_one_key(lines, ["bbb", "ccc"])
        self.assertEqual(idx, 1)

    def test_find_only_one_key_0(self):
        lines = ["aaa bbb\n", "bbb aaa\n", "ccc ddd\n"]
        with self.assertRaises(RuntimeError):
            idx = find_only_one_key(lines, ["ccc", "eee"])

    def test_find_only_one_key_2(self):
        lines = ["aaa bbb\n", "bbb ccc\n", "bbb ccc\n", "fff eee\n"]
        with self.assertRaises(RuntimeError):
            idx = find_only_one_key(lines, ["bbb", "ccc"])

    def test_revise_lmp_input_model_0(self):
        lines = ["foo\n", "pair_style deepmd  aaa ccc fff\n", "bar\n", "\n"]
        ref_lines = copy.deepcopy(lines)
        lines = revise_lmp_input_model(lines, ["model0", "model1"], 10, "0.1")
        for ii in [0, 2, 3]:
            self.assertEqual(lines[ii], ref_lines[ii])
        tmp = " ".join(lines[1].split())
        self.assertEqual(tmp, "pair_style deepmd model0 model1 10 model_devi.out")

    def test_revise_lmp_input_model_1(self):
        lines = ["foo\n", "pair_style deepmd aaa ccc fff\n", "bar\n", "\n"]
        ref_lines = copy.deepcopy(lines)
        lines = revise_lmp_input_model(lines, ["model0", "model1"], 10, "1")
        for ii in [0, 2, 3]:
            self.assertEqual(lines[ii], ref_lines[ii])
        tmp = " ".join(lines[1].split())
        self.assertEqual(
            tmp, "pair_style deepmd model0 model1 out_freq 10 out_file model_devi.out"
        )

    def test_revise_lmp_input_dump(self):
        lines = ["foo\n", "dump dpgen_dump ccc fff\n", "bar\n", "\n"]
        ref_lines = copy.deepcopy(lines)
        lines = revise_lmp_input_dump(lines, 10)
        for ii in [0, 2, 3]:
            self.assertEqual(lines[ii], ref_lines[ii])
        tmp = " ".join(lines[1].split())
        self.assertEqual(
            tmp, "dump dpgen_dump all custom 10 traj/*.lammpstrj id type x y z"
        )

    def test_revise_lmp_input_plm(self):
        lines = ["foo\n", "fix dpgen_plm ccc fff\n", "bar\n", "\n"]
        ref_lines = copy.deepcopy(lines)
        lines = revise_lmp_input_plm(lines, "input.plumed")
        for ii in [0, 2, 3]:
            self.assertEqual(lines[ii], ref_lines[ii])
        tmp = " ".join(lines[1].split())
        self.assertEqual(
            tmp,
            "fix dpgen_plm all plumed plumedfile input.plumed outfile output.plumed",
        )

    def test_revise_by_key(self):
        lines = ["foo\n", "aaa\n", "bar\n", "bbb\n", "\n"]
        ref_lines = copy.deepcopy(lines)
        lines = revise_by_keys(lines, ["aaa", "bbb"], ["ccc", "ddd"])
        for ii in [0, 2, 4]:
            self.assertEqual(lines[ii], ref_lines[ii])
        tmp = " ".join(lines[1].split())
        self.assertEqual(tmp, "ccc")
        tmp = " ".join(lines[3].split())
        self.assertEqual(tmp, "ddd")


class TestMakeMDAMBER(unittest.TestCase):
    def tearDown(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")

    def test_make_model_devi(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_amber_file) as fp:
            jdata = json.load(fp)
        with open(machine_file) as fp:
            mdata = json.load(fp)
        # TODO: these should be normalized in the main program
        jdata["sys_configs_prefix"] = os.path.abspath(jdata["sys_configs_prefix"])
        jdata["disang_prefix"] = os.path.abspath(jdata["disang_prefix"])
        jdata["mdin_prefix"] = os.path.abspath(jdata["mdin_prefix"])
        jdata["parm7_prefix"] = os.path.abspath(jdata["mdin_prefix"])
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        self._check_input(0)

    def test_restart_from_iter(self):
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        if os.path.isdir("iter.000001"):
            shutil.rmtree("iter.000001")
        with open(param_amber_file) as fp:
            jdata = json.load(fp)
        with open(machine_file) as fp:
            mdata = json.load(fp)
        jdata["model_devi_jobs"].append(
            {
                "sys_idx": [0],
                "restart_from_iter": 0,
            }
        )
        jdata["sys_configs_prefix"] = os.path.abspath(jdata["sys_configs_prefix"])
        jdata["disang_prefix"] = os.path.abspath(jdata["disang_prefix"])
        jdata["mdin_prefix"] = os.path.abspath(jdata["mdin_prefix"])
        jdata["parm7_prefix"] = os.path.abspath(jdata["mdin_prefix"])
        _make_fake_models(0, jdata["numb_models"])
        make_model_devi(0, jdata, mdata)
        _check_pb(self, 0)
        self._check_input(0)
        restart_text = "This is the fake restart file to test `restart_from_iter`"
        with open(
            os.path.join(
                "iter.%06d" % 0, "01.model_devi", "task.000.000000", "rc.rst7"
            ),
            "w",
        ) as fw:
            fw.write(restart_text)
        _make_fake_models(1, jdata["numb_models"])
        make_model_devi(1, jdata, mdata)
        _check_pb(self, 1)
        self._check_input(1)
        with open(
            os.path.join(
                "iter.%06d" % 1, "01.model_devi", "task.000.000000", "init.rst7"
            )
        ) as f:
            assert f.read() == restart_text

    def _check_input(self, iter_idx: int):
        md_dir = os.path.join("iter.%06d" % iter_idx, "01.model_devi")
        assert os.path.isfile(os.path.join(md_dir, "init0.mdin"))
        assert os.path.isfile(os.path.join(md_dir, "qmmm0.parm7"))
        tasks = glob.glob(os.path.join(md_dir, "task.*"))
        for tt in tasks:
            assert os.path.isfile(os.path.join(tt, "init.rst7"))
            assert os.path.isfile(os.path.join(tt, "TEMPLATE.disang"))


if __name__ == "__main__":
    unittest.main()

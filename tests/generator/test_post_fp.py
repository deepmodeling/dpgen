import json
import os
import shutil
import sys
import unittest

import dpdata
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "generator"
from .comp_sys import (
    CompLabeledSys,
)
from .context import (
    param_abacus_post_file,
    param_amber_file,
    param_cp2k_file,
    param_custom_fp_file,
    param_file,
    param_gaussian_file,
    param_pwmat_file,
    param_pwscf_file,
    param_siesta_file,
    post_fp,
    post_fp_vasp,
    setup_ele_temp,
    setUpModule,  # noqa: F401
)


class TestPostFPVasp(unittest.TestCase):
    def setUp(self):
        assert os.path.isdir(
            "out_data_post_fp_vasp"
        ), "out data for post fp vasp should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_vasp", "iter.000000")
        self.ref_coord = [[[0, 0, 0], [2.3, 2.3, 2.3]], [[0, 0, 0], [2.2, 2.3, 2.4]]]
        self.ref_cell = [4.6 * np.eye(3), 4.6 * np.eye(3)]
        # type_map = ["Mg", "Al"], Al OUTCAR provided
        self.ref_at = [1, 1]
        self.ref_e = [-1.90811235, -1.89718546]
        self.ref_f = [
            [[0.0, 0.0, 0.0], [-0.0, -0.0, -0.0]],
            [[-0.110216, 0.0, 0.110216], [0.110216, -0.0, -0.110216]],
        ]
        self.ref_v = [
            [[1.50816698, 0.0, -0.0], [0.0, 1.50816698, 0.0], [-0.0, 0.0, 1.50816795]],
            [
                [1.45208913, 0.0, 0.03036584],
                [0.0, 1.67640928, 0.0],
                [0.03036584, 0.0, 1.45208913],
            ],
        ]
        self.ref_coord = np.array(self.ref_coord)
        self.ref_cell = np.array(self.ref_cell)
        self.ref_at = np.array(self.ref_at, dtype=int)
        self.ref_e = np.array(self.ref_e)
        self.ref_f = np.array(self.ref_f)
        self.ref_v = np.array(self.ref_v)
        # backup dpdata system data type
        self._system_dtypes = dpdata.System.DTYPES
        self._labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

    def tearDown(self):
        shutil.rmtree("iter.000000")
        # recover
        dpdata.System.DTYPES = self._system_dtypes
        dpdata.LabeledSystem.DTYPES = self._labeled_system_dtypes

    def test_post_fp_vasp_0(self):
        with open(param_file) as fp:
            jdata = json.load(fp)
        jdata["use_ele_temp"] = 2
        setup_ele_temp(True)
        post_fp_vasp(0, jdata, rfailed=0.3)

        sys = dpdata.LabeledSystem("iter.000000/02.fp/data.000/", fmt="deepmd/raw")
        self.assertEqual(sys.get_nframes(), 2)

        if sys.data["coords"][0][1][0] < sys.data["coords"][1][1][0]:
            idx = [1, 0]
        else:
            idx = [0, 1]
        ref_coord = self.ref_coord[idx]
        ref_cell = self.ref_cell[idx]
        ref_e = self.ref_e[idx]
        ref_f = self.ref_f[idx]
        ref_v = self.ref_v[idx]
        ref_at = self.ref_at

        for ff in range(2):
            self.assertAlmostEqual(ref_e[ff], sys.data["energies"][ff])
        for ii in range(2):
            self.assertEqual(ref_at[ff], sys.data["atom_types"][ff])
        for ff in range(2):
            for ii in range(2):
                for dd in range(3):
                    self.assertAlmostEqual(
                        ref_coord[ff][ii][dd], sys.data["coords"][ff][ii][dd]
                    )
                    self.assertAlmostEqual(
                        ref_f[ff][ii][dd], sys.data["forces"][ff][ii][dd]
                    )
        for ff in range(2):
            for ii in range(3):
                for jj in range(3):
                    self.assertAlmostEqual(
                        ref_v[ff][ii][jj], sys.data["virials"][ff][ii][jj], places=5
                    )
                    self.assertAlmostEqual(
                        ref_cell[ff][ii][jj], sys.data["cells"][ff][ii][jj]
                    )

        self.assertTrue(os.path.isfile("iter.000000/02.fp/data.000/set.000/aparam.npy"))
        aparam = np.load("iter.000000/02.fp/data.000/set.000/aparam.npy")
        natoms = sys.get_natoms()
        self.assertEqual(natoms, 2)
        self.assertEqual(list(list(aparam)[0]), [0, 0])
        self.assertEqual(list(list(aparam)[1]), [1, 1])

    def test_post_fp_vasp_1(self):
        with open(param_file) as fp:
            jdata = json.load(fp)
        jdata["use_ele_temp"] = 1
        setup_ele_temp(False)
        post_fp_vasp(0, jdata, rfailed=0.3)

        sys = dpdata.LabeledSystem("iter.000000/02.fp/data.001/", fmt="deepmd/raw")
        self.assertEqual(sys.get_nframes(), 1)

        # if sys.data['coords'][0][1][0] < sys.data['coords'][1][1][0]:
        #     idx = [0]
        # else :
        idx = [1]
        ref_coord = self.ref_coord[idx]
        ref_cell = self.ref_cell[idx]
        ref_e = self.ref_e[idx]
        ref_f = self.ref_f[idx]
        ref_v = self.ref_v[idx]
        ref_at = self.ref_at

        for ff in range(1):
            self.assertAlmostEqual(ref_e[ff], sys.data["energies"][ff])
        for ii in range(2):
            self.assertEqual(ref_at[ff], sys.data["atom_types"][ff])
        for ff in range(1):
            for ii in range(2):
                for dd in range(3):
                    self.assertAlmostEqual(
                        ref_coord[ff][ii][dd], sys.data["coords"][ff][ii][dd]
                    )
                    self.assertAlmostEqual(
                        ref_f[ff][ii][dd], sys.data["forces"][ff][ii][dd]
                    )
        for ff in range(1):
            for ii in range(3):
                for jj in range(3):
                    self.assertAlmostEqual(
                        ref_v[ff][ii][jj], sys.data["virials"][ff][ii][jj], places=5
                    )
                    self.assertAlmostEqual(
                        ref_cell[ff][ii][jj], sys.data["cells"][ff][ii][jj]
                    )

        fparam = np.load("iter.000000/02.fp/data.001/set.000/fparam.npy")
        self.assertEqual(fparam.shape[0], 1)
        self.assertEqual(list(fparam), [100000])

    def test_post_fp_vasp_2(self):
        with open(param_file) as fp:
            jdata = json.load(fp)
        jdata["use_ele_temp"] = 1
        with self.assertRaises(RuntimeError):
            post_fp_vasp(0, jdata)


class TestPostFPPWSCF(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir(
            "out_data_post_fp_pwscf"
        ), "out data for post fp pwscf should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_pwscf", "iter.000000")
        with open(param_pwscf_file) as fp:
            jdata = json.load(fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem("iter.000000/orig", fmt="deepmd/raw")
        self.system_2 = dpdata.LabeledSystem(
            "iter.000000/02.fp/data.000", fmt="deepmd/raw"
        )


class TestPostFPABACUS(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir(
            "out_data_post_fp_abacus"
        ), "out data for post fp pwscf should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_abacus", "iter.000000")
        with open(param_abacus_post_file) as fp:
            jdata = json.load(fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem("iter.000000/orig", fmt="deepmd/raw")
        self.system_2 = dpdata.LabeledSystem(
            "iter.000000/02.fp/data.000", fmt="deepmd/raw"
        )

    def test_nframs_with_failed_job(self):
        self.assertEqual(self.system_2.get_nframes(), 2)


class TestPostFPSIESTA(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir(
            "out_data_post_fp_siesta"
        ), "out data for post fp siesta should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_siesta", "iter.000000")
        with open(param_siesta_file) as fp:
            jdata = json.load(fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem("iter.000000/orig", fmt="deepmd/raw")
        self.system_2 = dpdata.LabeledSystem(
            "iter.000000/02.fp/data.000", fmt="deepmd/raw"
        )


class TestPostGaussian(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir(
            "out_data_post_fp_gaussian"
        ), "out data for post fp gaussian should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_gaussian", "iter.000000")
        with open(param_gaussian_file) as fp:
            jdata = json.load(fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem("iter.000000/orig", fmt="deepmd/raw")
        self.system_2 = dpdata.LabeledSystem(
            "iter.000000/02.fp/data.000", fmt="deepmd/raw"
        )


class TestPostCP2K(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir(
            "out_data_post_fp_cp2k"
        ), "out data for post fp cp2k should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_cp2k", "iter.000000")
        with open(param_cp2k_file) as fp:
            jdata = json.load(fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem("iter.000000/orig", fmt="deepmd/raw")
        self.system_2 = dpdata.LabeledSystem(
            "iter.000000/02.fp/data.000", fmt="deepmd/raw"
        )


class TestPostFPPWmat(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir(
            "out_data_post_fp_pwmat"
        ), "out data for post fp pwmat should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        shutil.copytree("out_data_post_fp_pwmat", "iter.000000")
        with open(param_pwmat_file) as fp:
            jdata = json.load(fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem("iter.000000/orig", fmt="deepmd/raw")
        self.system_2 = dpdata.LabeledSystem(
            "iter.000000/02.fp/data.000", fmt="deepmd/raw"
        )


class TestPostAmberDiff(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5

        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        ms = dpdata.MultiSystems(
            dpdata.LabeledSystem(os.path.join("data", "deepmd"), fmt="deepmd/raw")
        )
        ms.to_deepmd_npy(
            os.path.join("iter.000000", "02.fp", "task.000.000000", "dataset")
        )
        self.system_1 = list(ms.systems.values())[0]
        with open(param_amber_file) as fp:
            jdata = json.load(fp)
        jdata["type_map"] = self.system_1.get_atom_names()
        post_fp(0, jdata)
        self.system_2 = list(
            dpdata.MultiSystems(type_map=jdata["type_map"])
            .from_deepmd_raw("iter.000000/02.fp/data.000")
            .systems.values()
        )[0]


class TestPostFPCustom(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir(
            "out_data_post_fp_pwmat"
        ), "out data for post fp pwmat should exist"
        if os.path.isdir("iter.000000"):
            shutil.rmtree("iter.000000")
        with open(param_custom_fp_file) as fp:
            jdata = json.load(fp)
        fp_params = jdata["fp_params"]
        output_fn = fp_params["output_fn"]
        output_fmt = fp_params["output_fmt"]
        type_map = jdata["type_map"] + ["Type_0"]
        ss = dpdata.LabeledSystem(
            os.path.join("data", "deepmd"), fmt="deepmd/raw", type_map=type_map
        )
        output_filename = os.path.join(
            "iter.000000", "02.fp", "task.000.000000", output_fn
        )
        os.makedirs(os.path.dirname(output_filename), exist_ok=True)
        ss.to(output_fmt, output_filename)
        post_fp(0, jdata)
        self.system_1 = ss
        self.system_2 = list(
            dpdata.MultiSystems(type_map=type_map)
            .from_deepmd_raw("iter.000000/02.fp/data.000")
            .systems.values()
        )[0]


if __name__ == "__main__":
    unittest.main()

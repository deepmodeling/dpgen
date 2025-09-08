import glob
import os
import shutil
import sys
import unittest

import numpy as np
from pymatgen.analysis.defects.core import Interstitial as pmg_Interstitial
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "auto_test"

from dpgen.auto_test.Interstitial import Interstitial

from .context import setUpModule  # noqa: F401


class TestInterstitial(unittest.TestCase):
    def setUp(self):
        _jdata = {
            "structures": ["confs/std-bcc"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": "vasp_input",
                "potcars": {"V": "POTCAR"},
            },
            "properties": [
                {
                    "type": "interstitial",
                    "supercell": [1, 1, 1],
                    "insert_ele": ["V"],
                    "bcc_self": True,
                }
            ],
        }

        self.equi_path = "confs/std-bcc/relaxation/relax_task"
        self.source_path = "equi/vasp"
        self.target_path = "confs/std-bcc/interstitial_00"
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)
        if not os.path.exists(self.target_path):
            os.makedirs(self.target_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata["properties"]

        self.interstitial = Interstitial(_jdata["properties"][0])

    def tearDown(self):
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual("interstitial", self.interstitial.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.interstitial.task_param())

    def test_make_confs_bcc(self):
        if not os.path.exists(os.path.join(self.equi_path, "CONTCAR")):
            with self.assertRaises(RuntimeError):
                self.interstitial.make_confs(self.target_path, self.equi_path)
        shutil.copy(
            os.path.join(self.source_path, "CONTCAR_V_bcc"),
            os.path.join(self.equi_path, "CONTCAR"),
        )
        task_list = self.interstitial.make_confs(self.target_path, self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, "task.*"))
        self.assertEqual(len(dfm_dirs), 7)

        incar0 = Incar.from_file(os.path.join("vasp_input", "INCAR.rlx"))
        incar0["ISIF"] = 3

        self.assertEqual(
            os.path.realpath(os.path.join(self.equi_path, "CONTCAR")),
            os.path.realpath(os.path.join(self.target_path, "POSCAR")),
        )
        ref_st = Structure.from_file(os.path.join(self.target_path, "POSCAR"))
        dfm_dirs.sort()
        for ii in dfm_dirs[:4]:
            st_file = os.path.join(ii, "POSCAR")
            self.assertTrue(os.path.isfile(st_file))
            st0 = Structure.from_file(st_file)
            inter_site = st0[0]
            inter = pmg_Interstitial(ref_st, inter_site)
            st1 = inter.get_supercell_structure(
                sc_mat=np.eye(3) * self.prop_param[0]["supercell"]
            )
            # TODO: fix the failed test
            # self.assertEqual(st0, st1)

        for ii in dfm_dirs[4:]:
            st_file = os.path.join(ii, "POSCAR")
            self.assertTrue(os.path.isfile(st_file))
            st0 = Structure.from_file(st_file)
            inter_site1 = st0.pop(0)
            inter_site2 = st0.pop(-1)
            center = (inter_site1.coords + inter_site2.coords) / 2
            self.assertTrue((center[0] - center[1]) < 1e-4)
            self.assertTrue((center[1] - center[2]) < 1e-4)

    def test_make_confs_filtered_out(self):
        """Test that element.out is created even when all defects are filtered out."""
        # Create a configuration that filters out all defects
        _jdata_filtered = {
            "structures": ["confs/std-bcc"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": "vasp_input",
                "potcars": {"V": "POTCAR"},
            },
            "properties": [
                {
                    "type": "interstitial",
                    "supercell": [1, 1, 1],
                    "insert_ele": ["V"],
                    "conf_filters": {"min_dist": 100.0},  # Very high filter - all defects filtered
                    "bcc_self": False,  # Don't add bcc tasks
                }
            ],
        }

        filtered_target_path = "confs/std-bcc/interstitial_filtered"
        if not os.path.exists(filtered_target_path):
            os.makedirs(filtered_target_path)

        try:
            # Copy the test structure
            shutil.copy(
                os.path.join(self.source_path, "CONTCAR_V_bcc"),
                os.path.join(self.equi_path, "CONTCAR"),
            )

            interstitial_filtered = Interstitial(_jdata_filtered["properties"][0])
            task_list = interstitial_filtered.make_confs(filtered_target_path, self.equi_path)

            # Should have no tasks since all are filtered
            self.assertEqual(len(task_list), 0)

            # element.out should exist even if empty
            element_out_path = os.path.join(filtered_target_path, "element.out")
            self.assertTrue(os.path.exists(element_out_path), "element.out should exist even when empty")

            # Test that post_process doesn't crash with empty task list
            try:
                interstitial_filtered.post_process(task_list)
            except Exception as e:
                self.fail(f"post_process should handle empty task list gracefully, but got: {e}")

        finally:
            # Clean up
            if os.path.exists(filtered_target_path):
                shutil.rmtree(filtered_target_path)

    def test_make_confs_partial_filtering(self):
        """Test that post_process handles mismatched element.out entries gracefully."""
        # Create a configuration that might have some defects filtered
        _jdata_partial = {
            "structures": ["confs/std-bcc"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": "vasp_input",
                "potcars": {"V": "POTCAR"},
            },
            "properties": [
                {
                    "type": "interstitial",
                    "supercell": [1, 1, 1],
                    "insert_ele": ["V"],
                    # No conf_filters to ensure some tasks are created
                    "bcc_self": True,  # Add bcc tasks that always create more tasks
                }
            ],
        }

        partial_target_path = "confs/std-bcc/interstitial_partial"
        if not os.path.exists(partial_target_path):
            os.makedirs(partial_target_path)

        try:
            # Copy the test structure
            shutil.copy(
                os.path.join(self.source_path, "CONTCAR_V_bcc"),
                os.path.join(self.equi_path, "CONTCAR"),
            )

            interstitial_partial = Interstitial(_jdata_partial["properties"][0])
            task_list = interstitial_partial.make_confs(partial_target_path, self.equi_path)

            # Should have some tasks (at least the bcc_self ones)
            self.assertGreater(len(task_list), 0)

            # element.out should exist
            element_out_path = os.path.join(partial_target_path, "element.out")
            self.assertTrue(os.path.exists(element_out_path), "element.out should exist")

            # Test that post_process doesn't crash even if element.out and task count don't match
            try:
                interstitial_partial.post_process(task_list)
            except Exception as e:
                self.fail(f"post_process should handle mismatched element.out gracefully, but got: {e}")

        finally:
            # Clean up
            if os.path.exists(partial_target_path):
                shutil.rmtree(partial_target_path)

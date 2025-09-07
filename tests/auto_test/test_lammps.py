import os
import shutil
import sys
import tempfile
import unittest

from monty.serialization import loadfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "auto_test"

from dpgen.auto_test.common_equi import make_equi, run_equi
from dpgen.auto_test.Lammps import Lammps
from dpgen.auto_test.lib.lammps import inter_deepmd

from .context import setUpModule  # noqa: F401


class TestLammps(unittest.TestCase):
    def setUp(self):
        self.jdata = {
            "structures": ["confs/std-fcc"],
            "interaction": {
                "type": "deepmd",
                "model": "lammps_input/frozen_model.pb",
                "deepmd_version": "1.1.0",
                "type_map": {"Al": 0},
            },
            "relaxation": {
                "cal_type": "relaxation",
                "cal_setting": {
                    "relax_pos": True,
                    "relax_shape": True,
                    "relax_vol": True,
                },
            },
        }

        self.equi_path = "confs/std-fcc/relaxation/relax_task"
        self.source_path = "equi/lammps"

        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        if not os.path.isfile(os.path.join(self.equi_path, "POSCAR")):
            shutil.copy(
                os.path.join(self.source_path, "Al-fcc.vasp"),
                os.path.join("confs/std-fcc", "POSCAR"),
            )

        self.confs = self.jdata["structures"]
        self.inter_param = self.jdata["interaction"]
        self.relax_param = self.jdata["relaxation"]
        self.Lammps = Lammps(
            self.inter_param, os.path.join(self.source_path, "Al-fcc.vasp")
        )

    def tearDown(self):
        if os.path.exists("confs/std-fcc/relaxation"):
            shutil.rmtree("confs/std-fcc/relaxation")

    def test_set_inter_type_func(self):
        self.Lammps.set_inter_type_func()
        self.assertEqual(inter_deepmd, self.Lammps.inter_func)

    def test_set_model_param(self):
        self.Lammps.set_model_param()
        model_param = {
            "model_name": ["frozen_model.pb"],
            "param_type": {"Al": 0},
            "deepmd_version": "1.1.0",
        }
        self.assertEqual(model_param, self.Lammps.model_param)

    def test_make_potential_files(self):
        cwd = os.getcwd()
        abs_equi_path = os.path.abspath(self.equi_path)
        self.Lammps.make_potential_files(abs_equi_path)
        self.assertTrue(os.path.islink(os.path.join(self.equi_path, "frozen_model.pb")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "inter.json")))
        ret = loadfn(os.path.join(self.equi_path, "inter.json"))
        self.assertTrue(self.inter_param, ret)
        os.chdir(cwd)

    def test_make_input_file(self):
        cwd = os.getcwd()
        abs_equi_path = os.path.abspath("confs/std-fcc/relaxation/relax_task")
        shutil.copy(
            os.path.join("confs/std-fcc", "POSCAR"),
            os.path.join(self.equi_path, "POSCAR"),
        )
        self.Lammps.make_input_file(abs_equi_path, "relaxation", self.relax_param)
        self.assertTrue(os.path.isfile(os.path.join(abs_equi_path, "conf.lmp")))
        self.assertTrue(os.path.islink(os.path.join(abs_equi_path, "in.lammps")))
        self.assertTrue(os.path.isfile(os.path.join(abs_equi_path, "task.json")))

    def test_forward_common_files(self):
        fc_files = ["in.lammps", "frozen_model.pb"]
        self.assertEqual(self.Lammps.forward_common_files(), fc_files)

    def test_backward_files(self):
        backward_files = ["log.lammps", "outlog", "dump.relax"]
        self.assertEqual(self.Lammps.backward_files(), backward_files)

    def test_forward_common_files_with_custom_in_lammps(self):
        """Test forward_common_files with custom in_lammps path (fixes #1757)."""
        # Test custom in_lammps path with different filename (not "in.lammps")
        custom_inter_param = {
            "type": "deepmd",
            "model": "lammps_input/frozen_model.pb",
            "in_lammps": "lammps_input/custom_input.lmp",  # Custom path with different filename
            "deepmd_version": "1.1.0",
            "type_map": {"Al": 0},
        }

        lammps_custom = Lammps(custom_inter_param, self.source_path + "/Al-fcc.vasp")
        fc_files_custom = lammps_custom.forward_common_files()
        expected_custom = ["lammps_input/custom_input.lmp", "frozen_model.pb"]
        self.assertEqual(fc_files_custom, expected_custom)

        # Test that forward_files also uses custom path
        forward_files_custom = lammps_custom.forward_files()
        expected_forward = [
            "conf.lmp",
            "in.lammps",
            "frozen_model.pb",
        ]  # Uses "in.lammps", not custom path
        self.assertEqual(forward_files_custom, expected_forward)

        # Test EOS property type (should not include in.lammps)
        fc_files_eos = lammps_custom.forward_common_files(property_type="eos")
        expected_eos = ["frozen_model.pb"]
        self.assertEqual(fc_files_eos, expected_eos)

        # Test with another completely different filename to ensure robustness
        alternate_inter_param = {
            "type": "deepmd",
            "model": "frozen_model.pb",
            "in_lammps": "input_files/my_lammps_script.in",  # Different extension and path
            "deepmd_version": "1.1.0",
            "type_map": {"Al": 0},
        }

        lammps_alternate = Lammps(
            alternate_inter_param, self.source_path + "/Al-fcc.vasp"
        )
        fc_files_alternate = lammps_alternate.forward_common_files()
        expected_alternate = ["input_files/my_lammps_script.in", "frozen_model.pb"]
        self.assertEqual(fc_files_alternate, expected_alternate)

        forward_files_alternate = lammps_alternate.forward_files()
        expected_forward_alternate = [
            "conf.lmp",
            "in.lammps",
            "frozen_model.pb",
        ]  # Uses "in.lammps", not custom path
        self.assertEqual(forward_files_alternate, expected_forward_alternate)

    def test_run_equi_with_custom_in_lammps(self):
        """Test run_equi locally with custom in_lammps path to verify fix works end-to-end."""
        # Create temporary directories for the test
        test_work_dir = "test_custom_lammps_run"
        if os.path.exists(test_work_dir):
            shutil.rmtree(test_work_dir)
        os.makedirs(test_work_dir)

        try:
            with tempfile.TemporaryDirectory() as remote_root:
                # Set up configuration structure - use existing fcc-Al conf
                conf_dir = os.path.join(test_work_dir, "confs", "fcc-Al")
                os.makedirs(conf_dir, exist_ok=True)

                # Copy structure file from existing conf (use POSCAR file for deepmd interaction)
                test_dir = os.path.dirname(os.path.abspath(__file__))
                poscar_path = os.path.join(test_dir, "equi", "lammps", "Al-fcc.vasp")
                shutil.copy(poscar_path, os.path.join(conf_dir, "POSCAR"))

                # Create the custom in_lammps file
                custom_lammps_dir = os.path.join(test_work_dir, "input_files")
                os.makedirs(custom_lammps_dir, exist_ok=True)
                custom_lammps_file = os.path.join(custom_lammps_dir, "my_custom.lmp")
                custom_lammps_source = os.path.join(
                    test_dir, "lammps_input", "custom_lammps_input.lmp"
                )
                shutil.copy(custom_lammps_source, custom_lammps_file)

                # Create a dummy frozen model
                model_file = os.path.join(test_work_dir, "frozen_model.pb")
                with open(model_file, "w") as f:
                    f.write("dummy model file for testing")

                # Change to test directory
                original_cwd = os.getcwd()
                os.chdir(test_work_dir)

                try:
                    # Test configuration with custom in_lammps path
                    jdata = {
                        "structures": ["confs/fcc-Al"],
                        "interaction": {
                            "type": "deepmd",
                            "model": "frozen_model.pb",
                            "in_lammps": "input_files/my_custom.lmp",  # Custom path with different filename
                            "deepmd_version": "1.1.0",
                            "type_map": {"Al": 0},
                        },
                        "relaxation": {
                            "cal_type": "relaxation",
                            "cal_setting": {
                                "relax_pos": True,
                                "relax_shape": True,
                                "relax_vol": True,
                            },
                        },
                    }

                    # Machine configuration for LocalContext
                    mdata = {
                        "model_devi_command": "touch log.lammps dump.relax; echo lmp",
                        "model_devi_machine": {
                            "context_type": "LocalContext",
                            "batch_type": "shell",
                            "local_root": "./",
                            "remote_root": remote_root,
                        },
                        "model_devi_resources": {
                            "group_size": 1,
                        },
                        "model_devi_group_size": 1,
                    }

                    # First create the equi structure
                    make_equi(
                        jdata["structures"], jdata["interaction"], jdata["relaxation"]
                    )

                    # Verify that the task was created correctly
                    task_dir = "confs/fcc-Al/relaxation/relax_task"
                    self.assertTrue(os.path.exists(task_dir))
                    self.assertTrue(
                        os.path.exists(os.path.join(task_dir, "inter.json"))
                    )

                    # Test the key functionality: verify that forward_common_files
                    # returns the custom path, which proves the fix works
                    virtual_calculator = Lammps(
                        jdata["interaction"], "confs/std-fcc/POSCAR"
                    )
                    forward_common = virtual_calculator.forward_common_files()

                    # This is the core test - the custom path should be included
                    self.assertIn("input_files/my_custom.lmp", forward_common)
                    self.assertIn("frozen_model.pb", forward_common)

                    # Verify that forward_files still uses "in.lammps" (for task directories)
                    forward_files = virtual_calculator.forward_files()
                    self.assertIn("in.lammps", forward_files)
                    self.assertNotIn("input_files/my_custom.lmp", forward_files)

                    # Verify that the custom file actually exists at the specified path
                    self.assertTrue(os.path.exists("input_files/my_custom.lmp"))

                    # Now actually run run_equi to test end-to-end functionality
                    # This will fail with FileNotFoundError if the fix doesn't work
                    run_equi(jdata["structures"], jdata["interaction"], mdata)

                    # If we reach here without exception, the fix works correctly
                    # The dispatcher was able to find the custom file because
                    # forward_common_files() returned the correct path

                finally:
                    os.chdir(original_cwd)

        finally:
            # Clean up test directory
            if os.path.exists(test_work_dir):
                shutil.rmtree(test_work_dir)

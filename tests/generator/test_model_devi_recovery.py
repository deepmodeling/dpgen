import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from dpgen.generator.run import (
    _normalize_model_devi_recovery,
    _read_model_devi_file,
    _recovery_task_report,
    _recovery_task_states,
    _select_by_model_devi_adaptive_trust_low,
    _select_by_model_devi_standard,
    post_model_devi,
)

DUMP = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 5
0 5
0 5
ITEM: ATOMS id type x y z fx fy fz
1 1 0 0 0 0 0 0
2 1 1 0 0 0 0 0
"""

MISMATCHED_TYPE_DUMP = DUMP.replace("2 1 1 0 0 0 0 0", "2 2 1 0 0 0 0 0")


class TestModelDeviRecovery(unittest.TestCase):
    def test_policy_defaults_and_validation(self):
        policy = _normalize_model_devi_recovery(
            {"model_devi_recovery": {"enabled": True}}
        )
        self.assertTrue(policy["enabled"])
        self.assertEqual(policy["max_failed_tasks"], 0)
        with self.assertRaises(ValueError):
            _normalize_model_devi_recovery(
                {"model_devi_recovery": {"max_failed_ratio": 2}}
            )
        with self.assertRaisesRegex(ValueError, "lammps"):
            _normalize_model_devi_recovery(
                {
                    "model_devi_engine": "gromacs",
                    "model_devi_recovery": {"enabled": True},
                }
            )

    def test_task_report_salvages_matching_prefix(self):
        with tempfile.TemporaryDirectory() as root:
            task = Path(root) / "task.000.000000"
            (task / "traj").mkdir(parents=True)
            (task / "traj" / "0.lammpstrj").write_text(DUMP)
            np.savetxt(
                task / "model_devi.out",
                np.asarray([[0, 0.1, 0, 0, 0.2, 0, 0], [10, 0.1, 0, 0, 0.2, 0, 0]]),
            )
            report = _recovery_task_report(task, False, True)
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["valid_steps"], [0])
            self.assertEqual(report["excluded_steps"], [10])

    def test_task_report_salvages_only_contiguous_prefix(self):
        with tempfile.TemporaryDirectory() as root:
            task = Path(root) / "task.000.000000"
            (task / "traj").mkdir(parents=True)
            for step in (0, 20):
                (task / "traj" / f"{step}.lammpstrj").write_text(
                    DUMP.replace("\n0\n", f"\n{step}\n", 1)
                )
            np.savetxt(
                task / "model_devi.out",
                np.asarray([[step, 0.1, 0, 0, 0.2, 0, 0] for step in (0, 10, 20)]),
            )
            report = _recovery_task_report(task, False, True)
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["valid_steps"], [0])
            self.assertEqual(report["excluded_steps"], [10, 20])

    def test_selection_reads_report_allow_list(self):
        with tempfile.TemporaryDirectory() as root:
            iteration = Path(root) / "iter.000000"
            task = iteration / "01.model_devi" / "task.000.000000"
            task.mkdir(parents=True)
            np.savetxt(
                task / "model_devi.out",
                np.asarray([[0, 0.1, 0, 0, 0.2, 0, 0], [10, 0.1, 0, 0, 0.2, 0, 0]]),
            )
            (iteration / "exploration_report.json").write_text(
                json.dumps({"tasks": [{"task": task.name, "valid_steps": [0]}]})
            )
            result = _read_model_devi_file(str(task))
            self.assertEqual(result.shape, (1, 7))
            self.assertEqual(int(result[0, 0]), 0)

    def test_selection_skips_empty_allow_list_before_loading_output(self):
        with tempfile.TemporaryDirectory() as root:
            iteration = Path(root) / "iter.000000"
            task = iteration / "01.model_devi" / "task.000.000000"
            task.mkdir(parents=True)
            (iteration / "exploration_report.json").write_text(
                json.dumps({"tasks": [{"task": task.name, "valid_steps": []}]})
            )
            standard = _select_by_model_devi_standard(
                [str(task)], 0.0, 1.0, 0.0, 1.0, None, "lammps"
            )
            adaptive = _select_by_model_devi_adaptive_trust_low(
                [str(task)], 1.0, 1, 0.0, 1.0, 1, 0.0
            )
            self.assertEqual(standard[-1]["candidate"], 0)
            self.assertEqual(adaptive[3]["candidate"], 0)

    def test_post_and_selection_tolerate_missing_task_output(self):
        with tempfile.TemporaryDirectory() as root:
            previous = Path.cwd()
            os.chdir(root)
            self.addCleanup(os.chdir, previous)
            try:
                work_path = Path("iter.000000") / "01.model_devi"
                valid_task = work_path / "task.000.000000"
                missing_task = work_path / "task.001.000000"
                (valid_task / "traj").mkdir(parents=True)
                missing_task.mkdir(parents=True)
                (valid_task / "traj" / "0.lammpstrj").write_text(DUMP)
                np.savetxt(
                    valid_task / "model_devi.out",
                    np.asarray([[0, 0.1, 0, 0, 0.2, 0, 0]]),
                )
                post_model_devi(
                    0,
                    {
                        "model_devi_recovery": {
                            "enabled": True,
                            "max_failed_tasks": 1,
                            "max_failed_ratio": 0.5,
                        }
                    },
                    {},
                )
                report = json.loads(
                    Path("iter.000000/exploration_report.json").read_text()
                )
                self.assertEqual(report["failed_tasks"], 1)
                selected = _select_by_model_devi_standard(
                    [str(valid_task), str(missing_task)],
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                    None,
                    "lammps",
                )
                self.assertEqual(selected[-1]["candidate"], 1)
            finally:
                os.chdir(previous)

    def test_task_report_preserves_terminal_failure_and_restart_state(self):
        with tempfile.TemporaryDirectory() as root:
            previous = Path.cwd()
            try:
                os.chdir(root)
                iteration = Path("iter.000000")
                work_path = iteration / "01.model_devi"
                task = work_path / "task.000.000000"
                (task / "traj").mkdir(parents=True)
                (task / "traj" / "0.lammpstrj").write_text(DUMP)
                np.savetxt(
                    task / "model_devi.out",
                    np.asarray(
                        [
                            [0, 0.1, 0, 0, 0.2, 0, 0],
                            [10, 0.1, 0, 0, 0.2, 0, 0],
                        ]
                    ),
                )
                (work_path / "cur_job.json").write_text(
                    json.dumps({"nsteps": 1000, "trj_freq": 10})
                )
                (work_path / "recovery_task_states.json").write_text(
                    json.dumps({"tasks": {task.name: 7}})
                )
                policy = {
                    "model_devi_recovery": {
                        "enabled": True,
                        "max_failed_tasks": 1,
                        "max_failed_ratio": 1.0,
                        "salvage_prefix": True,
                    }
                }
                post_model_devi(0, policy, {})
                report = json.loads((iteration / "exploration_report.json").read_text())
                self.assertEqual(report["tasks"][0]["status"], "failed")
                self.assertEqual(report["tasks"][0]["valid_steps"], [0])
                self.assertIn("failure state", report["tasks"][0]["reason"])
                self.assertIn("1000", report["tasks"][0]["reason"])
                post_model_devi(0, policy, {})
                repeated = json.loads(
                    (iteration / "exploration_report.json").read_text()
                )
                self.assertEqual(repeated["tasks"][0]["status"], "failed")
                self.assertEqual(repeated["tasks"][0]["valid_steps"], [0])
                self.assertEqual(repeated["tasks"][0]["excluded_steps"], [10])
            finally:
                os.chdir(previous)

    def test_task_report_uses_last_scheduled_dump_step(self):
        with tempfile.TemporaryDirectory() as root:
            task = Path(root) / "task.000.000000"
            (task / "traj").mkdir(parents=True)
            steps = [0, 30, 60, 90]
            for step in steps:
                (task / "traj" / f"{step}.lammpstrj").write_text(
                    DUMP.replace("\n0\n", f"\n{step}\n", 1)
                )
            np.savetxt(
                task / "model_devi.out",
                np.asarray([[step, 0.1, 0, 0, 0.2, 0, 0] for step in steps]),
            )
            report = _recovery_task_report(
                task, False, False, task_state=5, expected_last_step=90
            )
            self.assertEqual(report["status"], "completed")
            self.assertEqual(report["valid_steps"], steps)

    def test_task_report_validates_atom_mapping_across_split_files(self):
        with tempfile.TemporaryDirectory() as root:
            task = Path(root) / "task.000.000000"
            (task / "traj").mkdir(parents=True)
            (task / "traj" / "0.lammpstrj").write_text(DUMP)
            (task / "traj" / "10.lammpstrj").write_text(
                MISMATCHED_TYPE_DUMP.replace("\n0\n", "\n10\n", 1)
            )
            np.savetxt(
                task / "model_devi.out",
                np.asarray([[0, 0.1, 0, 0, 0.2, 0, 0], [10, 0.1, 0, 0, 0.2, 0, 0]]),
            )
            report = _recovery_task_report(task, False, True)
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["valid_steps"], [0])
            self.assertEqual(report["excluded_steps"], [10])
            self.assertEqual(report["invalid_frames"], 1)

    def test_task_report_rejects_zero_volume_box(self):
        with tempfile.TemporaryDirectory() as root:
            task = Path(root) / "task.000.000000"
            (task / "traj").mkdir(parents=True)
            (task / "traj" / "0.lammpstrj").write_text(
                DUMP.replace("0 5\n0 5\n0 5", "0 0\n0 5\n0 5")
            )
            np.savetxt(
                task / "model_devi.out",
                np.asarray([[0, 0.1, 0, 0, 0.2, 0, 0]]),
            )
            report = _recovery_task_report(task, False, True)
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["valid_steps"], [])
            self.assertEqual(report["invalid_frames"], 1)

    def test_extracts_dispatcher_task_states(self):
        result = _recovery_task_states(
            {
                "belonging_jobs": [
                    {
                        "job-hash": {
                            "job_state": 7,
                            "task_states": [5, 7],
                            "job_task_list": [
                                {"task_work_path": "task.000.000000"},
                                {"task_work_path": "task.001.000000"},
                            ],
                        }
                    }
                ]
            }
        )
        self.assertEqual(result, {"task.000.000000": 5, "task.001.000000": 7})


if __name__ == "__main__":
    unittest.main()

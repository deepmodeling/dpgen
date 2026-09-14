import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

test_dir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.insert(0, os.path.join(test_dir, ".."))
__package__ = "tools"
from .context import (
    setUpModule,  # noqa: F401
    stat_sys,
)


class TestRunReport(unittest.TestCase):
    def test_stat_sys(self):
        folder = "run_report_test_output"
        sys, sys_count, sys_all = stat_sys(
            os.path.join(test_dir, folder), verbose=False, mute=True
        )
        with open(os.path.join(test_dir, folder, "param.json")) as fp:
            jdata = json.load(fp)
        self.assertEqual(sys, jdata["sys_configs"])
        self.assertEqual(sys_count, [jdata["fp_task_max"], jdata["fp_task_max"]])
        ref_all = [
            [
                ["npt", 50.0, 1.0, 4],
                ["npt", 50.0, 2.0, 1],
                ["npt", 100.0, 1.0, 2],
                ["npt", 100.0, 2.0, 1],
            ],
            [["npt", 50.0, 1.0, 2], ["npt", 50.0, 2.0, 4], ["npt", 100.0, 1.0, 2]],
        ]
        self.assertEqual(sys_all, ref_all)

    def test_stat_sys_restores_cwd_on_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)
            (target / "param.json").write_text(json.dumps({"sys_configs": ["sys"]}))
            cwd = Path.cwd()
            with (
                patch(
                    "dpgen.tools.stat_sys.glob.glob", side_effect=RuntimeError("boom")
                ),
                self.assertRaisesRegex(RuntimeError, "boom"),
            ):
                stat_sys(target, verbose=False, mute=True)
            self.assertEqual(Path.cwd(), cwd)

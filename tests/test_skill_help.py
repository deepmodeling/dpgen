"""Discoverable workflow guidance belongs to help, not command output."""

import contextlib
import io
import unittest
from unittest.mock import patch

from dpgen.main import main_parser
from dpgen.skills import workflow_help


class TestSkillHelp(unittest.TestCase):
    def test_contextual_help(self) -> None:
        for command, skill in {
            "run": "dpgen-run",
            "simplify": "dpgen-simplify",
        }.items():
            with self.subTest(command=command):
                output = io.StringIO()
                with (
                    contextlib.redirect_stdout(output),
                    self.assertRaises(SystemExit) as result,
                ):
                    main_parser().parse_args(
                        ([command] if command else []) + ["--help"]
                    )
                self.assertEqual(result.exception.code, 0)
                self.assertIn(f"dpgen skills read {skill}", output.getvalue())

    def test_missing_guide_has_no_pointer(self) -> None:
        self.assertEqual(workflow_help("not-a-bundled-skill"), "")
        with patch("dpgen.skills._root", side_effect=ValueError("missing resources")):
            self.assertEqual(workflow_help("dpgen-run"), "")

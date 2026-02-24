import subprocess as sp
import sys
import unittest


class TestCLI(unittest.TestCase):
    def test_cli(self):
        sp.check_output([sys.executable, "-m", "dpgen.main", "-h"])
        for subcommand in (
            "run",
            "simplify",
            "init_surf",
            "init_bulk",
            "init_reaction",
            "autotest",
        ):
            output = sp.check_output(
                [sys.executable, "-m", "dpgen.main", subcommand, "-h"]
            )
            if subcommand in ("run", "simplify", "init_reaction"):
                self.assertIn(b"--allow-ref", output)

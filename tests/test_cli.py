import subprocess as sp
import unittest


class TestCLI(unittest.TestCase):
    def test_cli(self):
        sp.check_output(["dpgen", "-h"])
        for subcommand in (
            "run",
            "simplify",
            "init_surf",
            "init_bulk",
            "init_reaction",
            "autotest",
        ):
            output = sp.check_output(["dpgen", subcommand, "-h"])
            if subcommand in ("run", "simplify", "init_reaction"):
                self.assertIn(b"--allow-ref", output)
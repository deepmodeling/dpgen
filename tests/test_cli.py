import subprocess as sp
import sys
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
            sp.check_output(["dpgen", subcommand, "-h"])

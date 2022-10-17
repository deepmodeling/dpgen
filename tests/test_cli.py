import sys
import subprocess as sp


class TestCLI(unittest.TestCase):
    def test_cli(self):
        sp.check_output([sys.executable, "-m", "dpgen", "-h"])
        for subcommand in ('run', 'simplify', 'init_surf', 'init_bulk', 'init_reaction', 'autotest'):
            sp.check_output([sys.executable, "-m", "dpgen", subcommand, "-h"])

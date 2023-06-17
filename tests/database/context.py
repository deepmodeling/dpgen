import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import dpgen  # noqa: F401
from dpgen.database.entry import Entry  # noqa: F401
from dpgen.database.run import parsing_vasp  # noqa: F401
from dpgen.database.vasp import DPPotcar, VaspInput  # noqa: F401


def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

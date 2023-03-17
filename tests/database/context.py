import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import dpgen
from dpgen.database.entry import Entry
from dpgen.database.run import parsing_vasp
from dpgen.database.vasp import DPPotcar, VaspInput


def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

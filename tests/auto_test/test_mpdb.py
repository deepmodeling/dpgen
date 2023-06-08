import os
import sys
import unittest

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "auto_test"
from dpgen.auto_test.mpdb import get_structure
from .context import setUpModule  # noqa: F401

try:
    os.environ["MAPI_KEY"]
    exist_key = True
except Exception:
    exist_key = False


def fit(struct0, struct1):
    m = StructureMatcher()
    if m.fit(struct0, struct1):
        return True
    return False


@unittest.skipIf(not exist_key, "skip mpdb")
class TestMpdb(unittest.TestCase):
    def setUp(self):
        if "MAPI_KEY" in os.environ:
            self.key = os.environ["MAPI_KEY"]
        else:
            self.key = None
        self.mpid = "mp-141"
        self.st_file = self.mpid + ".vasp"
        self.st0_file = os.path.join("confs/", self.mpid, self.mpid + ".cif")

    def tearDown(self):
        if os.path.exists(self.st_file):
            os.remove(self.st_file)

    def test_get_structure(self):
        st1 = get_structure(self.mpid)
        st0 = Structure.from_file(self.st0_file)
        self.assertTrue(fit(st0, st1))

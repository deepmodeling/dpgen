import glob
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "auto_test"

from .context import make_kspacing_kpoints
from .context import setUpModule  # noqa: F401


class TestVASPMakeKpoint(unittest.TestCase):
    def test_gf_all(self):
        kspacing = 0.16
        gamma = False
        all_test = glob.glob(os.path.join("data.vasp.kp.gf", "test.*"))
        for ii in all_test:
            ret = make_kspacing_kpoints(os.path.join(ii, "POSCAR"), kspacing, gamma)
            kp = [int(jj) for jj in (ret.split("\n")[3].split())]
            kp_ref = list(np.loadtxt(os.path.join(ii, "kp.ref"), dtype=int))
            self.assertTrue(kp == kp_ref)

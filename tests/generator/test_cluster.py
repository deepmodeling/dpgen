import importlib
import os
import sys
import unittest

import dpdata
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "generator"
from .comp_sys import CompSys
from .context import _crd2frag, take_cluster


@unittest.skipIf(importlib.util.find_spec("openbabel") is None, "requires openbabel")
class Test_take_cluster(unittest.TestCase, CompSys):
    def setUp(self):
        type_map = ["C", "H"]
        jdata = {"cluster_cutoff": 3.5}
        self.system_1 = take_cluster("cluster/14400.lammpstrj", type_map, 1125, jdata)
        self.system_2 = dpdata.System.load("cluster/cluster1.json")
        self.places = 0


@unittest.skipIf(importlib.util.find_spec("openbabel") is None, "requires openbabel")
class Test_take_cluster_minify(unittest.TestCase, CompSys):
    def setUp(self):
        type_map = ["C", "H"]
        jdata = {"cluster_cutoff": 3.5, "cluster_minify": True}
        self.system_1 = take_cluster("cluster/14400.lammpstrj", type_map, 1125, jdata)
        self.system_2 = dpdata.LabeledSystem(
            "cluster/input0_new.gaussianlog", fmt="gaussian/log"
        )
        self.system_2.data["cells"] = self.system_1["cells"]
        self.places = 0


class TestCrd2Frag(unittest.TestCase):
    def test_crd2frag_pbc(self):
        crds = np.array([[0.0, 0.0, 0.0], [19.0, 19.0, 19.0]])
        symbols = ["O", "O"]
        cell = np.diag([20.0, 20.0, 20.0])
        frag_numb, _ = _crd2frag(symbols, crds, pbc=True, cell=cell)
        self.assertEqual(frag_numb, 1)

    def test_crd2frag_nopbc(self):
        crds = np.array([[0.0, 0.0, 0.0], [19.0, 19.0, 19.0]])
        symbols = ["O", "O"]
        frag_numb, _ = _crd2frag(symbols, crds, pbc=False)
        self.assertEqual(frag_numb, 2)


if __name__ == "__main__":
    unittest.main()

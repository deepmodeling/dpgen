import dpdata
import numpy as np
import unittest

from context import take_cluster
from comp_sys import CompSys


class Test_take_cluster(unittest.TestCase, CompSys):
    def setUp (self) :
        type_map = ['C', 'H']
        jdata={
            "cluster_cutoff": 3.5
        }
        self.system_1 = take_cluster("cluster/14400.lammpstrj", type_map, 1125, jdata)
        self.system_2 = dpdata.LabeledSystem("cluster/input0_old.gaussianlog", fmt="gaussian/log")
        self.system_2.data['cells'] = self.system_1['cells']
        self.places=0


class Test_take_cluster_minify(unittest.TestCase, CompSys):
    def setUp (self) :
        type_map = ['C', 'H']
        jdata={
            "cluster_cutoff": 3.5,
            "cluster_minify": True
        }
        self.system_1 = take_cluster("cluster/14400.lammpstrj", type_map, 1125, jdata)
        self.system_2 = dpdata.LabeledSystem("cluster/input0_new.gaussianlog", fmt="gaussian/log")
        self.system_2.data['cells'] = self.system_1['cells']
        self.places=0

if __name__ == '__main__':
    unittest.main()

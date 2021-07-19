import os,sys
import dpdata
import numpy as np
import unittest
import importlib

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import check_cluster
from .context import setUpModule


class Test_check_cluster(unittest.TestCase):
    def test (self) :
        conf_name='POSCAR_Au_cluster'
        fmt='POSCAR'
        sys = dpdata.System(conf_name, fmt = fmt)
        ret = check_cluster(sys, fp_cluster_vacuum=15)
        #bad cluster
        self.assertTrue(ret)
        #good cluster 
        ret = check_cluster(sys, fp_cluster_vacuum=10)
        self.assertFalse(ret)

if __name__ == '__main__':
    unittest.main()
         



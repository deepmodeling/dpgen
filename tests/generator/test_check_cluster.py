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
        ret = check_cluster(conf_name, fp_cluster_vacuum=15, fmt=fmt)
        #bad cluster
        self.assertTrue(ret)
        #good cluster 
        ret = check_cluster(conf_name, fp_cluster_vacuum=10, fmt=fmt)
        self.assertFalse(ret)

if __name__ == '__main__':
    unittest.main()
         



import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import check_bad_box

class TestCheckBadBox(unittest.TestCase):
    def test_length_ratio(self):
        dirname = os.path.dirname(__file__)
        conf_bad = os.path.join(dirname, 'check_bad_box', 'bad.length.lammpstrj')
        conf_good = os.path.join(dirname, 'check_bad_box', 'good.lammpstrj')
        sys_bad = dpdata.System(conf_bad, fmt = 'lammps/dump')
        sys_good = dpdata.System(conf_good, fmt = 'lammps/dump')
        self.assertTrue(check_bad_box(sys_bad, 'length_ratio:5'))
        self.assertFalse(check_bad_box(sys_good, 'length_ratio:5'))

    def test_height_ratio(self):
        dirname = os.path.dirname(__file__)
        conf_bad = os.path.join(dirname, 'check_bad_box', 'bad.height.POSCAR')
        sys_bad = dpdata.System(conf_bad, fmt = 'vasp/POSCAR')
        self.assertTrue(check_bad_box(sys_bad, 'height_ratio:5'))
        self.assertFalse(check_bad_box(sys_bad, 'length_ratio:5'))

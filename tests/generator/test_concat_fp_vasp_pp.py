import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import sys_link_fp_vasp_pp
from .context import make_iter_name
from .context import fp_name
from .context import setUpModule

class TestConcatVASPPP(unittest.TestCase):
    def test(self):
        iter_index = 0
        iter_name = make_iter_name(iter_index)
        work_path = os.path.join(iter_name, fp_name)

        if os.path.isdir(iter_name) :
            shutil.rmtree(iter_name)
        os.makedirs(iter_name, exist_ok = False)
        os.makedirs(work_path, exist_ok = False)
        os.makedirs(os.path.join(work_path, 'task.000.000000'), exist_ok = False)
        os.makedirs(os.path.join(work_path, 'task.000.000001'), exist_ok = False)
        os.makedirs(os.path.join(work_path, 'task.001.000000'), exist_ok = False)
        os.makedirs(os.path.join(work_path, 'task.001.000001'), exist_ok = False)
        shutil.copyfile(os.path.join('vasp', 'POSCAR.oh'), 
                        os.path.join(work_path, 'task.000.000000', 'POSCAR'))
        shutil.copyfile(os.path.join('vasp', 'POSCAR.oh'), 
                        os.path.join(work_path, 'task.000.000001', 'POSCAR'))
        shutil.copyfile(os.path.join('vasp', 'POSCAR.ch4'), 
                        os.path.join(work_path, 'task.001.000000', 'POSCAR'))
        shutil.copyfile(os.path.join('vasp', 'POSCAR.ch4'), 
                        os.path.join(work_path, 'task.001.000001', 'POSCAR'))
        sys_link_fp_vasp_pp(0, {'type_map' : ['H', 'C', 'O'], 
                                'fp_pp_path': os.path.join('vasp', 'potcars'), 
                                'fp_pp_files': ['POTCAR.H', 'POTCAR.C', 'POTCAR.O'],
        })
        self.assertTrue(os.path.isfile(os.path.join(work_path, 'POTCAR.000')))
        self.assertTrue(os.path.isfile(os.path.join(work_path, 'POTCAR.001')))        
        with open((os.path.join(work_path, 'POTCAR.000'))) as fp:
            pot = fp.read()
            self.assertEqual(pot, 'O\nH\n')
        with open((os.path.join(work_path, 'POTCAR.001'))) as fp:
            pot = fp.read()
            self.assertEqual(pot, 'H\nC\n')
        for ii in ['task.000.000000', 'task.000.000001'] :
            with open(os.path.join(work_path, ii, 'POTCAR')) as fp:
                pot = fp.read()
                self.assertEqual(pot, 'O\nH\n')
        for ii in ['task.001.000000', 'task.001.000001'] :
            with open(os.path.join(work_path, ii, 'POTCAR')) as fp:
                pot = fp.read()
                self.assertEqual(pot, 'H\nC\n')
                


if __name__ == '__main__':
    unittest.main()



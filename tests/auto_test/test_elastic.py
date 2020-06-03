import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import make_kspacing_kpoints
from .context import setUpModule

from pymatgen.io.vasp import Incar
from dpgen.auto_test.gen_02_elastic import make_vasp

class Test02(unittest.TestCase):
    def tearDown(self):
        if os.path.exists('02.elastic'):
            shutil.rmtree('02.elastic')

    def test_make_vasp (self):
        jdata = {
            'relax_incar' : 'vasp_input/INCAR.rlx',
            'potcar_map': {'Si': 'vasp_input/POTCAR' },
        }
        make_vasp(jdata, 'confs/si/mp-149')
        
        target_path = '02.elastic/si/mp-149/vasp-relax_incar'
        equi_path = '00.equi/si/mp-149/vasp-relax_incar'
        dfm_dirs = glob.glob(os.path.join(target_path, 'dfm*'))
       
        # check root INCAR
        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar1 = Incar.from_file(os.path.join(target_path, 'INCAR'))
        self.assertFalse(incar0 == incar1)
        incar0['ISIF'] = 2
        self.assertTrue(incar0 == incar1)
        # check root POTCAR
        with open(os.path.join('vasp_input', 'POTCAR')) as fp:
            pot0 = fp.read()
        with open(os.path.join(target_path, 'POTCAR')) as fp:
            pot1 = fp.read()
        self.assertEqual(pot0, pot1)
        # check root POSCAR
        self.assertEqual(os.path.realpath(os.path.join(target_path, 'POSCAR')), 
                         os.path.realpath(os.path.join(equi_path, 'CONTCAR')))
        # check subdir
        for ii in dfm_dirs:
            self.assertEqual(os.path.realpath(os.path.join(ii, 'INCAR')),
                             os.path.realpath(os.path.join(target_path, 'INCAR')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'KPOINTS')),
                             os.path.realpath(os.path.join(target_path, 'KPOINTS')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'POTCAR')),
                             os.path.realpath(os.path.join(target_path, 'POTCAR')))

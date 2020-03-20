import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import make_kspacing_kpoints
from .context import setUpModule

from pymatgen.io.vasp import Incar
from dpgen.auto_test.gen_01_eos import make_vasp

class Test01(unittest.TestCase):
    def tearDown(self):
        if os.path.exists('01.eos'):
            shutil.rmtree('01.eos')

    def test_make_vasp_rlx_cell_shape (self):
        jdata = {
            'relax_incar' : 'vasp_input/INCAR.rlx',
            'potcar_map': {'Si': 'vasp_input/POTCAR' },
            'vol_start': 15,
            'vol_end':  25,
            'vol_step':  1.0,
            'eos_relax_cell_shape':  True,
        }
        make_vasp(jdata, 'confs/si/mp-149')
        
        target_path = '01.eos/si/mp-149/vasp-relax_incar'
        equi_path = '00.equi/si/mp-149/vasp-relax_incar'
        dfm_dirs = glob.glob(os.path.join(target_path, 'vol*'))
       
        # check root INCAR
        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar1 = Incar.from_file(os.path.join(target_path, 'INCAR'))
        self.assertFalse(incar0 == incar1)
        incar0['ISIF'] = 4
        self.assertTrue(incar0 == incar1)
        # check root POTCAR
        with open(os.path.join('vasp_input', 'POTCAR')) as fp:
            pot0 = fp.read()
        with open(os.path.join(target_path, 'POTCAR')) as fp:
            pot1 = fp.read()
        self.assertEqual(pot0, pot1)
        # check subdir
        for ii in dfm_dirs:
            self.assertFalse(os.path.isfile(os.path.join(ii, 'KPOINTS')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'POSCAR.orig')),
                             os.path.realpath(os.path.join(equi_path, 'CONTCAR')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'INCAR')),
                             os.path.realpath(os.path.join(target_path, 'INCAR')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'POTCAR')),
                             os.path.realpath(os.path.join(target_path, 'POTCAR')))
            sys = dpdata.System(os.path.join(ii, 'POSCAR'))
            vol = float(ii.split('/')[-1].split('-')[1])
            natoms = sys.get_natoms()
            self.assertAlmostEqual(vol, np.linalg.det(sys['cells'][0]) / natoms)


    def test_make_vasp_norlx_cell_shape (self):
        jdata = {
            'relax_incar' : 'vasp_input/INCAR.rlx',
            'potcar_map': {'Si': 'vasp_input/POTCAR' },
            'vol_start': 15,
            'vol_end':  25,
            'vol_step':  1.0,
            'eos_relax_cell_shape':  False,
        }
        make_vasp(jdata, 'confs/si/mp-149')
        
        target_path = '01.eos/si/mp-149/vasp-relax_incar'
        equi_path = '00.equi/si/mp-149/vasp-relax_incar'
        dfm_dirs = glob.glob(os.path.join(target_path, 'vol*'))
       
        # check root INCAR
        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar1 = Incar.from_file(os.path.join(target_path, 'INCAR'))
        self.assertFalse(incar0 == incar1)
        incar0['ISIF'] = 2
        self.assertTrue(incar0 == incar1)

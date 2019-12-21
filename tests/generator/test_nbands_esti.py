import os,sys
import dpdata
import numpy as np
import unittest
import importlib

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import NBandsEsti

class TestNBandsEsti(unittest.TestCase):    
    def test_predict(self):
        self.nbe = NBandsEsti(['out_data_nbands_esti/md.010000K',
                               'out_data_nbands_esti/md.020000K', 
                               'out_data_nbands_esti/md.040000K', 
                               'out_data_nbands_esti/md.080000K', 
                               'out_data_nbands_esti/md.160000K', 
        ])
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.010000K'), 72)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.020000K'), 83)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.040000K'), 112)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.080000K'), 195)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.160000K'), 429)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.240000K'), 732)

    def test_save_load(self):
        self.nbe2 = NBandsEsti(['out_data_nbands_esti/md.010000K',
                               'out_data_nbands_esti/md.020000K', 
                               'out_data_nbands_esti/md.040000K', 
                               'out_data_nbands_esti/md.080000K', 
                               'out_data_nbands_esti/md.160000K', 
        ])
        self.nbe2.save('tmp.log')
        self.nbe = NBandsEsti('tmp.log')
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.010000K'), 72)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.020000K'), 83)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.040000K'), 112)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.080000K'), 195)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.160000K'), 429)
        self.assertEqual(self.nbe.predict('out_data_nbands_esti/md.240000K'), 732)
        os.remove('tmp.log')

    def test_get_default_nbands(self):
        res = NBandsEsti._get_res('out_data_nbands_esti/md.020000K/')
        nb = NBandsEsti._get_default_nbands(res)
        self.assertEqual(nb, 66)

    def test_get_default_nbands(self):
        res = NBandsEsti._get_res('out_data_nbands_esti/mgal/')
        nb = NBandsEsti._get_default_nbands(res)
        self.assertEqual(nb, 124)

    def test_potcar_nvalence (self) :
        res = NBandsEsti._get_potcar_nvalence('out_data_nbands_esti/POTCAR.dbl')
        self.assertEqual(res, [10., 3.])

    def test_incar_ele_temp (self) :
        res = NBandsEsti._get_incar_ele_temp('out_data_nbands_esti/md.000300K/INCAR')
        self.assertAlmostEqual(res, 0.025851991011651636)

    def test_incar_nbands (self) :
        res = NBandsEsti._get_incar_nbands('out_data_nbands_esti/md.020000K/INCAR')
        self.assertEqual(res, 81)

    def test_get_res(self):
        res = NBandsEsti._get_res('out_data_nbands_esti/md.020000K/')
        ref = {
            'natoms': [32],
            'vol': 138.55418502346618,
            'nvalence': [3.],
            'ele_temp': 20000.0,
            'nbands': 81
        }
        self.assertEqual(res['natoms'], ref['natoms'])
        self.assertAlmostEqual(res['vol'], ref['vol'])
        self.assertAlmostEqual(res['nvalence'][0], ref['nvalence'][0])
        self.assertEqual(len(res['nvalence']), len(ref['nvalence']))
        self.assertAlmostEqual(res['ele_temp'], ref['ele_temp'], places = 1)
        self.assertEqual(res['nbands'], ref['nbands'])    

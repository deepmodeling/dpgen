import sys, os
import unittest
from dpgen.data.tools.create_random_disturb import create_disturbs_abacus_dev
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'data'
from .context import *

class TestPertAbacus(unittest.TestCase):
    def setUp(self):
        create_disturbs_abacus_dev(abacus_stru_file,1)
    def tearDown(self):
        if os.path.isfile('STRU.hcp1.abacus'):
            os.remove('STRU.hcp1.abacus')
    def test_stru(self):
        if os.path.isfile('STRU.hcp1.abacus'):
            stru1 = get_abacus_STRU('STRU.hcp1.abacus')
        stru0 = get_abacus_STRU(abacus_stru_file)
        self.assertEqual(stru0['atom_names'],stru1['atom_names'])
        self.assertEqual(stru0['atom_numbs'],stru1['atom_numbs'])
        self.assertEqual(stru0['atom_masses'],stru1['atom_masses'])
        #print(stru0['atom_numbs'],stru1['atom_numbs'])
        #print(stru0['coords'],stru1['coords'])
        #print(stru0['cells'],stru1['cells'])
        #print(stru0['atom_types'],stru1['atom_types'])
    def test_FileExist(self):
        self.assertTrue(os.path.isfile('STRU.hcp1.abacus'))

if __name__ == '__main__':
    unittest.main()

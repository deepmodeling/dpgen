import sys, os
import unittest
import numpy as np
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
        coords = stru1['coords']
        cells = stru1['cells']
        mindis = 100
        minm = minn = mini = minj = mink = 0
        maxdis = 0
        volume = np.linalg.det(np.array(cells))
        for m in range(len(coords)):
            x1,y1,z1 = coords[m]

            for n in range(m,len(coords)):
                for i in range(-1,2):
                    for j in range(-1,2):
                        for k in range(-1,2):
                            if m==n and i==0 and j==0 and k==0:continue
                            x2 = coords[n][0] + i * cells[0][0] + j * cells[1][0] + k * cells[2][0]
                            y2 = coords[n][1] + i * cells[0][1] + j * cells[1][1] + k * cells[2][1]
                            z2 = coords[n][2] + i * cells[0][2] + j * cells[1][2] + k * cells[2][2]

                            distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
                            if distance < mindis:
                                mindis = distance
        self.assertTrue(volume>0.0)
        self.assertTrue(mindis>0.01)

    def test_FileExist(self):
        self.assertTrue(os.path.isfile('STRU.hcp1.abacus'))

if __name__ == '__main__':
    unittest.main()

import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'data'
from .context import coll_vasp_md
from .context import out_dir_name
from .context import param_file
from .context import setUpModule

class TestCollVasp(unittest.TestCase):
    def setUp(self):
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        self.odir = out_dir_name(jdata)
        assert os.path.isdir('out_data_02_md'), 'out data for post fp vasp should exist'
        if os.path.isdir(self.odir) :
            shutil.rmtree(self.odir)
        shutil.copytree('out_data_02_md', self.odir)
        self.ref_coord = [[[0, 0, 0], [2.3, 2.3, 2.3]],
                          [[0, 0, 0], [2.2, 2.3, 2.4]]]
        self.ref_cell = [4.6 * np.eye(3), 4.6 * np.eye(3)]
        self.ref_at = [1, 1]
        self.ref_e = [-1.90811235, -1.89718546]
        self.ref_f = [[[ 0.      ,  0.      ,  0.      ], \
                       [-0.      , -0.      , -0.      ]],\
                      [[-0.110216,  0.      ,  0.110216], \
                       [ 0.110216, -0.      , -0.110216]]]
        self.ref_v = [[[ 1.50816698,  0.        , -0.        ], \
                       [ 0.        ,  1.50816698,  0.        ], \
                       [-0.        ,  0.        ,  1.50816795]],\
                      [[ 1.45208913,  0.        ,  0.03036584], \
                       [ 0.        ,  1.67640928,  0.        ], \
                       [ 0.03036584,  0.        ,  1.45208913]]]
        self.ref_coord = np.array(self.ref_coord)
        self.ref_cell = np.array(self.ref_cell)
        self.ref_at = np.array(self.ref_at, dtype = int)
        self.ref_e = np.array(self.ref_e)
        self.ref_f = np.array(self.ref_f)
        self.ref_v = np.array(self.ref_v)

    def tearDown(self):
        shutil.rmtree(self.odir)

    def test_coll(self):

        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['out_dir'] = self.odir
        coll_vasp_md(jdata)
        
        sys = dpdata.LabeledSystem(self.odir + '/02.md/sys-004/deepmd//', fmt = 'deepmd/raw')
        self.assertEqual(sys.get_nframes(), 2)
        
        if sys.data['coords'][0][1][0] < sys.data['coords'][1][1][0]:
            idx = [1, 0]
        else :
            idx = [0, 1]
        ref_coord = self.ref_coord[idx]
        ref_cell = self.ref_cell[idx]
        ref_e = self.ref_e[idx]
        ref_f = self.ref_f[idx]
        ref_v = self.ref_v[idx]
        ref_at = self.ref_at
            
        for ff in range(2) :
            self.assertAlmostEqual(ref_e[ff], sys.data['energies'][ff])
        for ii in range(2) :
            self.assertEqual(ref_at[ff], sys.data['atom_types'][ff])
        for ff in range(2) :
            for ii in range(2) :
                for dd in range(3) :
                    self.assertAlmostEqual(ref_coord[ff][ii][dd], 
                                           sys.data['coords'][ff][ii][dd])
                    self.assertAlmostEqual(ref_f[ff][ii][dd], 
                                           sys.data['forces'][ff][ii][dd])
        for ff in range(2):
            for ii in range(3) :
                for jj in range(3) :
                    self.assertAlmostEqual(ref_v[ff][ii][jj], 
                                           sys.data['virials'][ff][ii][jj], places = 5)
                    self.assertAlmostEqual(ref_cell[ff][ii][jj], 
                                           sys.data['cells'][ff][ii][jj])

if __name__ == '__main__':
    unittest.main()
        

    

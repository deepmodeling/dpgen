import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'data'
from .context import coll_abacus_md
from .context import out_dir_name
from .context import abacus_param_file
from .context import setUpModule

class TestCollAbacus(unittest.TestCase):
    def setUp(self):
        with open (abacus_param_file, 'r') as fp :
            jdata = json.load (fp)
        self.odir = out_dir_name(jdata)
        assert os.path.isdir('out_data_02_md_abacus'), 'out data for post fp vasp should exist'
        if os.path.isdir(self.odir) :
            shutil.rmtree(self.odir)
        shutil.copytree('out_data_02_md_abacus', self.odir)
        self.ref_coord = np.reshape(np.genfromtxt("abacus.out/coord.raw"), [8, 5, 3])
        self.ref_cell = np.reshape(np.genfromtxt("abacus.out/box.raw"), [8, 3, 3])
        self.ref_e = np.reshape(np.genfromtxt("abacus.out/energy.raw"), [8, ])
        self.ref_f = np.reshape(np.genfromtxt("abacus.out/force.raw"), [8, 5, 3])
        self.ref_v = np.reshape(np.genfromtxt("abacus.out/virial.raw"), [8, 3, 3])
    def tearDown(self):
        #print("escape.")
        shutil.rmtree(self.odir)

    def test_coll(self):

        with open (abacus_param_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['out_dir'] = self.odir
        print(os.getcwd())
        coll_abacus_md(jdata)
        
        sys = dpdata.LabeledSystem(self.odir + '/02.md/sys-0004-0001/deepmd//', fmt = 'deepmd/raw')
        self.assertEqual(sys.get_nframes(), 8)
            
        for ff in range(8) :
            self.assertAlmostEqual(self.ref_e[ff], sys.data['energies'][ff])
        for ff in range(8) :
            for ii in range(5) :
                for dd in range(3) :
                    self.assertAlmostEqual(self.ref_coord[ff][ii][dd], 
                                           sys.data['coords'][ff][ii][dd])
                    self.assertAlmostEqual(self.ref_f[ff][ii][dd], 
                                           sys.data['forces'][ff][ii][dd])
        for ff in range(8):
            for ii in range(3) :
                for jj in range(3) :
                    self.assertAlmostEqual(self.ref_v[ff][ii][jj], 
                                           sys.data['virials'][ff][ii][jj], places = 5)
                    self.assertAlmostEqual(self.ref_cell[ff][ii][jj], 
                                           sys.data['cells'][ff][ii][jj])
                    
if __name__ == '__main__':
    unittest.main()
        

    

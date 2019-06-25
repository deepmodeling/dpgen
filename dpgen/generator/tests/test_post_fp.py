import os,json,glob,shutil,filecmp
import dpdata
import numpy as np
import unittest

from context import post_fp
from context import post_fp_pwscf
from context import post_fp_vasp
from context import param_file
from context import param_old_file
from context import param_pwscf_file
from context import param_pwscf_old_file
from context import machine_file
from comp_sys import test_atom_names
from comp_sys import test_atom_types
from comp_sys import test_coord
from comp_sys import test_cell
from comp_sys import CompLabeledSys


class TestPostFPVasp(unittest.TestCase):
    def setUp(self):
        assert os.path.isdir('out_data_post_fp_vasp'), 'out data for post fp vasp should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_vasp', 'iter.000000')
        self.ref_coord = [[[0, 0, 0], [2.3, 2.3, 2.3]],
                          [[0, 0, 0], [2.2, 2.3, 2.4]]]
        self.ref_cell = [4.6 * np.eye(3), 4.6 * np.eye(3)]
        self.ref_at = [0, 0]
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
        shutil.rmtree('iter.000000')

    def test_post_fp_vasp_0(self):

        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp_vasp(0, jdata)
        
        sys = dpdata.LabeledSystem('iter.000000/02.fp/data.000/', fmt = 'deepmd/raw')
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


    def test_post_fp_vasp_1(self):

        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        
        sys = dpdata.LabeledSystem('iter.000000/02.fp/data.001/', fmt = 'deepmd/raw')
        self.assertEqual(sys.get_nframes(), 1)
        
        # if sys.data['coords'][0][1][0] < sys.data['coords'][1][1][0]:
        #     idx = [0]
        # else :
        idx = [1]
        ref_coord = self.ref_coord[idx]
        ref_cell = self.ref_cell[idx]
        ref_e = self.ref_e[idx]
        ref_f = self.ref_f[idx]
        ref_v = self.ref_v[idx]
        ref_at = self.ref_at
            
        for ff in range(1) :
            self.assertAlmostEqual(ref_e[ff], sys.data['energies'][ff])
        for ii in range(2) :
            self.assertEqual(ref_at[ff], sys.data['atom_types'][ff])
        for ff in range(1) :
            for ii in range(2) :
                for dd in range(3) :
                    self.assertAlmostEqual(ref_coord[ff][ii][dd], 
                                           sys.data['coords'][ff][ii][dd])
                    self.assertAlmostEqual(ref_f[ff][ii][dd], 
                                           sys.data['forces'][ff][ii][dd])
        for ff in range(1):
            for ii in range(3) :
                for jj in range(3) :
                    self.assertAlmostEqual(ref_v[ff][ii][jj], 
                                           sys.data['virials'][ff][ii][jj], places = 5)
                    self.assertAlmostEqual(ref_cell[ff][ii][jj], 
                                           sys.data['cells'][ff][ii][jj])


class TestPostFPPWSCF(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir('out_data_post_fp_pwscf'), 'out data for post fp pwscf should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_pwscf', 'iter.000000')
        with open (param_pwscf_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')
        

    

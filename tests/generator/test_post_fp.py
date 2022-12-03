import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import post_fp
from .context import post_fp_pwscf
from .context import post_fp_abacus_scf
from .context import post_fp_siesta
from .context import post_fp_vasp
from .context import post_fp_gaussian
from .context import post_fp_cp2k
from .context import param_file
from .context import param_old_file
from .context import param_pwscf_file
from .context import param_pwscf_old_file
from .context import param_abacus_post_file
from .context import param_siesta_file
from .context import param_gaussian_file
from .context import param_cp2k_file
from .context import param_amber_file
from .context import machine_file
from .context import setUpModule
from .comp_sys import test_atom_names
from .comp_sys import test_atom_types
from .comp_sys import test_coord
from .comp_sys import test_cell
from .comp_sys import CompLabeledSys
from .context import param_pwmat_file


class TestPostFPVasp(unittest.TestCase):
    def setUp(self):
        assert os.path.isdir('out_data_post_fp_vasp'), 'out data for post fp vasp should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_vasp', 'iter.000000')
        self.ref_coord = [[[0, 0, 0], [2.3, 2.3, 2.3]],
                          [[0, 0, 0], [2.2, 2.3, 2.4]]]
        self.ref_cell = [4.6 * np.eye(3), 4.6 * np.eye(3)]
        # type_map = ["Mg", "Al"], Al OUTCAR provided
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
        shutil.rmtree('iter.000000')

    def test_post_fp_vasp_0(self):

        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['use_ele_temp'] = 2
        post_fp_vasp(0, jdata, rfailed=0.3)

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

        self.assertTrue(os.path.isfile('iter.000000/02.fp/data.000/set.000/aparam.npy'))
        aparam = np.load('iter.000000/02.fp/data.000/set.000/aparam.npy')
        natoms = sys.get_natoms()
        self.assertEqual(natoms, 2)
        self.assertEqual(list(list(aparam)[0]), [0,0])
        self.assertEqual(list(list(aparam)[1]), [1,1])


    def test_post_fp_vasp_1(self):

        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['use_ele_temp'] = 1
        post_fp_vasp(0, jdata, rfailed=0.3)

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

        fparam = np.load('iter.000000/02.fp/data.001/set.000/fparam.npy')
        self.assertEqual(fparam.shape[0], 1)
        self.assertEqual(list(fparam), [100000])


    def test_post_fp_vasp_2(self):
        with open (param_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['use_ele_temp'] = 1
        with self.assertRaises(RuntimeError):
            post_fp_vasp(0, jdata)


class TestPostFPPWSCF(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir('out_data_post_fp_pwscf'), 'out data for post fp pwscf should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_pwscf', 'iter.000000')
        with open (param_pwscf_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')

class TestPostFPABACUS(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir('out_data_post_fp_abacus'), 'out data for post fp pwscf should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_abacus', 'iter.000000')
        with open (param_abacus_post_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')

class TestPostFPSIESTA(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir('out_data_post_fp_siesta'), 'out data for post fp siesta should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_siesta', 'iter.000000')
        with open (param_siesta_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')

class TestPostGaussian(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir('out_data_post_fp_gaussian'), 'out data for post fp gaussian should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_gaussian', 'iter.000000')
        with open (param_gaussian_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')

class TestPostCP2K(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5
        assert os.path.isdir('out_data_post_fp_cp2k'), 'out data for post fp gaussian should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_cp2k', 'iter.000000')
        with open (param_cp2k_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')


class TestPostFPPWmat(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 2
        assert os.path.isdir('out_data_post_fp_pwmat'), 'out data for post fp pwmat should exist'
        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        shutil.copytree('out_data_post_fp_pwmat', 'iter.000000')
        with open (param_pwmat_file, 'r') as fp :
            jdata = json.load (fp)
        post_fp(0, jdata)
        self.system_1 = dpdata.LabeledSystem('iter.000000/orig', fmt = 'deepmd/raw')
        self.system_2 = dpdata.LabeledSystem('iter.000000/02.fp/data.000', fmt = 'deepmd/raw')


class TestPostAmberDiff(unittest.TestCase, CompLabeledSys):
    def setUp(self):
        self.places = 5
        self.e_places = 5
        self.f_places = 5
        self.v_places = 5

        if os.path.isdir('iter.000000') :
            shutil.rmtree('iter.000000')
        ms = dpdata.MultiSystems(dpdata.LabeledSystem(os.path.join('data', 'deepmd'), fmt="deepmd/raw"))
        ms.to_deepmd_npy(os.path.join('iter.000000', '02.fp', 'task.000.000000', 'dataset'))
        self.system_1 = list(ms.systems.values())[0]
        with open (param_amber_file, 'r') as fp :
            jdata = json.load (fp)
        jdata['type_map'] = self.system_1.get_atom_names()
        post_fp(0, jdata)
        self.system_2 = list(dpdata.MultiSystems(type_map = jdata['type_map']).from_deepmd_raw('iter.000000/02.fp/data.000').systems.values())[0]


if __name__ == '__main__':
    unittest.main()

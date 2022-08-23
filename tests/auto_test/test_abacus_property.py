import os, sys, shutil, glob
import numpy as np
import unittest
from monty.serialization import loadfn
from dpgen.generator.lib import abacus_scf
from dpgen.auto_test.ABACUS import  ABACUS

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import setUpModule

from dpgen.auto_test.EOS import EOS
from dpgen.auto_test.Elastic import Elastic
from dpgen.auto_test.Vacancy import Vacancy
from dpgen.auto_test.Interstitial import Interstitial
from dpgen.auto_test.Surface import Surface
from dpgen.auto_test.Gamma import Gamma
from dpgen.auto_test.common_prop import make_property

class TestABACUS(unittest.TestCase):

    def setUp(self):
        self.jdata = {
            "structures": ["confs/fcc-Al"],
            "interaction": {
                "type": "abacus",
                "incar": "abacus_input/INPUT",
                "potcar_prefix": "abacus_input",
                "potcars": {"Al": "Al_ONCV_PBE-1.0.upf"},
                "orb_files": {"Al":"Al_gga_9au_100Ry_4s4p1d.orb"}
            }
        }

        self.conf_path = 'confs/fcc-Al'
        self.equi_path = 'confs/fcc-Al/relaxation/relax_task'
        self.source_path = 'equi/abacus'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)
        if not os.path.exists(os.path.join(self.equi_path,'OUT.ABACUS')):
            os.makedirs(os.path.join(self.equi_path,'OUT.ABACUS'))
        for ifile in ['INPUT','STRU']:
            if not os.path.exists(os.path.join(self.equi_path,ifile)):
                shutil.copy(os.path.join(self.source_path,ifile),os.path.join(self.equi_path,ifile))
        for ifile in ['running_cell-relax.log','STRU_ION_D']:
            if not os.path.exists(os.path.join(self.equi_path,'OUT.ABACUS',ifile)):    
                shutil.copy(os.path.join(self.source_path,ifile),os.path.join(self.equi_path,'OUT.ABACUS',ifile))
        shutil.copy(os.path.join(self.source_path,'cell-relax.json'),os.path.join(self.equi_path,'result.json'))

        self.confs = self.jdata["structures"]
        self.inter_param = self.jdata["interaction"]
        self.ABACUS = ABACUS(self.inter_param, os.path.join(self.conf_path, 'STRU'))

    def tearDown(self):
        if os.path.exists('confs/fcc-Al/relaxation'):
            shutil.rmtree('confs/fcc-Al/relaxation')
        if os.path.exists('confs/fcc-Al/eos_00'):
            shutil.rmtree('confs/fcc-Al/eos_00')
        if os.path.exists('confs/fcc-Al/eos_02'):
            shutil.rmtree('confs/fcc-Al/eos_02')
        if os.path.exists('confs/fcc-Al/elastic_00'):
            shutil.rmtree('confs/fcc-Al/elastic_00')
        if os.path.exists('confs/fcc-Al/vacancy_00'):
            shutil.rmtree('confs/fcc-Al/vacancy_00')
        if os.path.exists('confs/fcc-Al/interstitial_00'):
            shutil.rmtree('confs/fcc-Al/interstitial_00')
        if os.path.exists('confs/fcc-Al/surface_00'):
            shutil.rmtree('confs/fcc-Al/surface_00')       

    def test_make_property(self):
        property = {"type":         "eos",
                    "vol_start":    0.85,
                    "vol_end":      1.15,
                    "vol_step":     0.01
                    }
        make_property(self.jdata["structures"], self.jdata["interaction"], [property])
        self.assertTrue(os.path.exists(os.path.join(self.conf_path,"eos_00")))
        self.assertTrue(os.path.exists(os.path.join(self.conf_path,"eos_00","INPUT")))
        for ii in glob.glob(os.path.join(self.conf_path,"eos_00", 'task.*')):
            self.assertTrue(os.path.exists(os.path.join(ii,"INPUT")))
            self.assertTrue(os.path.exists(os.path.join(ii,"pp_orb")))
            self.assertTrue(os.path.exists(os.path.join(ii,"KPT")))
            self.assertTrue(os.path.exists(os.path.join(ii,"STRU")))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'pp_orb', 'Al_ONCV_PBE-1.0.upf')),
                             os.path.realpath(os.path.join(self.jdata['interaction']['potcar_prefix'], 'Al_ONCV_PBE-1.0.upf')))

    def test_make_property_eos(self):
        property = {"type":         "eos",
                    "vol_start":    0.85,
                    "vol_end":      1.15,
                    "vol_step":     0.01
                    }
        work_path = os.path.join(self.conf_path,"eos_00")
        eos = EOS(property,self.inter_param)
        eos.make_confs(work_path, self.equi_path, refine=False)

        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'eos.json')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'STRU.orig')),
                             os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D')))

            eos_json = loadfn(os.path.join(ii, 'eos.json'))
            stru_data = abacus_scf.get_abacus_STRU(os.path.realpath(os.path.join(ii, 'STRU')))  
            vol_per_atom = abs(np.linalg.det(stru_data['cells'])) / np.array(stru_data['atom_numbs']).sum()
            self.assertAlmostEqual(eos_json['volume'], vol_per_atom)

    def test_make_property_elastic(self):
        property = {"type":         "elastic",
                    "norm_deform":  1e-2,
                    "shear_deform": 1e-2
                    }
        work_path = os.path.join(self.conf_path,"elastic_00")
        elastic = Elastic(property,self.inter_param)
        elastic.make_confs(work_path, self.equi_path, refine=False)

        self.assertEqual(os.path.realpath(os.path.join(work_path, 'STRU')),
                         os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D')))
        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'strain.json')))  

        os.remove(os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D')))
        with self.assertRaises(RuntimeError):
            elastic.make_confs(work_path, self.equi_path, refine=False)

    def test_make_property_elastic_post_process(self):
        property = {"type":         "elastic",
                    "norm_deform":  1e-2,
                    "shear_deform": 1e-2
                    }
        make_property(self.jdata["structures"], self.jdata["interaction"], [property])
        work_path = os.path.join(self.conf_path,"elastic_00")

        self.assertTrue(os.path.exists(os.path.join(work_path,"INPUT")))
        self.assertTrue(os.path.exists(os.path.join(work_path,"KPT")))

        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertEqual(os.path.realpath(os.path.join(work_path, 'KPT')),
                             os.path.realpath(os.path.join(ii, 'KPT')))
            self.assertEqual(os.path.realpath(os.path.join(work_path, 'INPUT')),
                             os.path.realpath(os.path.join(ii, 'INPUT')))  
 
    def test_make_property_vacancy(self):
        property = {"type":         "vacancy",
                    "supercell":    [1, 1, 1]
                    }
        work_path = os.path.join(self.conf_path,"vacancy_00")
        vacancy = Vacancy(property,self.inter_param)
        vacancy.make_confs(work_path, self.equi_path, refine=False)

        self.assertEqual(os.path.realpath(os.path.join(work_path, 'STRU')),
                         os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D')))

        stru_data = abacus_scf.get_abacus_STRU(os.path.realpath(os.path.join(work_path, 'STRU')))
        natom1 = np.array(stru_data['atom_numbs']).sum()
        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            stru_data = abacus_scf.get_abacus_STRU(os.path.realpath(os.path.join(ii, 'STRU'))) 
            natom2 = np.array(stru_data['atom_numbs']).sum()
            self.assertTrue(natom1==natom2+1)

    def test_make_property_interstitial(self):
        property = {"type":         "interstitial",
                    "supercell": [1, 1, 1],
                    "insert_ele": ["H"]
                    }
        self.inter_param['potcars']['H'] = 'H_ONCV_PBE-1.0.upf'
        self.inter_param['orb_files']['H'] = 'H_gga_8au_100Ry_2s1p.orb'

        work_path = os.path.join(self.conf_path,"interstitial_00")
        if os.path.exists(work_path):
            shutil.rmtree(work_path) 
        os.makedirs(work_path)
        interstitial = Interstitial(property,self.inter_param)
        interstitial.make_confs(work_path, self.equi_path, refine=False)  

        self.assertEqual(os.path.realpath(os.path.join(work_path, 'STRU')),
                         os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D'))) 
        stru_data = abacus_scf.get_abacus_STRU(os.path.realpath(os.path.join(work_path, 'STRU')))
        natom1 = np.array(stru_data['atom_numbs']).sum()         
        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            stru_data = abacus_scf.get_abacus_STRU(os.path.realpath(os.path.join(ii, 'STRU'))) 
            self.assertTrue('H' in stru_data['atom_names'])
            natom2 = np.array(stru_data['atom_numbs']).sum()
            self.assertTrue(natom1==natom2-1)

    def test_make_property_surface(self):
        property = {"type":         "surface",
                    "min_slab_size":  15,
                    "min_vacuum_size":11,
                    "pert_xz":        0.01,
                    "max_miller":     2,
                    "cal_type":       "static"
                    }
        work_path = os.path.join(self.conf_path,"surface_00")
        surface = Surface(property,self.inter_param)
        surface.make_confs(work_path, self.equi_path, refine=False)

        self.assertEqual(os.path.realpath(os.path.join(work_path, 'STRU')),
                         os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D')))
        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'miller.json')))       

    def test_make_property_gamma(self):
        property = {"type": "gamma",
                    "lattice_type": "fcc",
                    "miller_index": [1, 1, 1],
                    "displace_direction": [1, 1, 0],
                    "supercell_size": [1, 1, 10],
                    "min_vacuum_size": 10,
                    "add_fix": ["true", "true", "false"],
                    "n_steps": 20
                    }
        work_path = os.path.join(self.conf_path,"gamma_00")
        gamma = Gamma(property,self.inter_param)
        gamma.make_confs(work_path, self.equi_path, refine=False)

        dfm_dirs = glob.glob(os.path.join(work_path, 'task.*'))
        self.assertEqual(len(dfm_dirs), gamma.n_steps+1)

        self.assertEqual(os.path.realpath(os.path.join(work_path, 'STRU')),
                         os.path.realpath(os.path.join(self.equi_path, 'OUT.ABACUS', 'STRU_ION_D')))
        for ii in glob.glob(os.path.join(work_path, 'task.*')):
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'miller.json')))          

    def test_make_property_refine(self):
        property = {"type":         "eos",
                    "vol_start":    0.85,
                    "vol_end":      1.15,
                    "vol_step":     0.01
                    }
        pwd=os.getcwd()
        target_path_0 = 'confs/fcc-Al/eos_00'
        target_path_2 = 'confs/fcc-Al/eos_02'
        path_to_work = os.path.abspath(target_path_0)

        make_property(self.jdata["structures"], self.jdata["interaction"], [property])
        dfm_dirs_0 = glob.glob(os.path.join(target_path_0, 'task.*'))
        for ii in dfm_dirs_0:
            self.assertTrue(os.path.isfile(os.path.join(ii, 'STRU')))
            os.makedirs(os.path.join(ii,'OUT.ABACUS'))
            shutil.copy(os.path.join(ii, 'STRU'),os.path.join(ii, 'OUT.ABACUS', 'STRU_ION_D'))
        
        new_prop_list=[
        {
         "type":         "eos",
         "init_from_suffix": "00",
         "output_suffix": "02",
         "cal_setting": {
                           "relax_pos": True,
                            "relax_shape": True,
                            "relax_vol": False}
        }
        ]
        make_property(self.jdata["structures"], self.jdata["interaction"], new_prop_list)
        self.assertTrue(os.path.isdir(path_to_work.replace('00','02')))
        os.chdir(pwd)
        dfm_dirs_2 = glob.glob(os.path.join(target_path_2, 'task.*'))
        self.assertEqual(len(dfm_dirs_2),len(dfm_dirs_0))
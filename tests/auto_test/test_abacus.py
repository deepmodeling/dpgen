import os, sys, shutil
import numpy as np
import unittest
from monty.serialization import loadfn
from dpgen.generator.lib import abacus_scf
from dpgen.auto_test.lib import  abacus
from dpgen.auto_test.ABACUS import  ABACUS

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import setUpModule

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
            },
            "relaxation": {
                         "cal_type":     "relaxation",
                         "cal_setting":  {"relax_pos":True,
                                          "relax_shape":True,
                                          "relax_vol":True}
            }
        }

        self.conf_path = 'confs/fcc-Al'
        self.equi_path = 'confs/fcc-Al/relaxation/relax_task'
        self.source_path = 'equi/abacus'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        self.task_param = self.jdata["relaxation"]
        self.ABACUS = ABACUS(inter_param, os.path.join(self.conf_path, 'STRU'))

    def tearDown(self):
        if os.path.exists('confs/fcc-Al/relaxation'):
            shutil.rmtree('confs/fcc-Al/relaxation')

    def test_make_potential_files(self):
        if not os.path.exists(os.path.join(self.equi_path, 'STRU')):
            with self.assertRaises(FileNotFoundError):
               self.ABACUS.make_potential_files(self.equi_path)
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        self.ABACUS.make_potential_files(self.equi_path)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "pp_orb/Al_ONCV_PBE-1.0.upf")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "pp_orb/Al_gga_9au_100Ry_4s4p1d.orb")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, 'inter.json')))

    def test_make_input_file_1(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":True,
                         "relax_shape":True,
                         "relax_vol":False}
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))                 
        self.ABACUS.make_input_file(self.equi_path,'relaxation',param)       
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "task.json")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "KPT"))) 
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "INPUT")))
        abacus_input = abacus_scf.get_abacus_input_parameters(os.path.join(self.equi_path, "INPUT"))
        self.assertEqual(abacus_input['calculation'].lower(),"cell-relax")
        self.assertEqual(abacus_input['fixed_axes'].lower(),"volume")
        self.assertTrue(abacus.check_stru_fixed(os.path.join(self.equi_path, 'STRU'),fixed=False))

    def test_make_input_file_2(self):
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        self.ABACUS.make_input_file(self.equi_path,'relaxation',self.task_param)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "task.json")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "KPT")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "INPUT")))
        abacus_input = abacus_scf.get_abacus_input_parameters(os.path.join(self.equi_path, "INPUT"))   
        self.assertEqual(abacus_input['calculation'].lower(),"cell-relax")
        self.assertTrue('fixed_axes' not in abacus_input or abacus_input['fixed_axes'] == 'None')
        self.assertTrue(abacus.check_stru_fixed(os.path.join(self.equi_path, 'STRU'),fixed=False))

    def test_make_input_file_3(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":True,
                         "relax_shape":False,
                         "relax_vol":False}
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        self.ABACUS.make_input_file(self.equi_path,'relaxation',param)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "task.json")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "KPT")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "INPUT")))
        abacus_input = abacus_scf.get_abacus_input_parameters(os.path.join(self.equi_path, "INPUT"))   
        self.assertEqual(abacus_input['calculation'].lower(),"relax")
        self.assertTrue('fixed_axes' not in abacus_input or abacus_input['fixed_axes'] == 'None')
        self.assertTrue(abacus.check_stru_fixed(os.path.join(self.equi_path, 'STRU'),fixed=False))

    def test_make_input_file_4(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":False,
                         "relax_shape":True,
                         "relax_vol":False}
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        self.ABACUS.make_input_file(self.equi_path,'relaxation',param)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "task.json")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "KPT")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "INPUT")))
        abacus_input = abacus_scf.get_abacus_input_parameters(os.path.join(self.equi_path, "INPUT"))   
        self.assertEqual(abacus_input['calculation'].lower(),"cell-relax")
        self.assertEqual(abacus_input['fixed_axes'].lower(),"volume")
        self.assertTrue(abacus.check_stru_fixed(os.path.join(self.equi_path, 'STRU'),fixed=True))

    def test_make_input_file_5(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":False,
                         "relax_shape":True,
                         "relax_vol":True}
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        self.ABACUS.make_input_file(self.equi_path,'relaxation',param)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "task.json")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "KPT")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "INPUT")))
        abacus_input = abacus_scf.get_abacus_input_parameters(os.path.join(self.equi_path, "INPUT"))   
        self.assertEqual(abacus_input['calculation'].lower(),"cell-relax")
        self.assertTrue('fixed_axes' not in abacus_input or abacus_input['fixed_axes'] == 'None')
        self.assertTrue(abacus.check_stru_fixed(os.path.join(self.equi_path, 'STRU'),fixed=True))

    def test_make_input_file_kspacing(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":False,
                         "relax_shape":True,
                         "relax_vol":True,
                         "kspacing":0.1}
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        self.ABACUS.make_input_file(self.equi_path,'relaxation',param)
        with open(os.path.join(self.equi_path, "KPT")) as f1: kpt = f1.read().strip().split('\n')[-1].split()
        self.assertEqual(kpt,['9','9','9','0','0','0'])

    def test_compuate(self):
        ret=self.ABACUS.compute(os.path.join(self.equi_path))
        self.assertIsNone(ret)
        shutil.copy(os.path.join(self.source_path, 'INPUT'), os.path.join(self.equi_path, 'INPUT'))
        shutil.copy(os.path.join(self.conf_path, 'STRU'), os.path.join(self.equi_path, 'STRU'))
        os.mkdir(os.path.join(self.equi_path, 'OUT.ABACUS'))
        shutil.copy(os.path.join(self.source_path, 'running_cell-relax.log'), os.path.join(self.equi_path, 'OUT.ABACUS','running_cell-relax.log'))
        ret=self.ABACUS.compute(self.equi_path)
        ret_ref=loadfn(os.path.join(self.source_path, 'cell-relax.json'))

        def compare_dict(dict1,dict2):
            self.assertEqual(dict1.keys(),dict2.keys())
            for key in dict1:
                if type(dict1[key]) is dict:
                    compare_dict(dict1[key],dict2[key])
                else:
                    if type(dict1[key]) is np.ndarray:
                        np.testing.assert_almost_equal(dict1[key], dict2[key],decimal=5)
                    else:
                        self.assertTrue(dict1[key] == dict2[key])

        compare_dict(ret,ret_ref.as_dict())
        

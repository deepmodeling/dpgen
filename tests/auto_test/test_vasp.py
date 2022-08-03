import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
from dpdata import LabeledSystem
from pymatgen.io.vasp import Incar
from monty.serialization import loadfn, dumpfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.VASP import  VASP
from dpgen.generator.lib.vasp import incar_upper


class TestVASP(unittest.TestCase):

    def setUp(self):
        self.jdata = {
            "structures": ["confs/hp-*"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR",
                "potcar_prefix": ".",
                "potcars": {"Li": "vasp_input/POTCAR"}
            },
            "relaxation": {
                         "cal_type":     "relaxation",
                         "cal_setting":  {"relax_pos":True,
                                          "relax_shape":True,
                                          "relax_vol":True}
            }
        }

        self.conf_path = 'confs/hp-Li'
        self.equi_path = 'confs/hp-Li/relaxation/relax_task'
        self.source_path = 'equi/vasp'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        self.task_param = self.jdata["relaxation"]
        self.VASP = VASP(inter_param, os.path.join(self.conf_path, 'POSCAR'))

    def tearDown(self):
        if os.path.exists('confs/hp-Li/relaxation'):
            shutil.rmtree('confs/hp-Li/relaxation')
        if os.path.exists('inter.json'):
            os.remove('inter.json')
        if os.path.exists('POTCAR'):
            os.remove('POTCAR')

    def test_make_potential_files(self):
        if not os.path.exists(os.path.join(self.equi_path, 'POSCAR')):
            with self.assertRaises(FileNotFoundError):
               self.VASP.make_potential_files(self.equi_path)
        shutil.copy(os.path.join(self.conf_path, 'POSCAR'), os.path.join(self.equi_path, 'POSCAR'))
        self.VASP.make_potential_files(self.equi_path)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "POTCAR")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, 'inter.json')))

    def test_make_input_file_1(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":True,
                         "relax_shape":True,
                         "relax_vol":False}
        self.VASP.make_input_file(self.equi_path,'relaxation',param)
        incar=incar_upper(Incar.from_file(os.path.join(self.equi_path, "INCAR")))
        self.assertEqual(incar['ISIF'],4)

    def test_make_input_file_2(self):
        self.VASP.make_input_file(self.equi_path,'relaxation',self.task_param)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "task.json")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "KPOINTS")))   
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, "INCAR")))   
        incar=incar_upper(Incar.from_file(os.path.join(self.equi_path, "INCAR")))
        self.assertEqual(incar['ISIF'],3)

    def test_make_input_file_3(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":True,
                         "relax_shape":False,
                         "relax_vol":False}
        self.VASP.make_input_file(self.equi_path,'relaxation',param)
        incar=incar_upper(Incar.from_file(os.path.join(self.equi_path, "INCAR")))
        self.assertEqual(incar['ISIF'],2)

    def test_make_input_file_4(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":False,
                         "relax_shape":True,
                         "relax_vol":False}
        self.VASP.make_input_file(self.equi_path,'relaxation',param)
        incar=incar_upper(Incar.from_file(os.path.join(self.equi_path, "INCAR")))
        self.assertEqual(incar['ISIF'],5)

    def test_make_input_file_5(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":False,
                         "relax_shape":True,
                         "relax_vol":True}
        self.VASP.make_input_file(self.equi_path,'relaxation',param)
        incar=incar_upper(Incar.from_file(os.path.join(self.equi_path, "INCAR")))
        self.assertEqual(incar['ISIF'],6)

    def test_make_input_file_5(self):
        param=self.task_param.copy()
        param["cal_setting"]= {"relax_pos":False,
                         "relax_shape":True,
                         "relax_vol":True,
                         "kspacing":0.01}
        self.VASP.make_input_file(self.equi_path,'relaxation',param)
        incar=incar_upper(Incar.from_file(os.path.join(self.equi_path, "INCAR")))
        self.assertEqual(incar['ISIF'],6)
        self.assertEqual(incar['KSPACING'],0.01)

    def test_compuate(self):
        ret=self.VASP.compute(os.path.join(self.conf_path,'relaxation'))
        self.assertIsNone(ret)
        shutil.copy(os.path.join(self.source_path, 'OUTCAR'), os.path.join(self.equi_path, 'OUTCAR'))
        ret=self.VASP.compute(self.equi_path)
        ret_ref=loadfn(os.path.join(self.source_path, 'outcar.json'))
        
        def compare_dict(dict1,dict2):
            self.assertEqual(dict1.keys(),dict2.keys())
            for key in dict1:
                if key == 'stress':
                    self.assertTrue((np.array(dict1[key]['data']) == dict2[key]).all())
                elif type(dict1[key]) is dict:
                    compare_dict(dict1[key],dict2[key])
                else:
                    if type(dict1[key]) is np.ndarray:
                        self.assertTrue((dict1[key] == dict2[key]).all())
                    else:
                        self.assertTrue(dict1[key] == dict2[key])

        compare_dict(ret,ret_ref.as_dict())
        
        
    def test_backward_files(self):
        backward_files = ['OUTCAR', 'outlog', 'CONTCAR', 'OSZICAR', 'XDATCAR']
        self.assertEqual(self.VASP.backward_files(), backward_files)

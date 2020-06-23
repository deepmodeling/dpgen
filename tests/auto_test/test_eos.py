import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest
from monty.serialization import loadfn,dumpfn
from pymatgen.io.vasp import Incar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.EOS import EOS

class TestEOS(unittest.TestCase):

    def setUp(self):
        _jdata={
        "structures":    ["confs/hp-*"],
        "interaction": {
            "type":      "vasp",
            "incar":     "vasp_input/INCAR.rlx",
            "potcar_prefix":".",
            "type_map": ['Li'],
            "potcars":  {"Li": "vasp_input/POTCAR"}
        },
        "properties": [
            {
             "type":         "eos",
             "vol_start":    10,
             "vol_end":      30,
             "vol_step":     3,
             "change_box":   True
            }
        ]
        }
                 
        self.equi_path =   'confs/hp-Li/relaxation'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/hp-Li/eos_00'
        if not os.path.exists(self.equi_path):
           os.mkdir(self.equi_path)

        self.confs=_jdata["structures"]
        self.inter_param=_jdata["interaction"]
        self.prop_param=_jdata['properties']
        
        self.eos=EOS(_jdata['properties'][0])

    def tearDown(self):
        if os.path.exists('confs/hp-Li/relaxation'):
            shutil.rmtree('confs/hp-Li/relaxation')
        if os.path.exists('frozen_model.pb'):
           os.remove('frozen_model.pb')
        if os.path.exists('inter.json'):
           os.remove('inter.json')

    def test_task_type (self):
        self.assertEqual('eos',self.eos.task_type() )

    def test_task_param (self):
        self.assertEqual(self.prop_param[0],self.eos.task_param() )


    def test_make_confs(self):
        if not os.path.exists(os.path.join(self.equi_path,'CONTCAR')):
           with self.assertRaises(RuntimeError):
                self.eos.make_confs(self.target_path,self.equi_path)
        shutil.copy(os.path.join(self.source_path,'CONTCAR'),os.path.join(self.equi_path,'CONTCAR'))
        task_list=self.eos.make_confs(self.target_path,self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, 'task.*'))

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar0['ISIF'] = 4

        with open(os.path.join('vasp_input', 'POTCAR')) as fp:
             pot0 = fp.read()
        for ii in dfm_dirs:
            self.assertTrue(os.path.isfile(os.path.join(ii, 'POSCAR')))
            self.assertEqual(os.path.realpath(os.path.join(ii, 'POSCAR.orig')),
                             os.path.realpath(os.path.join(self.equi_path, 'CONTCAR')))

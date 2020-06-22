import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest
from monty.serialization import loadfn,dumpfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.VASP import VASP

class TestVASP(unittest.TestCase):

    def setUp(self):
        self.jdata= {
                      "structures":    ["confs/hp-*"],
                      "interaction": {
                          "type":      "vasp",
                          "incar":     "vasp_input/INCAR",
                          "potcar_prefix":".",
                          "potcars":    {"Li": "vasp_input/POTCAR"}
                      },
                      "relaxation": {
                                 "ediff": 1e-7,
                                 "ediffg": -0.01,
                                 "encut": 650,
                                 "kspacing": 0.1,
                                 "kgamma": False
                      }
                  }
                 
        self.equi_path =   'confs/hp-Li/relaxation'
        self.source_path = 'equi/vasp'
        if not os.path.exists(self.equi_path):
           os.mkdir(self.equi_path)

        self.confs=self.jdata["structures"]
        inter_param=self.jdata["interaction"]
        self.VASP=VASP(inter_param,os.path.join(self.source_path,'POSCAR'))

    def tearDown(self):
        if os.path.exists('confs/hp-Li/relaxation'):
            shutil.rmtree('confs/hp-Li/relaxation')
        if os.path.exists('inter.json'):
           os.remove('inter.json')
        if os.path.exists('POTCAR'):
           os.remove('POTCAR')

    def test_make_potential_files(self):
        self.VASP.make_potential_files(".")
        self.assertTrue(os.path.isfile("POTCAR"))
        self.assertTrue(os.path.isfile('inter.json'))

        self.VASP.make_potential_files(self.equi_path)
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path,"POTCAR")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path,'inter.json')))
 

    def test_backward_files(self):
        backward_files=['OUTCAR', 'outlog', 'CONTCAR', 'OSZICAR']
        self.assertEqual(self.VASP.backward_files(),backward_files)

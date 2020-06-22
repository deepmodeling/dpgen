import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest
from monty.serialization import loadfn,dumpfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.lib.lammps import inter_deepmd
from dpgen.auto_test.Lammps import Lammps

class TestLammps(unittest.TestCase):

    def setUp(self):
        self.jdata={
                     "structures":    ["confs/hp-*"],
                     "interaction": {
                         "type":      "deepmd",
                         "model":     "lammps_input/frozen_model.pb",
                         "deepmd_version":"1.1.0",
                         "type_map":    {"Li": 0}
                     },
                     "relaxation": {
                         "etol": 1e-12,
                         "ftol": 1e-6,
                         "maxiter": 5000,
                         "maximal": 500000,
                         "change_box": True
                     }
                 }
                 
        self.equi_path =   'confs/hp-Li/relaxation'
        self.source_path = 'equi/vasp'
        if not os.path.exists(self.equi_path):
           os.mkdir(self.equi_path)

        self.confs=self.jdata["structures"]
        inter_param=self.jdata["interaction"]
        self.Lammps=Lammps(inter_param,os.path.join(self.source_path,'POSCAR'))

    def tearDown(self):
        if os.path.exists('confs/hp-Li/relaxation'):
            shutil.rmtree('confs/hp-Li/relaxation')
        if os.path.exists('frozen_model.pb'):
           os.remove('frozen_model.pb')
        if os.path.exists('inter.json'):
           os.remove('inter.json')

    def test_set_inter_type_func (self):
        self.Lammps.set_inter_type_func()
        self.assertEqual(inter_deepmd,self.Lammps.inter_func )

    def test_set_model_param(self):
        self.Lammps.set_model_param()
        model_param = {'model_name': ['frozen_model.pb'],
                       'param_type':  {"Li": 0},
                       'deepmd_version': '1.1.0'}
        self.assertEqual(model_param,self.Lammps.model_param)

    def test_make_potential_files(self):
        self.Lammps.make_potential_files(".")
        self.assertTrue(os.path.islink("frozen_model.pb"))
        self.assertTrue(os.path.isfile('inter.json'))

        self.Lammps.make_potential_files(self.equi_path)
        self.assertTrue(os.path.islink(os.path.join(self.equi_path,"frozen_model.pb")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path,'inter.json')))

    def test_backward_files(self):
        backward_files=['log.lammps', 'outlog', 'dump.relax']
        self.assertEqual(self.Lammps.backward_files(),backward_files)

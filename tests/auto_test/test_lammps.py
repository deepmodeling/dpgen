import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
from monty.serialization import loadfn, dumpfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.lib.lammps import inter_deepmd
from dpgen.auto_test.Lammps import Lammps


class TestLammps(unittest.TestCase):

    def setUp(self):
        self.jdata = {
            "structures": ["confs/std-fcc"],
            "interaction": {
                "type": "deepmd",
                "model": "lammps_input/frozen_model.pb",
                "deepmd_version": "1.1.0",
                "type_map": {"Al": 0}
            },
            "relaxation": {
                "cal_type": "relaxation",
                "cal_setting": {"relax_pos": True,
                                "relax_shape": True,
                                "relax_vol": True}
            }
        }

        self.equi_path = 'confs/std-fcc/relaxation/relax_task'
        self.source_path = 'equi/lammps'

        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        if not os.path.isfile(os.path.join(self.equi_path,'POSCAR')):
           shutil.copy(os.path.join(self.source_path, 'Al-fcc.vasp'), os.path.join('confs/std-fcc','POSCAR'))
             

        self.confs = self.jdata["structures"]
        self.inter_param = self.jdata["interaction"]
        self.relax_param = self.jdata["relaxation"]
        self.Lammps = Lammps(self.inter_param, os.path.join(self.source_path, 'Al-fcc.vasp'))

    def tearDown(self):
        if os.path.exists('confs/std-fcc/relaxation'):
            shutil.rmtree('confs/std-fcc/relaxation')

    def test_set_inter_type_func(self):
        self.Lammps.set_inter_type_func()
        self.assertEqual(inter_deepmd, self.Lammps.inter_func)

    def test_set_model_param(self):
        self.Lammps.set_model_param()
        model_param = {'model_name': ['frozen_model.pb'],
                       'param_type': {"Al": 0},
                       'deepmd_version': '1.1.0'}
        self.assertEqual(model_param, self.Lammps.model_param)

    def test_make_potential_files(self):
        cwd=os.getcwd()
        abs_equi_path=os.path.abspath(self.equi_path)
        self.Lammps.make_potential_files(abs_equi_path)
        self.assertTrue(os.path.islink(os.path.join(self.equi_path, "frozen_model.pb")))
        self.assertTrue(os.path.isfile(os.path.join(self.equi_path, 'inter.json')))
        ret=loadfn(os.path.join(self.equi_path, 'inter.json'))
        self.assertTrue(self.inter_param,ret)
        os.chdir(cwd)

    def test_make_input_file(self):
        cwd=os.getcwd()
        abs_equi_path=os.path.abspath('confs/std-fcc/relaxation/relax_task')
        shutil.copy(os.path.join('confs/std-fcc','POSCAR'), os.path.join(self.equi_path, 'POSCAR'))
        self.Lammps.make_input_file(abs_equi_path,'relaxation', self.relax_param)
        self.assertTrue(os.path.isfile(os.path.join(abs_equi_path, "conf.lmp")))
        self.assertTrue(os.path.islink(os.path.join(abs_equi_path, "in.lammps")))
        self.assertTrue(os.path.isfile(os.path.join(abs_equi_path, "task.json")))

    def test_forward_common_files(self):
        fc_files = ['in.lammps', 'frozen_model.pb']
        self.assertEqual(self.Lammps.forward_common_files(), fc_files)

    def test_backward_files(self):
        backward_files = ['log.lammps', 'outlog', 'dump.relax']
        self.assertEqual(self.Lammps.backward_files(), backward_files)

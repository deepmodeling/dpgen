import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
import dpdata
from monty.serialization import loadfn, dumpfn
from pymatgen.io.vasp import Incar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.EOS import EOS


class TestEOS(unittest.TestCase):

    def setUp(self):
        _jdata = {
            "structures": ["confs/std-fcc"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": ".",
                "potcars": {"Li": "vasp_input/POTCAR"}
            },
            "properties": [
                 {
                  "type":         "eos",
                  "skip":  False,
                  "vol_start":    0.8,
                  "vol_end":      1.2,
                  "vol_step":     0.01,
                  "cal_setting": {
                                  "relax_pos": True, 
                                  "relax_shape": True, 
                                  "relax_vol": False, 
                                  "overwrite_interaction":{
                                              "type":    "deepmd", 
                                              "model":   "lammps_input/frozen_model.pb", 
                                              "type_map":{"Al": 0}
                                              }
                                 }
                 }
            ]
        }

        self.equi_path = 'confs/std-fcc/relaxation/relax_task'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/std-fcc/eos_00'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata['properties']

        self.eos = EOS(_jdata['properties'][0])

    def tearDown(self):
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual('eos', self.eos.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.eos.task_param())

    def test_make_confs_0(self):

        if not os.path.exists(os.path.join(self.equi_path, 'CONTCAR')):
            with self.assertRaises(RuntimeError):
                self.eos.make_confs(self.target_path, self.equi_path)
        shutil.copy(os.path.join(self.source_path, 'CONTCAR'), os.path.join(self.equi_path, 'CONTCAR'))
        task_list = self.eos.make_confs(self.target_path, self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, 'task.*'))

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar0['ISIF'] = 4

        for ii in dfm_dirs:
            self.assertTrue(os.path.isfile(os.path.join(ii, 'POSCAR')))
            eos_json_file = os.path.join(ii, 'eos.json')
            self.assertTrue(os.path.isfile(eos_json_file))
            eos_json = loadfn(eos_json_file)
            self.assertEqual(os.path.realpath(os.path.join(ii, 'POSCAR.orig')),
                             os.path.realpath(os.path.join(self.equi_path, 'CONTCAR')))
            sys = dpdata.System(os.path.join(ii, 'POSCAR'))
            natoms = sys.get_natoms()
            self.assertAlmostEqual(eos_json['volume'], np.linalg.det(sys['cells'][0]) / natoms)

    def test_make_confs_1(self):
        self.eos.reprod = True
        shutil.copy(os.path.join(self.source_path, 'CONTCAR'), os.path.join(self.equi_path, 'CONTCAR'))
        with self.assertRaises(RuntimeError):
            self.eos.make_confs(self.target_path, self.equi_path)

import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
import dpdata
from monty.serialization import loadfn, dumpfn
from pymatgen.analysis.elasticity.strain import Strain, Deformation
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.Elastic import Elastic


class TestElastic(unittest.TestCase):

    def setUp(self):
        _jdata = {
            "structures": ["confs/std-fcc"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": ".",
                "potcars": {"Al": "vasp_input/POT_Al"}
            },
            "properties": [
                {
                    "skip":False,
                    "type": "elastic",
                    "norm_deform": 2e-2,
                    "shear_deform": 5e-2
                }
            ]
        }

        self.equi_path = 'confs/std-fcc/relaxation/task_relax'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/std-fcc/elastic_00'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata['properties']

        self.elastic = Elastic(_jdata['properties'][0])

    def tearDown(self):
        if os.path.exists(os.path.join(self.equi_path,'..')):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual('elastic', self.elastic.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.elastic.task_param())

    def test_make_confs(self):

        shutil.copy(os.path.join(self.source_path, 'Al-fcc.json'), os.path.join(self.equi_path, 'result.json'))
        if not os.path.exists(os.path.join(self.equi_path, 'CONTCAR')):
            with self.assertRaises(RuntimeError):
                self.elastic.make_confs(self.target_path, self.equi_path)
        shutil.copy(os.path.join(self.source_path, 'CONTCAR_Al_fcc'), os.path.join(self.equi_path, 'CONTCAR'))
        task_list = self.elastic.make_confs(self.target_path, self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, 'task.*'))

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar0['ISIF'] = 4

        self.assertEqual(os.path.realpath(os.path.join(self.equi_path, 'CONTCAR')),
                         os.path.realpath(os.path.join(self.target_path, 'POSCAR')))
        ref_st = Structure.from_file(os.path.join(self.target_path, 'POSCAR'))
        dfm_dirs.sort()
        for ii in dfm_dirs:
            st_file = os.path.join(ii, 'POSCAR')
            self.assertTrue(os.path.isfile(st_file))
            strain_json_file = os.path.join(ii, 'strain.json')
            self.assertTrue(os.path.isfile(strain_json_file))

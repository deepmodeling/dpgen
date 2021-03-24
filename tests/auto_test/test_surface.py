import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
import dpdata
from monty.serialization import loadfn, dumpfn
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar
from pymatgen.core.surface import SlabGenerator

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.Surface import Surface


class TestSurface(unittest.TestCase):

    def setUp(self):
        _jdata = {
            "structures": ["confs/mp-141"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": ".",
                "potcars": {"Yb": "vasp_input/POTCAR"}
            },
            "properties": [
                {
                    "type": "surface",
                    "min_slab_size": 10,
                    "min_vacuum_size": 11,
                    "pert_xz": 0.01,
                    "max_miller": 1,
                    "cal_type": "relaxation"
                }
            ]
        }

        self.equi_path = 'confs/mp-141/relaxation/relax_task'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/mp-141/surface_00'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata['properties']

        self.surface = Surface(_jdata['properties'][0])

    def tearDown(self):
        if os.path.exists(os.path.abspath(os.path.join(self.equi_path,'..'))):
            shutil.rmtree(os.path.abspath(os.path.join(self.equi_path,'..')))
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual('surface', self.surface.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.surface.task_param())

    def test_make_confs_0(self):
        if not os.path.exists(os.path.join(self.equi_path, 'CONTCAR')):
            with self.assertRaises(RuntimeError):
                self.surface.make_confs(self.target_path, self.equi_path)
        shutil.copy(os.path.join(self.source_path, 'mp-141.vasp'), os.path.join(self.equi_path, 'CONTCAR'))
        task_list = self.surface.make_confs(self.target_path, self.equi_path)
        self.assertEqual(len(task_list), 7)
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
            st0 = Structure.from_file(st_file)
            st1_file = os.path.join(ii, 'POSCAR.tmp')
            self.assertTrue(os.path.isfile(st1_file))
            st1 = Structure.from_file(st1_file)
            miller_json_file = os.path.join(ii, 'miller.json')
            self.assertTrue(os.path.isfile(miller_json_file))
            miller_json = loadfn(miller_json_file)
            sl = SlabGenerator(ref_st,
                               miller_json,
                               self.prop_param[0]["min_slab_size"],
                               self.prop_param[0]["min_vacuum_size"])
            slb = sl.get_slab()
            st2 = Structure(slb.lattice, slb.species, slb.frac_coords)
            self.assertEqual(len(st1), len(st2))

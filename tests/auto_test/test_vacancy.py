import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
import dpdata
from monty.serialization import loadfn, dumpfn
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.core import Vacancy as pmg_Vacancy

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.Vacancy import Vacancy


class TestVacancy(unittest.TestCase):

    def setUp(self):
        _jdata = {
            "structures": ["confs/hp-Li"],
            "interaction": {
                "type": "vasp",
                "incar": "vasp_input/INCAR.rlx",
                "potcar_prefix": ".",
                "potcars": {"Yb": "vasp_input/POTCAR"}
            },
            "properties": [
                {
                    "type": "vacancy",
                    "supercell": [1, 1, 1]
                }
            ]
        }

        self.equi_path = 'confs/hp-Li/relaxation/relax_task'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/hp-Li/vacancy_00'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata['properties']

        self.vacancy = Vacancy(_jdata['properties'][0])

    def tearDown(self):
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual('vacancy', self.vacancy.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.vacancy.task_param())

    def test_make_confs_0(self):
        if not os.path.exists(os.path.join(self.equi_path, 'CONTCAR')):
            with self.assertRaises(RuntimeError):
                self.vacancy.make_confs(self.target_path, self.equi_path)
        shutil.copy(os.path.join(self.source_path, 'CONTCAR'), os.path.join(self.equi_path, 'CONTCAR'))
        task_list = self.vacancy.make_confs(self.target_path, self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, 'task.*'))
        self.assertEqual(len(dfm_dirs), 5)

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar0['ISIF'] = 4

        self.assertEqual(os.path.realpath(os.path.join(self.equi_path, 'CONTCAR')),
                         os.path.realpath(os.path.join(self.target_path, 'POSCAR')))
        ref_st = Structure.from_file(os.path.join(self.target_path, 'POSCAR'))
        sga = SpacegroupAnalyzer(ref_st)
        sym_st = sga.get_symmetrized_structure()
        equiv_site_seq = list(sym_st.equivalent_sites)
        dfm_dirs.sort()
        for ii in dfm_dirs:
            st_file = os.path.join(ii, 'POSCAR')
            self.assertTrue(os.path.isfile(st_file))
            st0 = Structure.from_file(st_file)
            vac_site = equiv_site_seq.pop(0)
            vac = pmg_Vacancy(ref_st, vac_site[0])
            st1 = vac.get_supercell_structure(sc_mat=np.eye(3)*self.prop_param[0]['supercell'])
            self.assertEqual(st0, st1)

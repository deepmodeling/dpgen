import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
import dpdata

from monty.serialization import loadfn, dumpfn
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Incar
from pymatgen.io.ase import AseAtomsAdaptor
from ase.lattice.cubic import BodyCenteredCubic as bcc
from ase.lattice.cubic import FaceCenteredCubic as fcc
from ase.lattice.hexagonal import HexagonalClosedPacked as hcp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.Gamma import Gamma

class TestGamma(unittest.TestCase):

    def setUp(self):
        _jdata = {
            "structures":    ["confs/hp-Mo"],
            "interaction": {
                "type":          "vasp",
                "incar":         "vasp_input/INCAR_Mo",
                "potcar_prefix": "vasp_input",
                "potcars":      {"Mo": "POTCAR_Mo"}
            },
            "properties": [
                {
                    "type": "gamma",
                    "lattice_type": "bcc",
                    "miller_index": [1, 1, 0],
                    "displace_direction": [1, 1, 1],
                    "min_supercell_size": [1, 1, 10],
                    "min_vacuum_size": 10,
                    "add_fix": ["true", "true", "false"],
                    "n_steps": 20
                }
            ]
        }

        self.equi_path = 'confs/hp-Mo/relaxation/relax_task'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/hp-Mo/gamma_00'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)
        if not os.path.exists(self.target_path):
            os.makedirs(self.target_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata['properties']

        self.gamma = Gamma(_jdata['properties'][0])

    def tearDown(self):
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual('gamma', self.gamma.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.gamma.task_param())

    def test_make_confs_bcc(self):
        if not os.path.exists(os.path.join(self.equi_path, 'CONTCAR')):
            with self.assertRaises(RuntimeError):
                self.gamma.make_confs(self.target_path, self.equi_path)
        shutil.copy(os.path.join(self.source_path, 'CONTCAR_Mo_bcc'), os.path.join(self.equi_path, 'CONTCAR'))
        task_list = self.gamma.make_confs(self.target_path, self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, 'task.*'))
        self.assertEqual(len(dfm_dirs), self.gamma.n_steps+1)

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
            with open(st1_file, mode='r') as f:
                z_coord_str = f.readlines()[-1].split()[-2]
                z_coord = float(z_coord_str)
            self.assertTrue(z_coord <= 1)

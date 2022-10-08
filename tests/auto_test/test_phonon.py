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

from dpgen.auto_test.Phonon import Phonon

class TestPhonon(unittest.TestCase):

    def setUp(self):
        _jdata = {
            "structures":       ["confs/std-fcc"],
            "interaction": {
                "type":          "vasp",
                "incar":"vasp_input/INCAR",
                "potcar_prefix" : "vasp_input",
                "potcars":{"Al": "POT_Al"}
            },
            "properties": [
                {
                    "type":         "phonon",
                    "band_path": "0.5000   0.5000   0.5000  0.0000   0.0000   0.0000 0.5000   0.0000   0.5000 0.5000   0.2500   0.7500",
                    "supercell_matrix":[5,5,5],
                    "primitive": True ,
                    "approach":"linear"
                }
            ]
        }

        self.equi_path = 'confs/std-fcc/relaxation/relax_task'
        self.source_path = 'equi/vasp'
        self.target_path = 'confs/std-fcc/phonon_00'
        if not os.path.exists(self.equi_path):
            os.makedirs(self.equi_path)

        self.confs = _jdata["structures"]
        self.inter_param = _jdata["interaction"]
        self.prop_param = _jdata['properties']

        self.phonon = Phonon(_jdata['properties'][0])
    

    def tearDown(self):
        if os.path.exists(self.equi_path):
            shutil.rmtree(self.equi_path)
        if os.path.exists(self.target_path):
            shutil.rmtree(self.target_path)

    def test_task_type(self):
        self.assertEqual('phonon',self.phonon.task_type())

    def test_task_param(self):
        self.assertEqual(self.prop_param[0], self.phonon.task_param())

    def test_make_confs(self):
        shutil.copy(os.path.join(self.source_path, 'CONTCAR_Al_fcc'), os.path.join(self.equi_path, 'CONTCAR'))
        task_list = self.phonon.make_confs(self.target_path, self.equi_path)
        dfm_dirs = glob.glob(os.path.join(self.target_path, 'task.*'))

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar0['ISIF'] = 4

        for ii in dfm_dirs:
            self.assertTrue(os.path.isfile(os.path.join(ii, 'POSCAR')))
            band_conf_file = os.path.join(ii, 'band.conf')
            self.assertTrue(os.path.isfile(band_conf_file))
            self.assertEqual(os.path.realpath(os.path.join(self.target_path, 'POSCAR.orig')),
                             os.path.realpath(os.path.join(self.equi_path, 'CONTCAR')))

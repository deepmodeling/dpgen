import os, sys, json, glob, shutil
import dpdata
import numpy as np
from monty.serialization import loadfn, dumpfn
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from pymatgen.io.vasp import Incar
from dpgen.auto_test.common_equi import make_equi, post_equi
from dpgen.auto_test.calculator import make_calculator


class TestEqui(unittest.TestCase):
    jdata = {
        "structures": ["confs/std-fcc"],
        "interaction": {
            "type": "vasp",
            "incar": "vasp_input/INCAR.rlx",
            "potcar_prefix": ".",
            "potcars": {"Al": "vasp_input/POT_Al"}
        },
        "relaxation": {
            "cal_type": "relaxation",
            "cal_setting":{"input_prop": "vasp_input/INCAR"}
        }
    }

    def tearDown(self):
        if os.path.isfile('confs/std-fcc/POSCAR'):
           os.remove('confs/std-fcc/POSCAR')
        if os.path.exists('confs/std-fcc/relaxation'):
            shutil.rmtree('confs/std-fcc/relaxation')

    def test_make_equi(self):
        confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        relax_param = self.jdata["relaxation"]
        make_equi(confs, inter_param, relax_param)

        self.assertTrue(os.path.isfile("confs/std-fcc/POSCAR"))

        target_path = 'confs/std-fcc/relaxation/relax_task'
        source_path = 'vasp_input'

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR'))
        incar1 = Incar.from_file(os.path.join(target_path, 'INCAR'))
        self.assertTrue(incar0 == incar1)

        with open(os.path.join('vasp_input', 'POT_Al')) as fp:
            pot0 = fp.read()
        with open(os.path.join(target_path, 'POTCAR')) as fp:
            pot1 = fp.read()
        self.assertEqual(pot0, pot1)

        self.assertTrue(os.path.isfile(os.path.join(target_path, 'KPOINTS')))

        task_json_file = os.path.join(target_path, 'task.json')
        self.assertTrue(os.path.isfile(task_json_file))
        task_json = loadfn(task_json_file)
        self.assertEqual(task_json, relax_param)

        inter_json_file = os.path.join(target_path, 'inter.json')
        self.assertTrue(os.path.isfile(inter_json_file))
        inter_json = loadfn(inter_json_file)
        self.assertEqual(inter_json, inter_param)

        self.assertTrue(os.path.islink(os.path.join(target_path, 'POSCAR')))


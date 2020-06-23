import os,sys,json,glob,shutil
import dpdata
import numpy as np
import unittest
from monty.serialization import loadfn,dumpfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import make_kspacing_kpoints
from .context import setUpModule

from pymatgen.io.vasp import Incar
from dpgen.auto_test.common_prop import make_property

class TestMakeProperty(unittest.TestCase):
    jdata={
    "structures":    ["confs/hp-*"],
    "interaction": {
        "type":      "vasp",
        "incar":     "vasp_input/INCAR.rlx",
        "potcar_prefix":".",
        "type_map": ['Li'],
        "potcars":  {"Li": "vasp_input/POTCAR"}
    },
    "properties": [
        {
         "type":         "eos",
         "vol_start":    10,
         "vol_end":      30,
         "vol_step":     3,
         "change_box":   True
        }
    ]
    }

    def tearDown(self):
        if os.path.exists('confs/hp-Li/eos_00'):
            shutil.rmtree('confs/hp-Li/eos_00')
        if os.path.exists('confs/hp-Li/relaxation'):
            shutil.rmtree('confs/hp-Li/relaxation')

    def test_make_eos (self):
        confs=self.jdata["structures"]
        inter_param=self.jdata["interaction"]
        property_list=self.jdata["properties"]
        
        target_path = 'confs/hp-Li/eos_00'
        equi_path =   'confs/hp-Li/relaxation'
        source_path = 'equi/vasp'

        if not os.path.exists(equi_path):
           os.mkdir(equi_path)
        shutil.copy(os.path.join(source_path,'CONTCAR'),os.path.join(equi_path,'CONTCAR'))

        make_property(confs,inter_param,property_list)

        dfm_dirs = glob.glob(os.path.join(target_path, 'task.*'))

        incar0 = Incar.from_file(os.path.join('vasp_input', 'INCAR.rlx'))
        incar0['ISIF'] = 4

        with open(os.path.join('vasp_input', 'POTCAR')) as fp:
             pot0 = fp.read()
        for ii in dfm_dirs:
            self.assertTrue(os.path.isfile(os.path.join(ii, 'KPOINTS')))
            incar1 = Incar.from_file(os.path.join(ii, 'INCAR'))
            self.assertTrue(incar0 == incar1)
            self.assertTrue(os.path.isfile(os.path.join(ii, 'INCAR')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'POSCAR')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'POTCAR')))
            self.assertTrue(os.path.isfile(os.path.join(ii, 'task.json')))
            inter_json_file=os.path.join(ii, 'inter.json')
            self.assertTrue(os.path.isfile(inter_json_file))
            inter_json=loadfn(inter_json_file)
            self.assertEqual(inter_json,inter_param)
            self.assertEqual(os.path.realpath(os.path.join(ii, 'POSCAR.orig')),
                             os.path.realpath(os.path.join(equi_path, 'CONTCAR')))
            with open(os.path.join(ii, 'POTCAR')) as fp:
                poti = fp.read()
            self.assertEqual(pot0, poti)
     #       sys = dpdata.System(os.path.join(ii, 'POSCAR'))
     #       vol = float(ii.split('/')[-1].split('-')[1])
     #       natoms = sys.get_natoms()
     #       self.assertAlmostEqual(vol, np.linalg.det(sys['cells'][0]) / natoms)

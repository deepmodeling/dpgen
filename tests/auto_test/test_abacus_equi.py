import os, sys, json, glob, shutil
from monty.serialization import loadfn
import unittest
from dpgen.generator.lib import abacus_scf
from dpgen.auto_test.common_equi import make_equi, post_equi

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'
from .context import setUpModule

class TestEqui(unittest.TestCase):
    jdata = {
        "structures": ["confs/fcc-Al"],
        "interaction": {
            "type": "abacus",
            "incar": "abacus_input/INPUT",
            "potcar_prefix": "abacus_input",
            "potcars": {"Al": "Al_ONCV_PBE-1.0.upf"},
            "orb_files": {"Al":"Al_gga_9au_100Ry_4s4p1d.orb"}
        },
        "relaxation": {
            "cal_type": "relaxation",
            "cal_setting":{"input_prop": "abacus_input/INPUT"}
        }
    }

    def tearDown(self):
        if os.path.exists('confs/fcc-Al/relaxation'):
            shutil.rmtree('confs/fcc-Al/relaxation')

    def test_make_equi(self):
        confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        relax_param = self.jdata["relaxation"]
        make_equi(confs, inter_param, relax_param)

        target_path = 'confs/fcc-Al/relaxation/relax_task'

        incar0 = abacus_scf.get_abacus_input_parameters(os.path.join('abacus_input', 'INPUT')) 
        incar1 = abacus_scf.get_abacus_input_parameters(os.path.join(target_path, 'INPUT'))
        self.assertTrue(incar0 == incar1)

        with open(os.path.join('abacus_input', 'Al_ONCV_PBE-1.0.upf')) as fp:
            pot0 = fp.read()
        with open(os.path.join(target_path,'pp_orb', 'Al_ONCV_PBE-1.0.upf')) as fp:
            pot1 = fp.read()
        self.assertEqual(pot0, pot1)

        with open(os.path.join('abacus_input', 'Al_gga_9au_100Ry_4s4p1d.orb')) as fp:
            pot0 = fp.read()
        with open(os.path.join(target_path,'pp_orb', 'Al_gga_9au_100Ry_4s4p1d.orb')) as fp:
            pot1 = fp.read()
        self.assertEqual(pot0, pot1)

        self.assertTrue(os.path.isfile(os.path.join(target_path, 'KPT')))

        task_json_file = os.path.join(target_path, 'task.json')
        self.assertTrue(os.path.isfile(task_json_file))
        task_json = loadfn(task_json_file)
        self.assertEqual(task_json, relax_param)

        inter_json_file = os.path.join(target_path, 'inter.json')
        self.assertTrue(os.path.isfile(inter_json_file))
        inter_json = loadfn(inter_json_file)
        self.assertEqual(inter_json, inter_param)

        self.assertTrue(os.path.islink(os.path.join(target_path, 'STRU')))

    def test_post_equi(self):
       confs = self.jdata["structures"]
       inter_param = self.jdata["interaction"]
       relax_param = self.jdata["relaxation"]
       target_path = 'confs/fcc-Al/relaxation/relax_task'
       source_path = 'equi/abacus'

       make_equi(confs, inter_param, relax_param)
       os.mkdir(os.path.join(target_path, 'OUT.ABACUS'))
       shutil.copy(os.path.join(source_path, 'INPUT'), os.path.join(target_path, 'INPUT'))
       shutil.copy(os.path.join(source_path, 'STRU'), os.path.join(target_path, 'STRU'))
       shutil.copy(os.path.join(source_path, 'running_cell-relax.log'), os.path.join(target_path, 'OUT.ABACUS','running_cell-relax.log'))
       post_equi(confs, inter_param)

       result_json_file = os.path.join(target_path, 'result.json')
       self.assertTrue(os.path.isfile(result_json_file))

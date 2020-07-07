import os, sys, json, glob, shutil
import dpdata
import numpy as np
import unittest
from monty.serialization import loadfn, dumpfn
from pymatgen.io.vasp import Incar

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'auto_test'

from .context import make_kspacing_kpoints
from .context import setUpModule

from dpgen.auto_test.common_prop import make_property
from dpgen.auto_test.refine import make_refine


class TestMakeProperty(unittest.TestCase):
    jdata = {
        "structures": ["confs/std-fcc"],
        "interaction": {
            "type": "vasp",
            "incar": "vasp_input/INCAR.rlx",
            "potcar_prefix": "vasp_input",
            "potcars": {"Al": "POT_Al"}
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
                                                  "type": "vasp",
                                                  "incar": "vasp_input/INCAR.rlx",
                                                  "potcar_prefix":"vasp_input",
                                                  "potcars": {"Al": "POT_Al"}
                                                  }
                                     }
                     }
        ]
    }

    def tearDown(self):
        if os.path.exists('confs/std-fcc/eos_00'):
            shutil.rmtree('confs/std-fcc/eos_00')
        if os.path.exists('confs/std-fcc/eos_02'):
            shutil.rmtree('confs/std-fcc/eos_02')
        if os.path.exists('confs/std-fcc/relaxation'):
            shutil.rmtree('confs/std-fcc/relaxation')

    def test_make_eos(self):

        pwd=os.getcwd()
        confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        property_list = self.jdata["properties"]

        target_path_0 = 'confs/std-fcc/eos_00'
        target_path_2 = 'confs/std-fcc/eos_02'
        equi_path = 'confs/std-fcc/relaxation/relax_task'
        source_path = 'equi/vasp'

        if not os.path.exists(equi_path):
            os.makedirs(equi_path)
        shutil.copy(os.path.join(source_path, 'CONTCAR_Al_fcc'), os.path.join(equi_path, 'CONTCAR'))

        make_property(confs, inter_param, property_list)

        dfm_dirs_0 = glob.glob(os.path.join(target_path_0, 'task.*'))

        for ii in dfm_dirs_0:
            self.assertTrue(os.path.isfile(os.path.join(ii, 'POSCAR')))
            shutil.copy(os.path.join(ii, 'POSCAR'),os.path.join(ii, 'CONTCAR'))

        path_to_work = os.path.abspath(target_path_0)
        new_prop_list=[
        {
         "type":         "eos",
         "init_from_suffix": "00",
         "output_suffix": "02",
         "cal_setting": {
                           "relax_pos": True,
                            "relax_shape": True,
                            "relax_vol": False,
         "input_prop":  "lammps_input/lammps_high"}
        }
        ]
        #ret=make_refine('00', '02', path_to_work)
        #self.assertEqual(len(ret),len(dfm_dirs_0))
        make_property(confs, inter_param, new_prop_list)
        self.assertTrue(os.path.isdir(path_to_work.replace('00','02')))
        os.chdir(pwd)
        dfm_dirs_2 = glob.glob(os.path.join(target_path_2, 'task.*'))
        self.assertEqual(len(dfm_dirs_2),len(dfm_dirs_0))

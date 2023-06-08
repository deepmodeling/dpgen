import os
import shutil
import sys
import unittest

from monty.serialization import loadfn

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "auto_test"

from pymatgen.io.vasp import Incar

from dpgen.auto_test.common_equi import make_equi, post_equi

from .context import setUpModule  # noqa: F401


class TestEqui(unittest.TestCase):
    jdata = {
        "structures": ["confs/hp-Li"],
        "interaction": {
            "type": "vasp",
            "incar": "vasp_input/INCAR.rlx",
            "potcar_prefix": ".",
            "potcars": {"Li": "vasp_input/POTCAR"},
        },
        "relaxation": {
            "cal_type": "relaxation",
            "cal_setting": {"input_prop": "vasp_input/INCAR"},
        },
    }

    def tearDown(self):
        if os.path.exists("confs/hp-Li/relaxation"):
            shutil.rmtree("confs/hp-Li/relaxation")

    def test_make_equi(self):
        confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        relax_param = self.jdata["relaxation"]
        make_equi(confs, inter_param, relax_param)

        target_path = "confs/hp-Li/relaxation/relax_task"
        source_path = "vasp_input"

        incar0 = Incar.from_file(os.path.join("vasp_input", "INCAR"))
        incar1 = Incar.from_file(os.path.join(target_path, "INCAR"))
        self.assertTrue(incar0 == incar1)

        with open(os.path.join("vasp_input", "POTCAR")) as fp:
            pot0 = fp.read()
        with open(os.path.join(target_path, "POTCAR")) as fp:
            pot1 = fp.read()
        self.assertEqual(pot0, pot1)

        self.assertTrue(os.path.isfile(os.path.join(target_path, "KPOINTS")))

        task_json_file = os.path.join(target_path, "task.json")
        self.assertTrue(os.path.isfile(task_json_file))
        task_json = loadfn(task_json_file)
        self.assertEqual(task_json, relax_param)

        inter_json_file = os.path.join(target_path, "inter.json")
        self.assertTrue(os.path.isfile(inter_json_file))
        inter_json = loadfn(inter_json_file)
        self.assertEqual(inter_json, inter_param)

        self.assertTrue(os.path.islink(os.path.join(target_path, "POSCAR")))

    def test_post_equi(self):
        confs = self.jdata["structures"]
        inter_param = self.jdata["interaction"]
        relax_param = self.jdata["relaxation"]
        target_path = "confs/hp-Li/relaxation/relax_task"
        source_path = "equi/vasp"

        poscar = os.path.join(source_path, "POSCAR")
        make_equi(confs, inter_param, relax_param)
        shutil.copy(
            os.path.join(source_path, "OUTCAR"), os.path.join(target_path, "OUTCAR")
        )
        shutil.copy(
            os.path.join(source_path, "CONTCAR"), os.path.join(target_path, "CONTCAR")
        )
        post_equi(confs, inter_param)

        result_json_file = os.path.join(target_path, "result.json")
        result_json = loadfn(result_json_file)
        self.assertTrue(os.path.isfile(result_json_file))

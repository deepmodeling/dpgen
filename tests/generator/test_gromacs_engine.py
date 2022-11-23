import os, sys, glob, shutil
import unittest
import json
import numpy as np
import importlib

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
dirname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "gromacs")

from .context import make_model_devi
from .context import make_fp_gaussian

def _make_fake_graphs(train_path):
    if not os.path.exists(train_path):
        os.mkdir(train_path)
    for ii in range(4):
        with open(os.path.join(train_path, f"graph.{ii:03}.pb"), 'w+') as f:
            f.write("Fake Model")

class TestGromacsModelDeviEngine(unittest.TestCase):
    def setUp(self):
        self.dirname = dirname
        self.jdata = {
            "type_map": ["H", "C", "N", "O", "Cl"],
            "mass_map": [2, 12, 14, 16, 35],
            "sys_configs_prefix": self.dirname,
            "sys_configs": [["model_devi_case"]],
            "sys_format": "gromacs/gro",
            "model_devi_engine": "gromacs",
            "gromacs_settings": {
                "mdp_filename":      "md.mdp",
                "topol_filename":    "processed.top",
                "conf_filename":     "npt.gro",
                "index_filename":    "index.raw",
                "ref_filename":      "em.tpr",
                "model_devi_script": "model_devi.py",
                "traj_filename":     "deepmd_traj.gro"
            },
            "model_devi_dt":         0.001,
            "model_devi_f_trust_lo": 0.05,
            "model_devi_f_trust_hi": 0.10,
            "model_devi_clean_traj": False,
            "model_devi_skip":       0,
            "model_devi_nopbc":      True,
            "model_devi_jobs": [
                {
                    "ensemble": "nvt",
                    "nsteps":   5000,
                    "sys_idx":  [0],
                    "trj_freq": 10
                }
            ],
            "shuffle_poscar": False,
            "fp_style": "gaussian",
            "shuffle_poscar": False,
            "fp_task_max": 20,
            "fp_task_min": 1,
            "fp_pp_path": "./",
            "fp_pp_files": [],
            "fp_params": {
                "keywords": "force m062x/6-31g(d) nosymm",
                "nproc": 16,
                "multiplicity": "auto"
            }
        }
        self.iter_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "iter.000000")
        if not os.path.exists(self.iter_path):
            os.mkdir(self.iter_path)
        self.train_path = os.path.join(self.iter_path, "00.train")
        self.model_devi_path = os.path.join(self.iter_path, "01.model_devi")
        self.model_devi_task_path = os.path.join(self.model_devi_path, "task.000.000000")
        self.fp_path = os.path.join(self.iter_path, "02.fp")
        _make_fake_graphs(self.train_path)
    
    def _check_dir(self, wdir, post=True):
        for key in self.jdata['gromacs_settings'].keys():
            if key != "traj_filename":
                self.assertTrue(os.path.exists(os.path.join(wdir, self.jdata['gromacs_settings'][key])))
            else:
                if post:
                    self.assertTrue(os.path.exists(os.path.join(wdir, self.jdata['gromacs_settings'][key])))
    
    def _copy_outputs(self, path_1, path_2):
        shutil.copy(os.path.join(path_1, "deepmd_traj.gro"), os.path.join(path_2, "deepmd_traj.gro"))
        shutil.copy(os.path.join(path_1, "model_devi.out"), os.path.join(path_2, "model_devi.out"))
        shutil.copytree(os.path.join(path_1, "traj"), os.path.join(path_2, "traj"))

    
    @unittest.skipIf(importlib.util.find_spec("openbabel") !=  None, "when openbabel is found, this test will be skipped. ")
    def test_make_model_devi_gromacs_without_openbabel(self):
        flag = make_model_devi(iter_index=0,
                               jdata=self.jdata,
                               mdata={"deepmd_version": "2.0"})
        self.assertTrue(flag)
        self.assertTrue(os.path.exists(self.model_devi_path))
        self.assertTrue(os.path.exists(self.model_devi_task_path))
        self._check_dir(self.model_devi_task_path, post=False)
        self._copy_outputs(os.path.join(self.dirname, "outputs"), self.model_devi_task_path)
        self._check_dir(self.model_devi_task_path, post=True)
    
    @unittest.skipIf(importlib.util.find_spec("openbabel") is None, "requires openbabel")
    def test_make_model_devi_gromacs_with_openbabel(self):
        flag = make_model_devi(iter_index=0,
                               jdata=self.jdata,
                               mdata={"deepmd_version": "2.0"})
        self._copy_outputs(os.path.join(self.dirname, "outputs"), self.model_devi_task_path)
        make_fp_gaussian(iter_index=0, jdata=self.jdata)
        candi = np.loadtxt(os.path.join(self.fp_path, "candidate.shuffled.000.out"), dtype=np.str)
        self.assertEqual(sorted([int(i) for i in candi[:,1]]), [0,10,20,30,50])
        
     
    def tearDown(self):
        # pass
        shutil.rmtree(self.iter_path)
if __name__ == '__main__':
    unittest.main()


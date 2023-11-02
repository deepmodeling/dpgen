import os
import shutil
import sys
import unittest
from pathlib import Path

import dpdata
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "simplify"
from .context import dpgen


class TestSimplifyModelDevi(unittest.TestCase):
    def setUp(self):
        self.work_path = Path("iter.000001") / dpgen.simplify.simplify.model_devi_name
        self.work_path.mkdir(exist_ok=True, parents=True)
        self.system = dpdata.System(
            data={
                "atom_names": ["H"],
                "atom_numbs": [1],
                "atom_types": np.zeros((1,), dtype=int),
                "coords": np.zeros((1, 1, 3), dtype=np.float32),
                "cells": np.zeros((1, 3, 3), dtype=np.float32),
                "orig": np.zeros(3, dtype=np.float32),
                "nopbc": True,
                "energies": np.zeros((1,), dtype=np.float32),
                "forces": np.zeros((1, 1, 3), dtype=np.float32),
            }
        )
        self.system.to_deepmd_npy(
            self.work_path / "data.rest.old" / self.system.formula
        )
        model_devi = np.array([[0, 0.2, 0.1, 0.15, 0.2, 0.1, 0.15, 0.2]])
        np.savetxt(
            self.work_path / "details",
            model_devi,
            fmt=["%12d"] + ["%19.6e" for _ in range(7)],
            header="data.rest.old/"
            + self.system.formula
            + "\n step max_devi_v min_devi_v avg_devi_v max_devi_f min_devi_f avg_devi_f devi_e",
        )
        np.savetxt(
            self.work_path / "true_error",
            model_devi,
            fmt=["%12d"] + ["%19.6e" for _ in range(7)],
            header="data.rest.old/"
            + self.system.formula
            + "\n step max_devi_v min_devi_v avg_devi_v max_devi_f min_devi_f avg_devi_f devi_e",
        )

    def tearDown(self):
        shutil.rmtree("iter.000001", ignore_errors=True)

    def test_post_model_devi_f_candidate(self):
        dpgen.simplify.simplify.post_model_devi(
            1,
            {
                "model_devi_f_trust_lo": 0.15,
                "model_devi_f_trust_hi": 0.25,
                "model_devi_e_trust_lo": float("inf"),
                "model_devi_e_trust_hi": float("inf"),
                "iter_pick_number": 1,
            },
            {},
        )
        assert (self.work_path / "data.picked" / self.system.formula).exists()

    def test_post_model_devi_e_candidate(self):
        dpgen.simplify.simplify.post_model_devi(
            1,
            {
                "model_devi_e_trust_lo": 0.15,
                "model_devi_e_trust_hi": 0.25,
                "model_devi_f_trust_lo": float("inf"),
                "model_devi_f_trust_hi": float("inf"),
                "iter_pick_number": 1,
            },
            {},
        )
        assert (self.work_path / "data.picked" / self.system.formula).exists()

    def test_post_model_devi_f_failed(self):
        with self.assertRaises(RuntimeError):
            dpgen.simplify.simplify.post_model_devi(
                1,
                {
                    "model_devi_f_trust_lo": 0.0,
                    "model_devi_f_trust_hi": 0.0,
                    "model_devi_e_trust_lo": float("inf"),
                    "model_devi_e_trust_hi": float("inf"),
                    "iter_pick_number": 1,
                },
                {},
            )

    def test_post_model_devi_e_failed(self):
        with self.assertRaises(RuntimeError):
            dpgen.simplify.simplify.post_model_devi(
                1,
                {
                    "model_devi_e_trust_lo": 0.0,
                    "model_devi_e_trust_hi": 0.0,
                    "model_devi_f_trust_lo": float("inf"),
                    "model_devi_f_trust_hi": float("inf"),
                    "iter_pick_number": 1,
                },
                {},
            )

    def test_post_model_devi_accurate(self):
        dpgen.simplify.simplify.post_model_devi(
            1,
            {
                "model_devi_e_trust_lo": 0.3,
                "model_devi_e_trust_hi": 0.4,
                "model_devi_f_trust_lo": 0.3,
                "model_devi_f_trust_hi": 0.4,
                "iter_pick_number": 1,
            },
            {},
        )
        assert (self.work_path / "data.accurate" / self.system.formula).exists()

    def test_post_model_devi_true_error_candidate(self):
        dpgen.simplify.simplify.post_model_devi(
            1,
            {
                "model_devi_e_trust_lo": 0.15,
                "model_devi_e_trust_hi": 0.25,
                "model_devi_f_trust_lo": float("inf"),
                "model_devi_f_trust_hi": float("inf"),
                "true_error_e_trust_lo": float("inf"),
                "true_error_e_trust_hi": float("inf"),
                "true_error_f_trust_lo": 0.15,
                "true_error_f_trust_hi": 0.25,
                "iter_pick_number": 1,
            },
            {},
        )
        assert (self.work_path / "data.picked" / self.system.formula).exists()

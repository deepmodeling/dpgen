"""Validate CALYPSO model-deviation parameter documentation."""

import unittest

from dargs import Argument

from dpgen.generator.arginfo import model_devi_args


class TestCalypsoArginfo(unittest.TestCase):
    def setUp(self):
        self.arginfo = Argument("model_devi", dict, sub_variants=model_devi_args())
        self.selection = {
            "model_devi_skip": 0,
            "model_devi_f_trust_lo": 0.05,
            "model_devi_f_trust_hi": 0.15,
            "model_devi_clean_traj": True,
        }

    def test_native_mode(self):
        data = {
            "model_devi_engine": "calypso",
            **self.selection,
            "model_devi_dt": 0.002,
            "shuffle_poscar": False,
            "model_devi_jobs": [
                {
                    "times": [0, 1],
                    "NameOfAtoms": ["Mg", "Al"],
                    "NumberOfAtoms": [1, 1],
                    "NumberOfFormula": [1, 2],
                    "Volume": 30.0,
                    "DistanceOfIon": [[1.4, 1.5], [1.5, 1.4]],
                    "PsoRatio": 0.6,
                    "PopSize": 30,
                    "MaxStep": 5,
                    "ICode": 1,
                    "Split": "T",
                    "VSC": "T",
                    "MaxNumAtom": 20,
                    "CtrlRange": [[1, 10], [1, 10]],
                    "PSTRESS": [0.0, 100.0],
                    "fmax": 0.01,
                    "task_min": 1,
                }
            ],
        }

        normalized = self.arginfo.normalize_value(data)
        self.arginfo.check_value(normalized, strict=True)

    def test_external_input_mode(self):
        data = {
            "model_devi_engine": "calypso",
            **self.selection,
            "model_devi_jobs": [],
            "calypso_input_path": "calypso_input",
            "model_devi_max_iter": 20,
            "vsc": True,
        }

        normalized = self.arginfo.normalize_value(data)
        self.arginfo.check_value(normalized, strict=True)


if __name__ == "__main__":
    unittest.main()

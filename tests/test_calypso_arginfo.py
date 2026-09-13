"""Validate CALYPSO model-deviation parameter documentation."""

import json
import unittest
from pathlib import Path

from dargs import Argument

from dpgen.generator.arginfo import model_devi_args


class TestCalypsoArginfo(unittest.TestCase):
    def test_pressure_forms_survive_strict_validation(self):
        """The schema accepts the pressure forms supported by input generation."""
        path = Path(__file__).parent.parent / "examples/run/dp-calypso-vasp/param.json"
        job = json.loads(path.read_text())["model_devi_jobs"][0]
        for pressures in (0.0, [0.0], [0.0, 2.5]):
            with self.subTest(pressures=pressures):
                data = {
                    "model_devi_engine": "calypso",
                    **self.selection,
                    "model_devi_jobs": [{**job, "PSTRESS": pressures}],
                }
                normalized = self.arginfo.normalize_value(data)
                self.arginfo.check_value(normalized, strict=True)
                self.assertEqual(normalized["model_devi_jobs"][0]["PSTRESS"], pressures)

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

    def test_legacy_singleton_scalars(self):
        """Checked-in CALYPSO inputs retain their historical list spelling."""
        data = {
            "model_devi_engine": "calypso",
            **self.selection,
            "model_devi_jobs": [
                {
                    "times": [0],
                    "NameOfAtoms": ["Mg"],
                    "NumberOfAtoms": [1],
                    "Volume": [30],
                    "DistanceOfIon": [[1.4]],
                    "PsoRatio": [0.6],
                    "PopSize": [30],
                    "MaxStep": [5],
                    "ICode": [1],
                    "VSC": "T",
                    "MaxNumAtom": [20],
                    "CtrlRange": [[1, 20]],
                    "fmax": [0.01],
                }
            ],
        }

        normalized = self.arginfo.normalize_value(data)
        self.arginfo.check_value(normalized, strict=True)

    def test_selection_arguments_are_accepted(self):
        """CALYPSO exposes the controls consumed by downstream FP selection."""
        arginfo = Argument("model_devi", dict, sub_variants=model_devi_args())
        data = {
            "model_devi_engine": "calypso",
            "model_devi_jobs": [],
            "calypso_input_path": "calypso-input",
            "model_devi_skip": 0,
            "model_devi_f_trust_lo": 0.05,
            "model_devi_f_trust_hi": 0.15,
            "model_devi_clean_traj": True,
            "model_devi_numb_candi_f": 10,
            "model_devi_numb_candi_v": 10,
            "shuffle_poscar": False,
        }

        normalized = arginfo.normalize_value(data)
        arginfo.check_value(normalized, strict=True)

    def test_job_task_min_survives_strict_normalization(self):
        """FP selection consumes task_min from the selected cur_job.json."""
        path = Path(__file__).parent.parent / "examples/run/dp-calypso-vasp/param.json"
        example = json.loads(path.read_text())
        job = example["model_devi_jobs"][0]
        job["task_min"] = 3
        data = {
            "model_devi_engine": "calypso",
            "model_devi_jobs": [job],
            "model_devi_skip": 0,
            "model_devi_f_trust_lo": 0.05,
            "model_devi_f_trust_hi": 0.15,
        }
        arginfo = Argument("model_devi", dict, sub_variants=model_devi_args())
        normalized = arginfo.normalize_value(data)
        arginfo.check_value(normalized, strict=True)
        self.assertEqual(normalized["model_devi_jobs"][0]["task_min"], 3)

    def test_checked_in_example_is_schema_compatible(self):
        """The maintained example's CALYPSO section passes strict validation."""
        param_file = (
            Path(__file__).parent.parent
            / "examples"
            / "run"
            / "dp-calypso-vasp"
            / "param.json"
        )
        with open(param_file) as fp:
            example = json.load(fp)

        model_devi_keys = {
            "model_devi_engine",
            "model_devi_jobs",
            "model_devi_dt",
            "model_devi_skip",
            "model_devi_f_trust_lo",
            "model_devi_f_trust_hi",
            "model_devi_clean_traj",
            "shuffle_poscar",
            "vsc",
        }
        data = {key: example[key] for key in model_devi_keys}
        arginfo = Argument("model_devi", dict, sub_variants=model_devi_args())
        normalized = arginfo.normalize_value(data)
        arginfo.check_value(normalized, strict=True)


if __name__ == "__main__":
    unittest.main()

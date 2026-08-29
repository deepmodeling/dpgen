import json
import unittest
from pathlib import Path

from dargs import Argument

from dpgen.generator.arginfo import model_devi_args
from dpgen.generator.lib.make_calypso import _unwrap_calypso_scalar


class TestCalypsoArginfo(unittest.TestCase):
    def test_legacy_singleton_scalars(self):
        """The checked-in CALYPSO list spelling remains schema-compatible."""
        arginfo = Argument("model_devi", dict, sub_variants=model_devi_args())
        data = {
            "model_devi_engine": "calypso",
            "model_devi_skip": 0,
            "model_devi_f_trust_lo": 0.05,
            "model_devi_f_trust_hi": 0.15,
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

        normalized = arginfo.normalize_value(data)
        arginfo.check_value(normalized, strict=True)

    def test_runtime_unwraps_singleton_scalars(self):
        self.assertEqual(_unwrap_calypso_scalar([30], "Volume"), 30)
        self.assertEqual(_unwrap_calypso_scalar(0.6, "PsoRatio"), 0.6)
        with self.assertRaisesRegex(ValueError, "one-item list"):
            _unwrap_calypso_scalar([1, 2], "PopSize")

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

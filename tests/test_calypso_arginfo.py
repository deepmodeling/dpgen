import json
import unittest
from pathlib import Path

from dargs import Argument

from dpgen.generator.arginfo import model_devi_args


class TestCalypsoArginfo(unittest.TestCase):
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

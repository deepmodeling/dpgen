import unittest

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


if __name__ == "__main__":
    unittest.main()

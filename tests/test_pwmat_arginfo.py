"""Validate PWmat first-principles parameter documentation."""

import unittest

from dargs import Argument

from dpgen.generator.arginfo import fp_style_variant_type_args


class TestPWmatArginfo(unittest.TestCase):
    def setUp(self):
        self.arginfo = Argument("fp", dict, sub_variants=[fp_style_variant_type_args()])
        self.pseudopotentials = {
            "fp_style": "pwmat",
            "fp_pp_path": ".",
            "fp_pp_files": ["C.UPF", "H.UPF"],
        }

    def check(self, data):
        normalized = self.arginfo.normalize_value(data)
        self.arginfo.check_value(normalized, strict=True)

    def test_existing_input_file(self):
        self.check({**self.pseudopotentials, "fp_incar": "etot.input"})

    def test_generated_fp_params(self):
        self.check(
            {
                **self.pseudopotentials,
                "fp_params": {
                    "node1": 4,
                    "node2": 1,
                    "in.atom": "atom.config",
                    "ecut": 50,
                    "e_error": 1e-4,
                    "rho_error": 1e-4,
                    "kspacing": 0.1,
                    "flag_symm": "NONE",
                    "icmix": 1.0,
                    "smearing": 2,
                    "sigma": 0.025,
                    "user_pwmat_params": {"job": "SCF", "out.wg": False},
                },
            }
        )

    def test_compatibility_user_fp_params(self):
        self.check(
            {
                **self.pseudopotentials,
                "user_fp_params": {
                    "node1": 4,
                    "node2": 1,
                    "job": "SCF",
                    "in.atom": "atom.config",
                    "in.psp1": "C.UPF",
                    "in.psp2": "H.UPF",
                    "ecut": 50,
                    "flag_symm": 2,
                    "e_error": 1e-4,
                    "rho_error": 1e-4,
                    "scf_iter0_1": "6 4 3 0.0000 0.025 2",
                    "xcfunctional": "PBE",
                    "kspacing": 0.1,
                },
            }
        )


if __name__ == "__main__":
    unittest.main()

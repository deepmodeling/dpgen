import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dpgen.tools import relabel


class TestRelabelInputs(unittest.TestCase):
    def test_create_tasks_copies_vasp_incar_to_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            job = root / "job"
            output = root / "output"
            pp_dir = job / "pp"
            job.mkdir()
            pp_dir.mkdir()
            (pp_dir / "H.POTCAR").write_text("potential")
            (job / "INCAR.source").write_text("ENCUT = 500\n")
            (job / "param.json").write_text(
                json.dumps(
                    {
                        "mass_map": [1.0],
                        "type_map": ["H"],
                        "sys_configs": [],
                    }
                )
            )
            fp_json = root / "fp.json"
            fp_json.write_text(
                json.dumps(
                    {
                        "fp_style": "vasp",
                        "fp_pp_path": "pp",
                        "fp_pp_files": ["H.POTCAR"],
                        "fp_incar": "INCAR.source",
                    }
                )
            )

            relabel.create_tasks(job, "param.json", output, fp_json, verbose=False)

            self.assertEqual((output / "INCAR").read_text(), "ENCUT = 500\n")
            self.assertEqual((output / "H.POTCAR").read_text(), "potential")

    @patch("dpgen.tools.relabel.make_pwscf")
    def test_pwscf_arguments_match_helper_signature(self, make_pwscf):
        params = {"ecut": 100}
        relabel.make_non_vasp_input(
            "task",
            "pwscf",
            {"user_fp_params": params},
            [1.0],
            "/pp",
            ["H.UPF"],
        )

        make_pwscf.assert_called_once_with(
            "task", params, [1.0], "/pp", ["H.UPF"], True
        )

    def test_pwscf_writes_inputs_for_both_parameter_modes(self):
        """Exercise the final formatter without mocking either PWSCF helper."""
        for params in (
            {"fp_params": {"ecut": 50, "ediff": 1e-6, "kspacing": 0.2}},
            {
                "user_fp_params": {
                    "control": {"calculation": "scf"},
                    "system": {"ecutwfc": 50},
                    "electrons": {"conv_thr": 1e-6},
                    "kspacing": 0.2,
                }
            },
        ):
            with self.subTest(params=params), tempfile.TemporaryDirectory() as tmp:
                task = Path(tmp)
                (task / "POSCAR").write_text(
                    "H\n1.0\n4 0 0\n0 4 0\n0 0 4\nH\n1\nDirect\n0 0 0\n"
                )
                cwd = os.getcwd()
                try:
                    relabel.make_non_vasp_input(
                        task, "pwscf", params, [1.0], tmp, ["H.UPF"]
                    )
                    self.assertEqual(os.getcwd(), cwd)
                finally:
                    os.chdir(cwd)
                text = (task / "input").read_text()
                self.assertIn("&control", text)
                self.assertIn("ecutwfc=50", text.replace(" ", ""))
                self.assertIn("H.UPF", text)
                self.assertIn("K_POINTS", text)

    @patch("dpgen.tools.relabel.make_siesta")
    def test_siesta_arguments_match_helper_signature(self, make_siesta):
        params = {"ecut": 100}
        relabel.make_non_vasp_input(
            "task",
            "siesta",
            {"fp_params": params},
            [1.0],
            "/pp",
            ["H.psf"],
        )

        make_siesta.assert_called_once_with("task", params, "/pp", ["H.psf"])


if __name__ == "__main__":
    unittest.main()

import json
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

"""Tests for carrying VASP relaxation charge densities into init MD tasks."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from dpgen.data.arginfo import init_bulk_jdata_arginfo
from dpgen.data.gen import make_vasp_md, run_vasp_md, run_vasp_relax


class TestInitChgcar(unittest.TestCase):
    """Validate local setup and remote staging for CHGCAR reuse."""

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.out_dir = self.root / "init"
        self.md_incar = self.root / "INCAR.md"
        self.potcar = self.root / "POTCAR"
        self.md_incar.write_text("ISTART = 0\nICHARG = 1\n")
        self.potcar.write_text("pseudopotential\n")
        self.jdata = {
            "out_dir": str(self.out_dir),
            "potcars": [str(self.potcar)],
            "scale": [1.0],
            "pert_numb": 0,
            "md_nstep": 1,
            "md_incar": str(self.md_incar),
            "reuse_relax_chgcar": True,
        }
        self.mdata = {
            "fp_command": "vasp_std",
            "fp_group_size": 1,
            "fp_resources": {},
            "fp_machine": {},
        }

    def tearDown(self):
        self.tempdir.cleanup()

    def _make_md_source(self):
        source = self.out_dir / "01.scale_pert" / "sys-0001" / "scale-1.000" / "000000"
        source.mkdir(parents=True)
        (source / "POSCAR").write_text("structure\n")

    def test_vasp_schema_accepts_chgcar_reuse(self):
        arginfo = init_bulk_jdata_arginfo()
        data = {
            "stages": [1, 2, 3],
            "elements": ["Al"],
            "potcars": ["POTCAR"],
            "cell_type": "fcc",
            "super_cell": [2, 2, 2],
            "relax_incar": "INCAR.relax",
            "md_incar": "INCAR.md",
            "scale": [1.0],
            "skip_relax": False,
            "pert_numb": 1,
            "pert_box": 0.01,
            "pert_atom": 0.01,
            "md_nstep": 10,
            "coll_ndata": 10,
            "init_fp_style": "VASP",
            "reuse_relax_chgcar": True,
        }

        normalized = arginfo.normalize_value(data)
        arginfo.check_value(normalized, strict=True)

    def test_make_vasp_md_links_system_relax_chgcar(self):
        self._make_md_source()
        relax_dir = self.out_dir / "00.place_ele" / "sys-0001"
        relax_dir.mkdir(parents=True)
        relax_chgcar = relax_dir / "CHGCAR"
        relax_chgcar.write_text("charge density\n")

        make_vasp_md(self.jdata, {"fp_resources": {}})

        task_chgcar = (
            self.out_dir / "02.md" / "sys-0001" / "scale-1.000" / "000000" / "CHGCAR"
        )
        self.assertTrue(task_chgcar.is_symlink())
        self.assertTrue(os.path.samefile(task_chgcar, relax_chgcar))

    def test_make_vasp_md_requires_relax_chgcar(self):
        self._make_md_source()

        with self.assertRaisesRegex(RuntimeError, "stage-1 VASP relaxation"):
            make_vasp_md(self.jdata, {"fp_resources": {}})

    @patch("dpgen.data.gen.check_api_version")
    @patch("dpgen.data.gen.make_submission")
    def test_dispatcher_transfers_chgcar(self, make_submission, _check_api_version):
        submission = MagicMock()
        make_submission.return_value = submission

        relax_task = self.out_dir / "00.place_ele" / "sys-0001"
        relax_task.mkdir(parents=True)
        run_vasp_relax(self.jdata, self.mdata)
        relax_call = make_submission.call_args.kwargs
        self.assertIn("CHGCAR", relax_call["backward_files"])

        make_submission.reset_mock()
        md_task = self.out_dir / "02.md" / "sys-0001" / "scale-1.000" / "000000"
        md_task.mkdir(parents=True)
        run_vasp_md(self.jdata, self.mdata)
        md_call = make_submission.call_args.kwargs
        self.assertIn("CHGCAR", md_call["forward_files"])


if __name__ == "__main__":
    unittest.main()

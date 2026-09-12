"""Validate the self-contained CH4 Lebesgue example."""

import json
import unittest
from pathlib import Path

from pymatgen.io.vasp import Incar, Poscar

from dpgen.data.arginfo import init_bulk_jdata_arginfo
from dpgen.generator.arginfo import run_jdata_arginfo, run_mdata_arginfo
from dpgen.util import normalize


class TestCH4LebesgueExample(unittest.TestCase):
    """Ensure the example contains usable inputs instead of path stubs."""

    example = Path(__file__).parent.parent / "examples" / "CH4-lebesgue"

    def test_json_inputs_match_schemas(self):
        """Normalize each JSON input with the schema used by its CLI command."""
        cases = (
            (init_bulk_jdata_arginfo(), "init.json"),
            (run_jdata_arginfo(), "param_CH4_deepmd-kit-2.0.1.json"),
            (run_mdata_arginfo(), "lebesgue_v2_machine.json"),
        )
        for arginfo, name in cases:
            with self.subTest(name=name):
                with (self.example / name).open() as file:
                    normalize(arginfo, json.load(file))

    def test_vasp_inputs_are_parseable(self):
        """Parse the POSCAR and all INCAR variants with pymatgen."""
        poscar = Poscar.from_file(self.example / "CH4.POSCAR")
        self.assertEqual(poscar.natoms, [4, 1])

        expected_ibrion = {
            "INCAR_methane": -1,
            "INCAR_methane.md": 0,
            "INCAR_methane.rlx": 2,
        }
        for name, ibrion in expected_ibrion.items():
            with self.subTest(name=name):
                incar = Incar.from_file(self.example / name)
                self.assertEqual(incar["IBRION"], ibrion)

    def test_readme_is_real_documentation(self):
        """Reject a regression to the former one-line README path stub."""
        readme = (self.example / "README.md").read_text()
        self.assertTrue(readme.startswith("# CH4 workflow"))

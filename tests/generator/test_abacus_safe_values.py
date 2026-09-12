import tempfile
import unittest
from pathlib import Path

from dpgen.generator.lib.abacus_scf import (
    _parse_abacus_binary,
    make_abacus_scf_input,
)


class TestAbacusSafeValues(unittest.TestCase):
    def test_binary_string_values(self):
        cases = {
            "1": 1,
            "true": 1,
            "T": 1,
            "0": 0,
            "false": 0,
            "F": 0,
        }
        for value, expected in cases.items():
            with self.subTest(value=value):
                self.assertEqual(_parse_abacus_binary(value, "option"), expected)

    def test_python_expression_is_rejected_without_execution(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            marker = Path(tmpdir) / "eval-ran"
            payload = (
                f"__import__('pathlib').Path({str(marker)!r}).write_text('ran') or 1"
            )

            with self.assertRaisesRegex(ValueError, "gamma_only"):
                make_abacus_scf_input({"gamma_only": payload})

            self.assertFalse(marker.exists())


if __name__ == "__main__":
    unittest.main()

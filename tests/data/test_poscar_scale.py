import tempfile
import unittest
from pathlib import Path

from dpgen.data.gen import poscar_scale


class TestPoscarScale(unittest.TestCase):
    def test_selective_cartesian_coordinates_and_flags(self):
        """Selective Dynamics moves the coordinate start and appends T/F flags."""
        poscar = """Mg
1.0
1 0 0
0 1 0
0 0 1
Mg
1
Selective dynamics
Cartesian
0.25 0.50 0.75 T F T
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "POSCAR.in"
            target = Path(tmpdir) / "POSCAR.out"
            source.write_text(poscar)

            poscar_scale(source, target, 2.0)
            scaled = target.read_text().splitlines()

        self.assertEqual(scaled[8], "Cartesian")
        self.assertEqual(
            scaled[9],
            "5.0000000000000000e-01 1.0000000000000000e+00 "
            "1.5000000000000000e+00 T F T",
        )


if __name__ == "__main__":
    unittest.main()

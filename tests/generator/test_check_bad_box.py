import os
import sys
import tempfile
import unittest

from ase import Atoms
from ase.io import write

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "generator"
from .context import (
    check_bad_box,
    setUpModule,  # noqa: F401
)


class TestCheckBadBox(unittest.TestCase):
    def test_length_ratio(self):
        dirname = os.path.dirname(__file__)
        conf_bad = os.path.join(dirname, "check_bad_box", "bad.length.lammpstrj")
        conf_good = os.path.join(dirname, "check_bad_box", "good.lammpstrj")
        self.assertTrue(check_bad_box(conf_bad, "length_ratio:5"))
        self.assertFalse(check_bad_box(conf_good, "length_ratio:5"))

    def test_height_ratio(self):
        dirname = os.path.dirname(__file__)
        conf_bad = os.path.join(dirname, "check_bad_box", "bad.height.POSCAR")
        self.assertTrue(check_bad_box(conf_bad, "height_ratio:5", fmt="vasp/POSCAR"))
        self.assertFalse(check_bad_box(conf_bad, "length_ratio:5", fmt="vasp/POSCAR"))

    def test_min_distance(self):
        dirname = os.path.dirname(__file__)
        conf = os.path.join(dirname, "check_bad_box", "good.lammpstrj")
        self.assertFalse(check_bad_box(conf, "min_distance:2.0"))
        self.assertTrue(check_bad_box(conf, "min_distance:2.2"))
        self.assertTrue(check_bad_box(conf, "min_dist:2.2"))

    def test_min_distance_includes_periodic_self_images(self):
        """Compressed single-atom and skewed cells have real periodic neighbors."""
        for cell, positions, expected in (
            ([0.4, 4, 4], [[0, 0, 0]], True),
            ([4, 4, 4], [[0, 0, 0]], False),
            ([1, 4, 4], [[0, 0, 0]], False),
            ([[4, 0, 0], [3.8, 0.4, 0], [0, 0, 4]], [[0, 0, 0]], True),
            ([4, 4, 4], [[0, 0, 0], [3.8, 0, 0]], True),
            ([4, 4, 4], [[0, 0, 0], [0, 0, 0]], True),
        ):
            with self.subTest(cell=cell, positions=positions):
                atoms = Atoms(
                    "H" * len(positions), positions=positions, cell=cell, pbc=True
                )
                with tempfile.TemporaryDirectory() as tmp:
                    path = os.path.join(tmp, "POSCAR")
                    write(path, atoms, format="vasp")
                    for criterion in ("min_distance:1.0", "min_dist:1.0"):
                        self.assertEqual(
                            check_bad_box(path, criterion, fmt="vasp/poscar"), expected
                        )

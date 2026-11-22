import unittest
from unittest.mock import patch

import numpy as np

from dpgen.generator.lib import gaussian


class TestGaussianInput(unittest.TestCase):
    def test_fragment_guesses_select_keywords_from_final_multiplicity(self):
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.9572, 0.0, 0.0],
                [-0.239, 0.927, 0.0],
                [5.0, 0.0, 0.0],
                [6.21, 0.0, 0.0],
            ]
        )
        atom_types = np.array([0, 1, 1, 0, 0])
        fp_params = {
            "nproc": 1,
            "keywords": "mn15/6-31g**",
            "multiplicity": "auto",
            "fragment_guesses": True,
            "keywords_high_multiplicity": ["mn15/6-31g** stable=opt"],
        }
        cases = (
            (
                slice(None),
                [0, 0, 0, 1, 1],
                "#mn15/6-31g** stable=opt guess=fragment=2",
                "0 3 0 1 0 3",
            ),
            (slice(3, None), [0, 0], "#mn15/6-31g** stable=opt", "0 3"),
            (slice(3), [0, 0, 0], "#mn15/6-31g**", "0 1"),
        )
        for atoms, fragment_indices, route, charge_multiplicity in cases:
            with self.subTest(fragment_indices=fragment_indices):
                sys_data = {
                    "atom_names": ["O", "H"],
                    "atom_types": atom_types[atoms],
                    "coords": coords[atoms][np.newaxis],
                }
                fragment_count = max(fragment_indices) + 1
                with (
                    patch.object(
                        gaussian,
                        "_crd2frag",
                        return_value=(fragment_count, np.array(fragment_indices)),
                    ),
                    patch.object(gaussian, "_crd2mul", return_value=1) as estimate,
                ):
                    lines = gaussian.make_gaussian_input(
                        sys_data, fp_params
                    ).splitlines()
                self.assertEqual(lines[1], route)
                self.assertEqual(lines[5], charge_multiplicity)
                estimate.assert_not_called()

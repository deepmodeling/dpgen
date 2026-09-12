import unittest
from copy import deepcopy

import numpy as np

from dpgen.generator.lib.cp2k import default_config, make_cp2k_input


class TestCp2kInput(unittest.TestCase):
    def test_user_config_does_not_mutate_defaults(self):
        system = {
            "cells": np.eye(3).reshape(1, 3, 3),
            "atom_names": ["H"],
            "atom_types": np.array([0]),
        }
        original_defaults = deepcopy(default_config)

        first_input = make_cp2k_input(system, {"GLOBAL": {"PROJECT": "FIRST"}})
        second_input = make_cp2k_input(system, {})

        self.assertIn("PROJECT FIRST", first_input)
        self.assertIn("PROJECT DPGEN", second_input)
        self.assertEqual(default_config, original_defaults)


if __name__ == "__main__":
    unittest.main()

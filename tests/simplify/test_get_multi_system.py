import unittest
import os
import shutil

import numpy as np
import dpdata

from context import dpgen


class TestGetMultiSystem(unittest.TestCase):
    def setUp(self) -> None:
        system = dpdata.System(data={
            "atom_names": ["H"],
            "atom_numbs": [1],
            "atom_types": np.zeros((1,), dtype=int),
            "coords": np.zeros((1, 1, 3), dtype=np.float32),
            "cells": np.zeros((1, 3, 3), dtype=np.float32),
            "orig": np.zeros(3, dtype=np.float32),
            "nopbc": True,
            "energies": np.zeros((1,), dtype=np.float32),
            "forces": np.zeros((1, 1, 3), dtype=np.float32),
        })
        system.to_deepmd_npy("data0")
        system.to_deepmd_npy("data1")
        system.to_deepmd_hdf5("data2.hdf5")
        self.data = [
            "data0",
            "data1",
            "data2.hdf5",
        ]

    def tearDown(self) -> None:
        for dd in self.data:
            if dd.endswith(".hdf5"):
                os.remove(dd)
            else:
                shutil.rmtree(dd)

    def test_get_multi_system(self):
        for list_data in (True, False):
            for labeled in (True, False):
                with self.subTest(list_data=list_data, labeled=labeled):
                    ms = dpgen.simplify.simplify.get_multi_system(
                        self.data if list_data else self.data[0],
                        {"labeled": labeled},
                    )
                    assert isinstance(ms, dpdata.MultiSystems)
                    for ss in ms.systems.values():
                        assert isinstance(ss, dpdata.LabeledSystem if labeled else dpdata.System)
                    assert ms.get_nframes() == len(self.data) if list_data else 1

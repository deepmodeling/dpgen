import unittest

import numpy as np

from dpgen.auto_test.Interstitial import _smallest_nonzero_distance


class TestInterstitialDistance(unittest.TestCase):
    def test_smallest_nonzero_distance(self):
        distance_matrix = np.array(
            [
                [0.0, 2.5, 1.25],
                [2.5, 0.0, 3.0],
                [1.25, 3.0, 0.0],
            ]
        )
        self.assertEqual(_smallest_nonzero_distance(distance_matrix), 1.25)


if __name__ == "__main__":
    unittest.main()

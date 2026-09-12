import unittest
from unittest.mock import patch

import numpy as np

from dpgen.data.tools.create_random_disturb import gen_random_disturb


class TestGenRandomDisturb(unittest.TestCase):
    def test_normal_disturbance_uses_scaled_normal_variate(self):
        with (
            patch(
                "dpgen.data.tools.create_random_disturb.np.random.rand",
                return_value=np.array([1.0, 0.5, 0.5]),
            ),
            patch(
                "dpgen.data.tools.create_random_disturb.np.random.normal",
                return_value=0.25,
            ) as normal,
        ):
            displacement = gen_random_disturb(2.0, -0.5, 0.5, dstyle="normal")

        normal.assert_called_once_with(loc=0.0, scale=0.5)
        np.testing.assert_allclose(displacement, [0.5, 0.0, 0.0])


if __name__ == "__main__":
    unittest.main()

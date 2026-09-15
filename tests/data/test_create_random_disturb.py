import unittest
from unittest.mock import call, patch

import numpy as np

from dpgen.data.tools.create_random_disturb import (
    create_disturbs_atomsk,
    gen_random_disturb,
)


class TestGenRandomDisturb(unittest.TestCase):
    @patch("dpgen.data.tools.create_random_disturb.subprocess.run")
    @patch("dpgen.data.tools.create_random_disturb.glob.glob", return_value=[])
    def test_atomsk_uses_argument_list_and_checks_exit_status(self, _glob, run):
        create_disturbs_atomsk("input file", 2, dmax=0.2, ofmt="lmp")

        self.assertEqual(
            run.call_args_list,
            [
                call(
                    [
                        "atomsk",
                        "input file",
                        "-disturb",
                        "0.2",
                        "-wrap",
                        "-ow",
                        "input file1.lmp",
                    ],
                    check=True,
                ),
                call(
                    [
                        "atomsk",
                        "input file",
                        "-disturb",
                        "0.2",
                        "-wrap",
                        "-ow",
                        "input file2.lmp",
                    ],
                    check=True,
                ),
            ],
        )

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

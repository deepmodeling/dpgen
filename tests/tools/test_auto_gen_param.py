import unittest

from dpgen.tools.auto_gen_param import Iteration


class TestIterationIndex(unittest.TestCase):
    def test_index_iteration_round_trip(self):
        iteration = Iteration([50])
        registered_index = iteration.index_iteration

        self.assertIsInstance(registered_index, int)
        iteration.index_iteration = 7
        self.assertEqual(iteration.index_iteration, 7)


if __name__ == "__main__":
    unittest.main()

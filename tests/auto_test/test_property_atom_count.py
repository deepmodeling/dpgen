import unittest

from dpgen.auto_test.Property import _total_atom_count


class TestPropertyAtomCount(unittest.TestCase):
    def test_counts_all_species(self):
        self.assertEqual(_total_atom_count({"atom_numbs": [2, 3, 4]}), 9)


if __name__ == "__main__":
    unittest.main()

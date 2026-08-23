import unittest
from unittest.mock import patch

from dpgen.tools.collect_data import collect_data


class TestLegacyCollectData(unittest.TestCase):
    @patch("dpgen.tools.collect_data.collect_current_data", return_value="result")
    def test_delegates_to_maintained_collector(self, current_collect_data):
        with self.assertWarnsRegex(DeprecationWarning, "dpgen collect"):
            result = collect_data("job", "param.json", "output", verbose=False)

        self.assertEqual(result, "result")
        current_collect_data.assert_called_once_with(
            "job",
            "param.json",
            "output",
            verbose=False,
            shuffle=True,
            merge=False,
        )


if __name__ == "__main__":
    unittest.main()

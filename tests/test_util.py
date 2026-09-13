import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dpgen.util import convert_training_data_to_hdf5


class TestConvertTrainingDataToHdf5(unittest.TestCase):
    def test_rewritten_json_is_truncated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_file = tmp_path / "input.json"
            input_file.write_text(
                json.dumps(
                    {
                        "training": {"systems": []},
                        "padding": "content that must not survive the rewrite",
                    }
                )
            )

            with patch(
                "dpgen.util.json.dump",
                side_effect=lambda _data, file, indent: file.write("{}"),
            ):
                convert_training_data_to_hdf5(
                    [str(input_file)], str(tmp_path / "data.hdf5")
                )

            self.assertEqual(input_file.read_text(), "{}")


if __name__ == "__main__":
    unittest.main()

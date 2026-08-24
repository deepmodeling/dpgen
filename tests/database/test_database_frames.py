import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from dpgen.database.run import _parsing_vasp


class TestVaspDatabaseFrames(unittest.TestCase):
    @patch("dpgen.database.run.Entry", side_effect=lambda *args, **kwargs: kwargs)
    @patch("dpgen.database.run.LabeledSystem")
    @patch("dpgen.database.run.VaspInput.from_directory")
    def test_each_outcar_frame_becomes_an_entry(
        self, from_directory, labeled_system, _entry
    ):
        vasp_input = MagicMock()
        vasp_input.__getitem__.return_value = SimpleNamespace(
            structure=SimpleNamespace(composition="H2")
        )
        vasp_input.as_dict.return_value = {"input": "vasp"}
        from_directory.return_value = vasp_input

        frames = [MagicMock(), MagicMock()]
        frames[0].as_dict.return_value = {"frame": 0}
        frames[1].as_dict.return_value = {"frame": 1}
        labeled_system.return_value.to_list.return_value = frames

        entries = _parsing_vasp(
            ["init/system/02.md/task.000"], {}, "dpgen", iters=False
        )

        self.assertEqual(len(entries), 2)
        self.assertEqual(
            [entry["entry_id"] for entry in entries], ["dpgen_0", "dpgen_1"]
        )


if __name__ == "__main__":
    unittest.main()

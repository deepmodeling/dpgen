import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, mock_open, patch

from dpgen.database.run import _parsing_vasp, parsing_vasp


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

    @patch("dpgen.database.run.dumpfn")
    @patch("builtins.open", new_callable=mock_open)
    @patch("dpgen.database.run.glob")
    @patch("dpgen.database.run.Entry", side_effect=lambda *args, **kwargs: kwargs)
    @patch("dpgen.database.run.LabeledSystem")
    @patch("dpgen.database.run.VaspInput.from_directory")
    def test_combined_initial_and_iteration_ids_are_unique(
        self, from_directory, labeled_system, _entry, glob_paths, _open, dump
    ):
        """The public collector must share one ID sequence across both sources."""
        vasp_input = MagicMock()
        from_directory.return_value = vasp_input
        frames = [MagicMock(), MagicMock()]
        labeled_system.return_value.to_list.return_value = frames
        glob_paths.side_effect = [
            ["run/iter.000000/02.fp/task.000.000000"],
            ["init/system/02.md/task.000"],
        ]
        parsing_vasp("run", {}, skip_init=False, id_prefix="dpgen")
        entries = dump.call_args.args[0]
        self.assertEqual(
            [entry["entry_id"] for entry in entries],
            ["dpgen_0", "dpgen_1", "dpgen_2", "dpgen_3"],
        )


if __name__ == "__main__":
    unittest.main()

import unittest
from unittest.mock import MagicMock, patch

from dpgen.auto_test.common_prop import run_property


class TestRunProperty(unittest.TestCase):
    @patch(
        "dpgen.auto_test.common_prop.convert_mdata", side_effect=lambda data, _: data
    )
    @patch("dpgen.auto_test.common_prop.make_calculator")
    @patch("dpgen.auto_test.common_prop.util.collect_task")
    @patch("dpgen.auto_test.common_prop.glob.glob")
    @patch("dpgen.auto_test.common_prop.Pool")
    def test_only_unfinished_tasks_are_sent_to_worker(
        self, pool_class, glob_files, collect_task, make_calculator, _convert_mdata
    ):
        task_paths = ["/work/task.000000", "/work/task.000001"]
        glob_files.side_effect = lambda pattern: (
            task_paths if "task.[0-9]*[0-9]" in pattern else ["conf"]
        )
        collect_task.return_value = ["task.000001"]

        calculator = MagicMock()
        calculator.forward_files.return_value = []
        calculator.forward_common_files.return_value = []
        calculator.backward_files.return_value = []
        make_calculator.return_value = calculator

        pool = MagicMock()
        result = MagicMock()
        result.successful.return_value = True
        pool.apply_async.return_value = result
        pool_class.return_value = pool

        run_property(
            ["conf"],
            {"type": "vasp"},
            [{"type": "eos"}],
            {},
        )

        worker_args = pool.apply_async.call_args.args[1]
        self.assertEqual(worker_args[1], ["task.000001"])


if __name__ == "__main__":
    unittest.main()

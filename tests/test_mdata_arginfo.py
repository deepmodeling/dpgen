import copy
import unittest

from dargs.dargs import ArgumentTypeError, ArgumentValueError

from dpgen.generator.arginfo import run_mdata_arginfo
from dpgen.remote.decide_machine import convert_mdata
from dpgen.util import normalize


def _task(command):
    return {
        "command": command,
        "machine": {
            "batch_type": "Shell",
            "context_type": "local",
            "local_root": "./",
            "remote_root": "/tmp/dpgen",
        },
        "resources": {
            "batch_type": "Shell",
            "number_node": 1,
            "cpu_per_node": 4,
            "gpu_per_node": 0,
            "group_size": 1,
        },
    }


def _mdata(section_wrapper):
    return {
        "api_version": "1.0",
        "train": section_wrapper(_task("dp")),
        "model_devi": section_wrapper(_task("lmp")),
        "fp": section_wrapper(_task("vasp")),
    }


class TestRunMdataArginfo(unittest.TestCase):
    def test_accepts_dict_sections(self):
        data = _mdata(lambda item: item)
        normalize(run_mdata_arginfo(), data, strict_check=False)

    def test_accepts_list_sections_used_by_runtime(self):
        data = _mdata(lambda item: [item])
        normalize(run_mdata_arginfo(), data, strict_check=False)

        converted = convert_mdata(copy.deepcopy(data))
        self.assertEqual(converted["train_command"], "dp")
        self.assertEqual(converted["model_devi_command"], "lmp")
        self.assertEqual(converted["fp_command"], "vasp")

    def test_rejects_list_with_non_dict_item(self):
        data = _mdata(lambda item: item)
        data["train"] = ["not-a-task-dict"]
        with self.assertRaises(ArgumentValueError):
            normalize(run_mdata_arginfo(), data, strict_check=False)

    def test_rejects_empty_list_section(self):
        data = _mdata(lambda item: item)
        data["train"] = []
        with self.assertRaises(ArgumentValueError):
            normalize(run_mdata_arginfo(), data, strict_check=False)

    def test_rejects_scalar_section(self):
        data = _mdata(lambda item: item)
        data["train"] = "not-a-task-section"
        with self.assertRaises(ArgumentTypeError):
            normalize(run_mdata_arginfo(), data, strict_check=False)

import unittest
from pathlib import Path

from dpgen.util import load_file

this_directory = Path(__file__).parent


class TestLoadFile(unittest.TestCase):
    def test_load_json_file(self):
        ref = {"aa": "bb"}
        jdata = load_file(this_directory / "sample.json")
        self.assertEqual(jdata, ref)

    def test_load_yaml_file(self):
        ref = {"aa": "bb"}
        jdata = load_file(this_directory / "sample.yaml")
        self.assertEqual(jdata, ref)

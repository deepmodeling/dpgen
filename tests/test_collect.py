import json
import unittest
import tempfile
import dpdata

from pathlib import Path
from dpgen.collect.collect import collect_data


class TestCollectData(unittest.TestCase):
    def setUp(self):
        self.data = dpdata.LabeledSystem(Path(__file__).parent / "generator" / "data" / "deepmd", fmt="deepmd/npy")

    def test_collect_data(self):
        with tempfile.TemporaryDirectory() as inpdir, \
             tempfile.TemporaryDirectory() as outdir, \
             tempfile.NamedTemporaryFile() as param_file:
            self.data.to_deepmd_npy(Path(inpdir) / "iter.000000" / "02.fp" / "data.000")
            self.data.to_deepmd_npy(Path(inpdir) / "iter.000001" / "02.fp" / "data.000" / "aa")
            self.data.to_deepmd_npy(Path(inpdir) / "iter.000001" / "02.fp" / "data.000" / "bb")
            with open(param_file.name, "w") as fp:
                json.dump({"sys_configs": ["sys1"], "model_devi_jobs": [{}, {}, {}]}, fp)

            collect_data(inpdir, param_file.name, outdir, verbose=True)
            ms = dpdata.MultiSystems().from_deepmd_npy(outdir)
            self.assertEqual(ms.get_nframes(), self.data.get_nframes() * 3)

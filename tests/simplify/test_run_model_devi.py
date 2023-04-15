import shutil
import textwrap
import unittest
import tempfile
import dpdata
from pathlib import Path
from context import dpgen


class TestOneH5(unittest.TestCase):
    def setUp(self):
        work_path = Path("iter.000000") / "01.model_devi"
        work_path.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(Path(tmpdir) / "test.xyz", 'w') as f:
                f.write(textwrap.dedent("""\
                    2

                    H 0.0 0.0 0.0
                    H 0.0 0.0 1.0
                    """))
            dpdata.MultiSystems(dpdata.System(Path(tmpdir) / "test.xyz", fmt="xyz")).to_deepmd_npy(work_path / (dpgen.simplify.simplify.rest_data_name + ".old"))

    def tearDown(self) -> None:
        shutil.rmtree("iter.000000")

    def test_npy(self):
        jdata = {
            "type_map": ["H"],
        }
        with tempfile.TemporaryDirectory() as remote_root:
            mdata = {
                "model_devi_command": (
                    f"test -d {dpgen.simplify.simplify.rest_data_name}.old"
                    f"&& touch {dpgen.simplify.simplify.detail_file_name_prefix}"
                    "&& echo dp"
                ),
                "model_devi_machine": {
                    "context_type": "LocalContext",
                    "batch_type": "shell",
                    "local_root": "./",
                    "remote_root": remote_root,
                },
                "model_devi_resources": {
                    "group_size": 1,
                }
            }
            dpgen.simplify.simplify.run_model_devi(0, jdata=jdata, mdata=mdata)

    def test_one_h5(self):
        jdata = {
            "type_map": ["H"],
            "one_h5": True,
        }
        with tempfile.TemporaryDirectory() as remote_root:
            mdata = {
                "model_devi_command": (
                    f"test -f {dpgen.simplify.simplify.rest_data_name}.old.hdf5"
                    f"&& touch {dpgen.simplify.simplify.detail_file_name_prefix}"
                    "&& echo dp"
                ),
                "model_devi_machine": {
                    "context_type": "LocalContext",
                    "batch_type": "shell",
                    "local_root": "./",
                    "remote_root": remote_root,
                },
                "model_devi_resources": {
                    "group_size": 1,
                }
            }
            dpgen.simplify.simplify.run_model_devi(0, jdata=jdata, mdata=mdata)

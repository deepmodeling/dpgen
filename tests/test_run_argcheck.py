import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from dargs.dargs import ArgumentKeyError

from dpgen.generator.run import run_iter


class TestRunArgcheck(unittest.TestCase):
    def test_run_rejects_unknown_parameter(self):
        repository_root = Path(__file__).resolve().parents[1]
        example = (
            repository_root
            / "examples"
            / "run"
            / "dp2.x-lammps-cp2k"
            / "param_CH4_deepmd-kit-2.0.1.json"
        )
        parameter_data = json.loads(example.read_text())
        parameter_data["unknown_parameter"] = True

        with TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            param_file = tmpdir / "param.json"
            machine_file = tmpdir / "machine.json"
            param_file.write_text(json.dumps(parameter_data))
            machine_file.write_text("{}")

            with self.assertRaises(ArgumentKeyError):
                run_iter(str(param_file), str(machine_file))


if __name__ == "__main__":
    unittest.main()

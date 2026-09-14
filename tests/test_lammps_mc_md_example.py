import json
import unittest
from pathlib import Path


class TestLammpsMcMdExample(unittest.TestCase):
    def test_template_placeholders_and_dump_command(self):
        """The example must become runnable after applying its revision matrix."""
        example = Path(__file__).parents[1] / "examples" / "run" / "dp-lammps-MC+MD"
        template = (example / "template" / "lammps.in").read_text()
        parameters = json.loads((example / "param.json").read_text())
        replacements = parameters["model_devi_jobs"][0]["rev_mat"]["lmp"]

        rendered = template
        for placeholder, values in replacements.items():
            self.assertIn(placeholder, template)
            rendered = rendered.replace(placeholder, str(values[0]))

        variable_lines = [
            line for line in rendered.splitlines() if line.startswith("variable")
        ]
        self.assertFalse(any("${" in line for line in variable_lines))
        self.assertNotIn('RANDOM_SEED}"', rendered)

        dump_line = next(
            line for line in template.splitlines() if line.startswith("dump")
        )
        self.assertGreaterEqual(len(dump_line.split()), 9)
        self.assertIn("traj/*.lammpstrj", dump_line)


if __name__ == "__main__":
    unittest.main()

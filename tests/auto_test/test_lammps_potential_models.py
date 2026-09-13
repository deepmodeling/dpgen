import unittest

from dpgen.auto_test.Lammps import Lammps
from dpgen.auto_test.lib.lammps import inter_eam_alloy, inter_meam


class TestLammpsPotentialModels(unittest.TestCase):
    def test_meam_model_parameters_generate_two_file_command(self):
        interaction = {
            "type": "meam",
            "model": ["lammps_input/meam.lib", "lammps_input/Al.meam"],
            "type_map": {"Al": 0},
        }
        calculator = Lammps(interaction, "unused.vasp")
        calculator.set_model_param()

        self.assertEqual(calculator.model_param["model_name"], ["meam.lib", "Al.meam"])
        command = inter_meam(calculator.model_param)
        self.assertIn("pair_coeff      * * meam.lib Al Al.meam Al", command)

    def test_eam_alloy_command_uses_plain_filename(self):
        interaction = {
            "type": "eam_alloy",
            "model": "lammps_input/Al.eam.alloy",
            "type_map": {"Al": 0},
        }
        calculator = Lammps(interaction, "unused.vasp")
        calculator.set_model_param()

        command = inter_eam_alloy(calculator.model_param)
        self.assertIn("pair_coeff      * * Al.eam.alloy Al", command)
        self.assertNotIn("['Al.eam.alloy']", command)


if __name__ == "__main__":
    unittest.main()

import unittest
import sys
import os

# Add the parent directory to sys.path to ensure imports work
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "auto_test"

from dpgen.auto_test.VASP import VASP
from dpgen.auto_test.ABACUS import ABACUS

from .context import setUpModule  # noqa: F401


class TestForwardCommonFiles(unittest.TestCase):
    """Test that forward_common_files returns appropriate files for different calculators."""

    def test_vasp_forward_common_files_empty(self):
        """Test that VASP forward_common_files returns empty list.
        
        This test ensures that VASP doesn't return INCAR/POTCAR as common files
        since they are created per-task and symlinked, not as true common files
        in the work_path. This prevents upload errors in dpdispatcher.
        """
        inter_parameter = {
            "type": "vasp",
            "incar": "vasp_input/INCAR.rlx",
            "potcar_prefix": ".",
            "potcars": {"Li": "vasp_input/POTCAR"},
        }
        
        vasp_calc = VASP(inter_parameter, "POSCAR")
        
        # Test different property types
        property_types = ["relaxation", "static", "elastic", "vacancy", "interstitial"]
        
        for prop_type in property_types:
            with self.subTest(property_type=prop_type):
                common_files = vasp_calc.forward_common_files(prop_type)
                self.assertEqual(common_files, [], 
                               f"VASP forward_common_files should return empty list for {prop_type}")

    def test_abacus_forward_common_files_consistency(self):
        """Test that ABACUS also returns empty list, showing consistency."""
        inter_parameter = {
            "type": "abacus",
            "potcar_prefix": ".", 
            "potcars": {"Al": "POT_Al"},
            "orb_files": {"Al": "Al_gga_7au_60Ry_2s2p1d.orb"},
            "dpks_descriptor": "jle.dat"
        }
        
        abacus_calc = ABACUS(inter_parameter, "POSCAR")
        common_files = abacus_calc.forward_common_files("relaxation")
        
        self.assertEqual(common_files, [], 
                        "ABACUS forward_common_files should return empty list")

    def test_vasp_forward_files_still_works(self):
        """Test that forward_files still returns required files for each task."""
        inter_parameter = {
            "type": "vasp",
            "incar": "vasp_input/INCAR.rlx",
            "potcar_prefix": ".",
            "potcars": {"Li": "vasp_input/POTCAR"},
        }
        
        vasp_calc = VASP(inter_parameter, "POSCAR")
        forward_files = vasp_calc.forward_files("relaxation")
        
        expected_files = ["INCAR", "POSCAR", "KPOINTS", "POTCAR"]
        self.assertEqual(forward_files, expected_files,
                        "VASP forward_files should still return task-specific files")


if __name__ == "__main__":
    unittest.main()
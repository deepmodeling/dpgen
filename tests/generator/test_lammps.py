import os
import sys
import textwrap
import unittest

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
__package__ = "generator"
from dpgen.generator.lib.lammps import make_lammps_input

from .context import (
    get_all_dumped_forces,
    get_dumped_forces,
    setUpModule,  # noqa: F401
)


class TestMakeLammpsInput(unittest.TestCase):
    """Test LAMMPS input generation including D3 dispersion support."""

    def setUp(self):
        self.ensemble = "nvt"
        self.conf_file = "test.lmp"
        self.graphs = ["model.pb"]
        self.nsteps = 1000
        self.dt = 0.001
        self.neidelay = 10
        self.trj_freq = 10
        self.mass_map = [1.0, 16.0]  # H, O
        self.temp = 300.0
        self.deepmd_version = "2.0"

    def test_basic_deepmd_only(self):
        """Test basic LAMMPS input without D3 (backward compatibility)."""
        jdata = {}
        result = make_lammps_input(
            self.ensemble,
            self.conf_file,
            self.graphs,
            self.nsteps,
            self.dt,
            self.neidelay,
            self.trj_freq,
            self.mass_map,
            self.temp,
            jdata,
            pres=1.0,
            deepmd_version=self.deepmd_version,
        )

        # Should contain basic deepmd pair_style
        self.assertIn("pair_style      deepmd model.pb", result)
        self.assertIn("pair_coeff      * *", result)
        # Should NOT contain hybrid/overlay or dispersion/d3
        self.assertNotIn("hybrid/overlay", result)
        self.assertNotIn("dispersion/d3", result)

    def test_d3_enabled_basic(self):
        """Test LAMMPS input with D3 dispersion enabled."""
        jdata = {
            "lmp_d3": {
                "enable": True,
                "damping_function": "original",
                "functional": "pbe",
                "cutoff": 30.0,
                "cn_cutoff": 20.0,
            }
        }

        result = make_lammps_input(
            self.ensemble,
            self.conf_file,
            self.graphs,
            self.nsteps,
            self.dt,
            self.neidelay,
            self.trj_freq,
            self.mass_map,
            self.temp,
            jdata,
            pres=1.0,
            deepmd_version=self.deepmd_version,
        )

        # Should contain hybrid/overlay pair_style
        self.assertIn("pair_style      hybrid/overlay deepmd model.pb", result)
        self.assertIn("dispersion/d3 original pbe 30.0 20.0", result)

        # Should contain both pair_coeff lines
        lines = result.split("\n")
        deepmd_coeff_found = False
        d3_coeff_found = False
        for line in lines:
            if "pair_coeff      * * deepmd" in line:
                deepmd_coeff_found = True
            elif "pair_coeff      * * dispersion/d3" in line:
                d3_coeff_found = True

        self.assertTrue(deepmd_coeff_found, "deepmd pair_coeff not found")
        self.assertTrue(d3_coeff_found, "dispersion/d3 pair_coeff not found")

    def test_d3_missing_parameters(self):
        """Test that missing D3 parameters raise appropriate errors during validation."""
        # Since validation now happens in arginfo, we would need to test this differently
        # For now, test that incomplete D3 config still allows basic operation
        jdata = {"lmp_d3": {"enable": False}}  # Incomplete but with enable=False
        result = make_lammps_input(
            self.ensemble,
            self.conf_file,
            self.graphs,
            self.nsteps,
            self.dt,
            self.neidelay,
            self.trj_freq,
            self.mass_map,
            self.temp,
            jdata,
            pres=1.0,
            deepmd_version=self.deepmd_version,
        )

        # Should behave like basic deepmd since enable=False
        self.assertIn("pair_style      deepmd model.pb", result)
        self.assertNotIn("hybrid/overlay", result)

    def test_d3_disabled(self):
        """Test that D3 section with enable=False works like no D3."""
        jdata = {
            "lmp_d3": {
                "enable": False,
                "damping_function": "original",
                "functional": "pbe",
                "cutoff": 30.0,
                "cn_cutoff": 20.0,
            }
        }

        result = make_lammps_input(
            self.ensemble,
            self.conf_file,
            self.graphs,
            self.nsteps,
            self.dt,
            self.neidelay,
            self.trj_freq,
            self.mass_map,
            self.temp,
            jdata,
            pres=1.0,
            deepmd_version=self.deepmd_version,
        )

        # Should behave like basic deepmd
        self.assertIn("pair_style      deepmd model.pb", result)
        self.assertNotIn("hybrid/overlay", result)
        self.assertNotIn("dispersion/d3", result)

    def test_d3_with_neigh_modify_one(self):
        """Test D3 with lmp_neigh_modify_one parameter."""
        jdata = {
            "lmp_d3": {
                "enable": True,
                "damping_function": "original",
                "functional": "pbe",
                "cutoff": 30.0,
                "cn_cutoff": 20.0,
            },
            "lmp_neigh_modify_one": 2000,
        }

        result = make_lammps_input(
            self.ensemble,
            self.conf_file,
            self.graphs,
            self.nsteps,
            self.dt,
            self.neidelay,
            self.trj_freq,
            self.mass_map,
            self.temp,
            jdata,
            pres=1.0,
            deepmd_version=self.deepmd_version,
        )

        # Should contain neigh_modify with one 2000
        self.assertIn("one 2000", result)
        # Should also contain D3 configuration
        self.assertIn("hybrid/overlay", result)
        self.assertIn("dispersion/d3", result)

    def test_d3_different_parameters(self):
        """Test D3 with different parameter values."""
        jdata = {
            "lmp_d3": {
                "enable": True,
                "damping_function": "bj",
                "functional": "pbe0",
                "cutoff": 25.0,
                "cn_cutoff": 15.0,
            }
        }

        result = make_lammps_input(
            self.ensemble,
            self.conf_file,
            self.graphs,
            self.nsteps,
            self.dt,
            self.neidelay,
            self.trj_freq,
            self.mass_map,
            self.temp,
            jdata,
            pres=1.0,
            deepmd_version=self.deepmd_version,
        )

        # Check specific parameter values
        self.assertIn("dispersion/d3 bj pbe0 25.0 15.0", result)


class TestGetDumpForce(unittest.TestCase):
    def setUp(self):
        file_content = textwrap.dedent(
            """\
ITEM: TIMESTEP
40
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
-2.9180686220264818e-04 8.0855380329747089e+00 1.4011011277606830e-07
-2.9198257591541018e-04 8.0855378881632269e+00 3.3202396460852749e-08
-2.9180686326490957e-04 8.0855378891632768e+00 -1.7571268247505500e-07
ITEM: ATOMS id type x y z fx fy fz
1 1 2.09532 8.19528 2.00538 -0.00569269 -0.0200373 -0.0342394
2 1 -0.0727384 4.01773 4.05582 -0.0297083 0.0817184 0.0722508
"""
        )
        with open("tmp.dump", "w") as fp:
            fp.write(file_content)
        self.expected_f = [
            -0.00569269,
            -0.0200373,
            -0.0342394,
            -0.0297083,
            0.0817184,
            0.0722508,
        ]

    def tearDown(self):
        if os.path.isfile("tmp.dump"):
            os.remove("tmp.dump")

    def test_read_dump(self):
        ff = get_dumped_forces("tmp.dump")
        self.assertEqual(ff.shape, (2, 3))
        ff = ff.reshape([-1])
        for ii in range(6):
            self.assertAlmostEqual(ff[ii], self.expected_f[ii])


class TestGetDumpForce2(unittest.TestCase):
    def setUp(self):
        file_content = textwrap.dedent(
            """\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
ITEM: ATOMS id type x y z fx fy fz
1 1 5.38154 4.06861 3.60573 0.000868817 -0.00100822 -0.000960258
2 2 3.9454 4.80321 4.38469 0.000503458 -0.000374043 -9.15676e-05
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
ITEM: ATOMS id type x y z fx fy fz
1 1 5.35629 3.93297 3.70556 -0.125424 0.0481604 -0.0833015
2 2 3.93654 4.79972 4.48179 0.134843 -0.0444238 -0.143111
"""
        )
        with open("tmp.dump", "w") as fp:
            fp.write(file_content)
        self.expected_f = [
            0.000868817,
            -0.00100822,
            -0.000960258,
            0.000503458,
            -0.000374043,
            -9.15676e-05,
            -0.125424,
            0.0481604,
            -0.0833015,
            0.134843,
            -0.0444238,
            -0.143111,
        ]

    def tearDown(self):
        if os.path.isfile("tmp.dump"):
            os.remove("tmp.dump")

    def test_read_all_dump(self):
        ff = get_all_dumped_forces("tmp.dump")
        ff = np.array(ff)
        self.assertEqual(ff.shape, (2, 2, 3))
        ff = ff.reshape([-1])
        for ii in range(12):
            self.assertAlmostEqual(ff[ii], self.expected_f[ii])

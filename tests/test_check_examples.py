"""This module ensures input in the examples directory
could pass the argument checking.
"""
import unittest
import json
from pathlib import Path

from dpgen.util import normalize
from dpgen.data.arginfo import (
    init_reaction_jdata_arginfo,
)
from dpgen.simplify.arginfo import (
    simplify_jdata_arginfo,
)
from dpgen.generator.arginfo import (
    run_jdata_arginfo,
)

init_reaction_jdata = init_reaction_jdata_arginfo()
simplify_jdata = simplify_jdata_arginfo()
run_jdata = run_jdata_arginfo()

# directory of examples
p_examples = Path(__file__).parent.parent / "examples"

# input_files : tuple[tuple[Argument, Path]]
#   tuple of example list
input_files = (
    (init_reaction_jdata, p_examples / "init" / "reaction.json"),
    (simplify_jdata, p_examples / "simplify" / "qm7.json"),
    (simplify_jdata, p_examples / "simplify-MAPbI3-scan-lebesgue" / "simplify_example" / "simplify.json"),
    (run_jdata, p_examples / "run" / "dp2.x-lammps-vasp" / "param_CH4_deepmd-kit-2.0.1.json"),
    (run_jdata, p_examples / "run" / "dp2.x-lammps-cp2k" / "param_CH4_deepmd-kit-2.0.1.json"),
    #(run_jdata, p_examples / "run" / "dp2.x-gromacs-gaussian" / "param.json"),
    (run_jdata, p_examples / "run" / "dp1.x-lammps-vasp" / "CH4" / "param_CH4_deepmd-kit-1.1.0.json"),
    (run_jdata, p_examples / "run" / "dp1.x-lammps-vasp" / "Al" / "param_al_all_gpu-deepmd-kit-1.1.0.json"),
    (run_jdata, p_examples / "run" / "dp1.x-lammps-vasp-et" / "param_elet.json"),
    (run_jdata, p_examples / "run" / "dp2.x-lammps-ABACUS-lcao" / "fcc-al" / "run_param.json"),
    (run_jdata, p_examples / "run" / "dp2.x-lammps-ABACUS-pw" / "fcc-al" / "run_param.json"),
    #(run_jdata, p_examples / "run" / "dp1.x-lammps-cp2k" / "methane" / "param-ch4.json"),
    #(run_jdata, p_examples / "run" / "dp1.x-lammps-ABACUS-pw" / "methane" / "param.json"),
    #(run_jdata, p_examples / "run" / "dp1.x-lammps-ABACUS-lcao-dpks" / "methane" / "param.json"),
    #(run_jdata, p_examples / "run" / "dp1.x_lammps_gaussian" / "dodecane" / "dodecane.json"),
    #(run_jdata, p_examples / "run" / "dp-lammps-enhance_sampling" / "param.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "param-mg-vasp.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "param-mg-vasp-ucloud.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "param-pyridine-pwscf.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "param-h2oscan-vasp.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "dp0.12-lammps-cp2k" / "CH4" / "param_CH4.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "dp0.12-lammps-pwmat" / "param_CH4.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "dp0.12-lammps-siesta" / "dp-lammps-siesta" / "CH4" / "param_CH4.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "dp0.12-lammps-vasp" / "Al" / "param_al_all_gpu.json"),
    #(run_jdata, p_examples / "run" / "deprecated" / "dp0.12-lammps-vasp" / "CH4" / "param_CH4.json"),
)


class TestExamples(unittest.TestCase):
    def test_arguments(self):
        for arginfo, fn in input_files:
            fn = str(fn)
            with self.subTest(fn=fn):
                with open(fn) as f:
                    data = json.load(f)
                normalize(arginfo, data)

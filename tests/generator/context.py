import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from dpgen.generator.lib.ele_temp import NBandsEsti  # noqa: F401
from dpgen.generator.lib.gaussian import _crd2frag, detect_multiplicity  # noqa: F401
from dpgen.generator.lib.lammps import (
    get_all_dumped_forces,  # noqa: F401
    get_dumped_forces,  # noqa: F401
)
from dpgen.generator.lib.make_calypso import (
    make_calypso_input,  # noqa: F401
    write_model_devi_out,  # noqa: F401
)
from dpgen.generator.lib.parse_calypso import (
    _parse_calypso_dis_mtx,  # noqa: F401
    _parse_calypso_input,  # noqa: F401
)
from dpgen.generator.run import *  # noqa: F403
from dpgen.util import setup_ele_temp  # noqa: F401

param_file = "param-mg-vasp.json"
param_pimd_file = "param-mg-pimd-vasp.json"
param_file_merge_traj = "param-mg-vasp_merge_traj.json"
param_file_v1 = "param-mg-vasp-v1.json"
param_file_v1_et = "param-mg-vasp-v1-et.json"
param_old_file = "param-mg-vasp-old.json"
param_pwscf_file = "param-pyridine-pwscf.json"
param_pwscf_old_file = "param-pyridine-pwscf-old.json"
param_gaussian_file = "param-pyridine-gaussian.json"
param_siesta_file = "param-pyridine-siesta.json"
param_cp2k_file = "param-pyridine-cp2k.json"
param_cp2k_file_exinput = "param-mgo-cp2k-exinput.json"
ref_cp2k_file_input = "cp2k_test_ref.inp"
ref_cp2k_file_exinput = "cp2k_test_exref.inp"
machine_file = "machine-local.json"
machine_file_v1 = "machine-local-v1.json"
param_diy_file = "param-mg-vasp-diy.json"
param_pwmat_file = "param-pyridine-pwmat.json"
param_abacus_file = "param-pyridine-abacus.json"
param_abacus_post_file = "param-methane-abacus.json"
param_diy_abacus_post_file = "param-methane-abacus-diy.json"
param_amber_file = "param-amber.json"
param_multiple_trust_file = "param-mg-vasp-multi-trust.json"
param_custom_fp_file = "param-custom-fp.json"


def my_file_cmp(test, f0, f1):
    with open(f0) as fp0:
        with open(f1) as fp1:
            test.assertTrue(fp0.read() == fp1.read())


def setUpModule():
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

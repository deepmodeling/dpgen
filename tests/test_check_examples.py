"""This module ensures input in the examples directory
could pass the argument checking.
"""

import json
import unittest
from pathlib import Path

from dargs import Argument

from dpgen.data.arginfo import (
    init_bulk_jdata_arginfo,
    init_reaction_jdata_arginfo,
    init_surf_jdata_arginfo,
)
from dpgen.generator.arginfo import (
    model_devi_args,
    run_jdata_arginfo,
    run_mdata_arginfo,
)
from dpgen.simplify.arginfo import simplify_jdata_arginfo, simplify_mdata_arginfo
from dpgen.util import normalize

init_bulk_jdata = init_bulk_jdata_arginfo()
init_surf_jdata = init_surf_jdata_arginfo()
init_reaction_jdata = init_reaction_jdata_arginfo()
simplify_jdata = simplify_jdata_arginfo()
simplify_mdata = simplify_mdata_arginfo()
run_jdata = run_jdata_arginfo()
run_mdata = run_mdata_arginfo()

# directory of examples
p_examples = Path(__file__).parent.parent / "examples"

# input_files : tuple[tuple[Argument, Path]]
#   tuple of example list
input_files = (
    (init_bulk_jdata, p_examples / "init" / "ch4.json"),
    (init_surf_jdata, p_examples / "init" / "surf.json"),
    # (init_surf_jdata, p_examples / "init" / "al.json"),
    # (init_surf_jdata, p_examples / "init" / "cu.surf.hcp.111.json"),
    (init_reaction_jdata, p_examples / "init" / "reaction.json"),
    (simplify_jdata, p_examples / "simplify" / "qm7.json"),
    (
        simplify_jdata,
        p_examples
        / "simplify-MAPbI3-scan-lebesgue"
        / "simplify_example"
        / "simplify.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-vasp" / "param_CH4_deepmd-kit-2.0.1.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-pwscf" / "param_CH4_deepmd-kit-2.0.1.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-cp2k" / "param_CH4_deepmd-kit-2.0.1.json",
    ),
    # (run_jdata, p_examples / "run" / "dp2.x-gromacs-gaussian" / "param.json"),
    (
        run_jdata,
        p_examples
        / "run"
        / "dp2.x-lammps-vasp"
        / "CH4"
        / "param_CH4_deepmd-kit-2.x.json",
    ),
    (
        run_jdata,
        p_examples
        / "run"
        / "dp2.x-lammps-vasp"
        / "Al"
        / "param_al_all_gpu-deepmd-kit-2.x.json",
    ),
    (
        run_jdata,
        p_examples
        / "run"
        / "dp2.x-lammps-pwscf"
        / "CH4"
        / "param_CH4_deepmd-kit-2.x.json",
    ),
    (
        run_jdata,
        p_examples
        / "run"
        / "dp2.x-lammps-pwscf"
        / "Al"
        / "param_al_all_gpu-deepmd-kit-2.x.json",
    ),
    (run_jdata, p_examples / "run" / "dp2.x-lammps-vasp-et" / "param_elet.json"),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-ABACUS-lcao" / "fcc-al" / "run_param.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-ABACUS-pw" / "fcc-al" / "run_param.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-cp2k" / "methane" / "param-ch4.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-ABACUS-pw" / "methane" / "param.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x-lammps-ABACUS-lcao-dpks" / "methane" / "param.json",
    ),
    (
        run_jdata,
        p_examples / "run" / "dp2.x_lammps_gaussian" / "dodecane" / "dodecane.json",
    ),
    (run_jdata, p_examples / "run" / "dp-lammps-enhance_sampling" / "param.json"),
    # (run_jdata, p_examples / "run" / "deprecated" / "param-mg-vasp.json"),
    # (run_jdata, p_examples / "run" / "deprecated" / "param-mg-vasp-ucloud.json"),
    # (run_jdata, p_examples / "run" / "deprecated" / "param-pyridine-pwscf.json"),
    # (run_jdata, p_examples / "run" / "deprecated" / "param-h2oscan-vasp.json"),
    (
        run_jdata,
        p_examples
        / "run"
        / "deprecated"
        / "dp2.x-lammps-cp2k"
        / "CH4"
        / "param_CH4.json",
    ),
    # (run_jdata, p_examples / "run" / "deprecated" / "dp2.x-lammps-pwmat" / "param_CH4.json"),
    (
        run_jdata,
        p_examples
        / "run"
        / "deprecated"
        / "dp2.x-lammps-siesta"
        / "dp-lammps-siesta"
        / "CH4"
        / "param_CH4.json",
    ),
    (
        run_jdata,
        p_examples
        / "run"
        / "deprecated"
        / "dp2.x-lammps-vasp"
        / "Al"
        / "param_al_all_gpu.json",
    ),
    (
        run_jdata,
        p_examples
        / "run"
        / "deprecated"
        / "dp2.x-lammps-vasp"
        / "CH4"
        / "param_CH4.json",
    ),
    (
        run_jdata,
        p_examples
        / "run"
        / "dp2.x-lammps-gaussian"
        / "param_C4H16N4_deepmd-kit-2.0.1.json",
    ),
    (run_jdata, p_examples / "run" / "dprc" / "generator.json"),
    (run_jdata, p_examples / "qe-cpx" / "run_param.json"),
    # machines
    (run_mdata, p_examples / "machine" / "DeePMD-kit-2.x" / "lebesgue_v2_machine.json"),
    (run_mdata, p_examples / "machine" / "DeePMD-kit-1.x" / "machine-local.json"),
    (
        run_mdata,
        p_examples / "machine" / "DeePMD-kit-1.x" / "machine-lsf-slurm-cp2k.json",
    ),
    (
        run_mdata,
        p_examples / "machine" / "DeePMD-kit-1.x" / "machine-pbs-gaussian.json",
    ),
    (run_mdata, p_examples / "machine" / "DeePMD-kit-1.x" / "machine-slurm-qe.json"),
    (run_mdata, p_examples / "machine" / "DeePMD-kit-1.0" / "machine-local-4GPU.json"),
    (run_mdata, p_examples / "CH4-refact-dpdispatcher" / "machine-ali-ehpc.json"),
    (run_mdata, p_examples / "CH4-refact-dpdispatcher" / "machine-dpcloudserver.json"),
    (
        run_mdata,
        p_examples / "run" / "dp2.x-lammps-ABACUS-lcao" / "fcc-al" / "machine.json",
    ),
    (
        run_mdata,
        p_examples / "run" / "dp2.x-lammps-ABACUS-pw" / "fcc-al" / "machine.json",
    ),
    (run_mdata, p_examples / "run" / "dp2.x-lammps-gaussian" / "machine.json"),
    (
        simplify_mdata,
        p_examples
        / "simplify-MAPbI3-scan-lebesgue"
        / "simplify_example"
        / "machine.json",
    ),
    (run_mdata, p_examples / "qe-cpx" / "run_machine.json"),
)


class TestExamples(unittest.TestCase):
    def test_arguments(self):
        for arginfo, fn in input_files:
            fn = str(fn)
            with self.subTest(fn=fn):
                with open(fn) as f:
                    data = json.load(f)
                normalize(arginfo, data)

    def test_gromacs_model_devi_arguments(self):
        """The GROMACS schema should accept every setting used by its runner."""
        arginfo = Argument("model_devi", dict, sub_variants=model_devi_args())
        data = {
            "model_devi_engine": "gromacs",
            "model_devi_jobs": [
                {
                    "sys_idx": [0],
                    "temps": [300.0],
                    "press": [],
                    "trj_freq": 10,
                    "nsteps": 100,
                    "ensemble": "nvt",
                    "lambdas": [0.5, 1.0],
                    "dt": 0.001,
                }
            ],
            "model_devi_dt": 0.002,
            "model_devi_skip": 0,
            "model_devi_f_trust_lo": 0.2,
            "model_devi_f_trust_hi": 0.6,
            "model_devi_v_trust_lo": 1e10,
            "model_devi_v_trust_hi": 1e10,
            "model_devi_clean_traj": False,
            "model_devi_nopbc": True,
            "gromacs_settings": {
                "mdp_filename": "md.mdp",
                "topol_filename": "processed.top",
                "conf_filename": "conf.gro",
                "index_filename": "index.raw",
                "type_filename": "type.raw",
                "ref_filename": "em.tpr",
                "ndx_filename": "index.ndx",
                "model_devi_script": "model_devi.py",
                "deffnm": "deepmd",
                "maxwarn": 1,
                "traj_filename": "deepmd_traj.gro",
                "group_name": "Other",
            },
        }

        normalized = arginfo.normalize_value(data)
        arginfo.check_value(normalized, strict=True)

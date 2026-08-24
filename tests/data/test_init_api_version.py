import unittest
from unittest.mock import patch

from dpgen.data import reaction, surf


class TestInitApiVersionDefaults(unittest.TestCase):
    def test_reaction_stages_default_to_supported_api(self):
        mdata = {
            "reaxff_command": "lmp",
            "reaxff_machine": {},
            "reaxff_resources": {},
            "build_command": "build",
            "build_machine": {},
            "build_resources": {"cpu_per_node": 1},
            "fp_command": "fp",
            "fp_machine": {},
            "fp_resources": {"cpu_per_node": 1},
            "fp_group_size": 1,
        }
        jdata = {
            "type_map": ["H"],
            "cutoff": 3.0,
            "dataset_size": 1,
            "qmkeywords": "force",
        }

        with (
            patch.object(reaction.glob, "glob", return_value=["task.000"]),
            patch.object(reaction, "make_submission_compat") as submit,
        ):
            reaction.run_reaxff(jdata, mdata)
            reaction.run_build_dataset(jdata, mdata)
            reaction.run_fp(jdata, mdata)

        self.assertEqual(submit.call_count, 3)
        for call in submit.call_args_list:
            self.assertEqual(call.kwargs["api_version"], "1.0")

    def test_surface_relaxation_defaults_to_supported_api(self):
        jdata = {"out_dir": "out"}
        mdata = {
            "fp_command": "vasp",
            "fp_group_size": 1,
            "fp_resources": {},
            "fp_machine": {},
        }

        with (
            patch.object(
                surf.glob,
                "glob",
                side_effect=[[], ["out/02.md/sys-0000/scale-1.000/000000"]],
            ),
            patch.object(surf, "_vasp_check_fin", return_value=False),
            patch.object(surf, "make_submission_compat") as submit,
        ):
            surf.run_vasp_relax(jdata, mdata)

        self.assertEqual(submit.call_args.kwargs["api_version"], "1.0")


if __name__ == "__main__":
    unittest.main()

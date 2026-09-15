import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dpgen.generator.arginfo import run_mdata_arginfo
from dpgen.generator.lib.run_calypso import (
    _find_models,
    _make_calypso_check_command,
    _make_calypso_opt_command,
)
from dpgen.generator.run import (
    _get_checkpoint_suffix,
    _get_export_command,
    _get_input_model_suffix,
    _get_model_suffix,
    _get_train_backend_flag,
    _normalize_training_params,
    _prepare_training_input,
    _validate_dpa_training_config,
    _validate_pt2_template_atom_map,
    make_model_devi,
    make_train_dp,
    post_train_dp,
    run_md_model_devi,
    run_model_devi,
    run_train_dp,
    training_complete_file,
)


class TestDeepmdBackendConfig(unittest.TestCase):
    def test_per_member_training_params_are_independent(self):
        jdata = {
            "numb_models": 2,
            "default_training_param": [
                {"model": {"descriptor": {"type": "dpa2"}}},
                {"model": {"descriptor": {"type": "dpa3"}}},
            ],
        }
        params = _normalize_training_params(jdata)
        self.assertEqual(params[0]["model"]["descriptor"]["type"], "dpa2")
        params[0]["model"]["descriptor"]["type"] = "changed"
        self.assertEqual(params[1]["model"]["descriptor"]["type"], "dpa3")

    def test_per_member_training_params_require_numb_models(self):
        with self.assertRaisesRegex(ValueError, "exactly numb_models"):
            _normalize_training_params(
                {
                    "numb_models": 2,
                    "default_training_param": [{"model": {}}],
                }
            )

    def test_cross_architecture_committee_compatibility(self):
        _validate_dpa_training_config(
            {
                "numb_models": 3,
                "train_backend": "pytorch",
                "model_devi_engine": "calypso",
                "type_map": ["H", "O"],
                "default_training_param": [
                    {"model": {"descriptor": {"type": "dpa2"}}},
                    {"model": {"descriptor": {"type": "dpa3"}}},
                    {"model": {"descriptor": {"type": "dpa4"}}},
                ],
            }
        )

    def test_cross_architecture_committee_rejects_cutoff_mismatch(self):
        with self.assertRaisesRegex(ValueError, "incompatible cutoff"):
            _validate_dpa_training_config(
                {
                    "numb_models": 2,
                    "type_map": ["H"],
                    "default_training_param": [
                        {
                            "model": {
                                "descriptor": {
                                    "type": "dpa2",
                                    "repinit": {"rcut": 6.0},
                                }
                            }
                        },
                        {
                            "model": {
                                "descriptor": {
                                    "type": "dpa3",
                                    "repflow": {"e_rcut": 8.0},
                                }
                            }
                        },
                    ],
                }
            )

    def test_cross_architecture_committee_rejects_parameter_mismatch(self):
        with self.assertRaisesRegex(ValueError, "incompatible fparam dimensions"):
            _validate_dpa_training_config(
                {
                    "numb_models": 2,
                    "type_map": ["H"],
                    "default_training_param": [
                        {
                            "model": {
                                "descriptor": {"type": "dpa2"},
                                "fitting_net": {
                                    "type": "ener",
                                    "numb_fparam": 1,
                                },
                            }
                        },
                        {
                            "model": {
                                "descriptor": {"type": "dpa3"},
                                "fitting_net": {
                                    "type": "ener",
                                    "numb_fparam": 2,
                                },
                            }
                        },
                    ],
                }
            )

    def test_model_dict_cutoffs_are_all_part_of_committee_signature(self):
        with self.assertRaisesRegex(ValueError, "incompatible cutoff"):
            _validate_dpa_training_config(
                {
                    "numb_models": 2,
                    "type_map": ["H"],
                    "default_training_param": [
                        {
                            "model": {
                                "model_dict": {
                                    "a": {
                                        "descriptor": {
                                            "type": "dpa2",
                                            "repinit": {"rcut": 6.0},
                                        }
                                    },
                                    "b": {
                                        "descriptor": {
                                            "type": "dpa3",
                                            "repflow": {"e_rcut": 8.0},
                                        }
                                    },
                                }
                            }
                        },
                        {
                            "model": {
                                "model_dict": {
                                    "a": {
                                        "descriptor": {
                                            "type": "dpa2",
                                            "repinit": {"rcut": 6.0},
                                        }
                                    },
                                    "b": {
                                        "descriptor": {
                                            "type": "dpa3",
                                            "repflow": {"e_rcut": 9.0},
                                        }
                                    },
                                }
                            }
                        },
                    ],
                }
            )

    def test_shared_references_are_part_of_committee_signature(self):
        def member(rcut, numb_fparam):
            return {
                "model": {
                    "shared_dict": {
                        "types": ["H"],
                        "descriptor": {
                            "type": "dpa2",
                            "repinit": {"rcut": rcut},
                        },
                        "fitting": {
                            "type": "ener",
                            "numb_fparam": numb_fparam,
                        },
                    },
                    "model_dict": {
                        "a": {
                            "type_map": "types",
                            "descriptor": "descriptor:0",
                            "fitting_net": "fitting",
                        }
                    },
                }
            }

        for changed, mismatch in (((8.0, 1), "cutoff"), ((6.0, 2), "fparam")):
            with self.subTest(mismatch=mismatch):
                with self.assertRaisesRegex(ValueError, f"incompatible {mismatch}"):
                    _validate_dpa_training_config(
                        {
                            "numb_models": 2,
                            "type_map": ["H"],
                            "default_training_param": [
                                member(6.0, 1),
                                member(*changed),
                            ],
                        }
                    )

    def test_omitted_fitting_type_defaults_to_energy(self):
        _validate_dpa_training_config(
            {
                "numb_models": 2,
                "type_map": ["H"],
                "default_training_param": [
                    {"model": {"descriptor": {"type": "se_e2_a"}, "fitting_net": {}}},
                    {
                        "model": {
                            "descriptor": {"type": "se_e2_a"},
                            "fitting_net": {"type": "ener"},
                        }
                    },
                ],
            }
        )

    def test_cross_architecture_pt2_rejects_lower_kind_mismatch(self):
        with self.assertRaisesRegex(ValueError, "same export lower kind"):
            _validate_dpa_training_config(
                {
                    "numb_models": 2,
                    "train_backend": "pt-expt",
                    "model_format": "pt2",
                    "default_training_param": [
                        {"model": {"descriptor": {"type": "dpa2"}}},
                        {"model": {"descriptor": {"type": "dpa4c"}}},
                    ],
                }
            )

    def test_legacy_defaults(self):
        cases = [
            ({}, ".pb", ".index", ""),
            ({"train_backend": "pytorch"}, ".pth", ".pt", "--pt"),
            ({"train_backend": "jax"}, ".savedmodel", ".jax", "--jax"),
        ]
        for jdata, model_suffix, checkpoint_suffix, backend_flag in cases:
            with self.subTest(jdata=jdata):
                self.assertEqual(_get_model_suffix(jdata), model_suffix)
                self.assertEqual(_get_checkpoint_suffix(jdata), checkpoint_suffix)
                self.assertEqual(_get_train_backend_flag(jdata), backend_flag)

    def test_pytorch_exportable_aliases(self):
        for backend in ("pytorch-exportable", "pt-expt"):
            with self.subTest(backend=backend):
                jdata = {"train_backend": backend}
                self.assertEqual(_get_model_suffix(jdata), ".pt2")
                self.assertEqual(_get_checkpoint_suffix(jdata), ".pt")
                self.assertEqual(_get_train_backend_flag(jdata), "--pt-expt")

    def test_explicit_pt2_formats(self):
        cases = [
            {"train_backend": "pytorch", "model_format": "pt2"},
            {"train_backend": "pytorch-exportable", "model_format": "pt2"},
            {"train_backend": "pt-expt", "model_format": "pt2"},
        ]
        for jdata in cases:
            with self.subTest(jdata=jdata):
                self.assertEqual(_get_model_suffix(jdata), ".pt2")
                self.assertEqual(_get_checkpoint_suffix(jdata), ".pt")

    def test_pte_is_rejected_for_lammps(self):
        with self.assertRaisesRegex(ValueError, "not supported by LAMMPS"):
            _get_model_suffix(
                {"train_backend": "pytorch-exportable", "model_format": "pte"}
            )
        self.assertEqual(
            _get_model_suffix(
                {"train_backend": "pt-expt", "model_devi_engine": "calypso"}
            ),
            ".pte",
        )

    def test_rejects_incompatible_model_format(self):
        with self.assertRaisesRegex(ValueError, "not available for backend"):
            _get_model_suffix({"train_backend": "tensorflow", "model_format": "pt2"})

    def test_input_models_must_share_suffix(self):
        self.assertEqual(_get_input_model_suffix(["a.pt", "b.pt"]), ".pt")
        with self.assertRaisesRegex(ValueError, "same non-empty file suffix"):
            _get_input_model_suffix(["a.pte", "b.pt2"])

    def test_pt2_template_requires_atom_map_before_read(self):
        _validate_pt2_template_atom_map(
            ["atom_modify map yes\n", "read_data conf.lmp\n"]
        )
        for atom_map in ("array", "hash"):
            with self.subTest(atom_map=atom_map):
                _validate_pt2_template_atom_map(
                    [f"atom_modify map {atom_map}\n", "read_data conf.lmp\n"]
                )
        for lines in (
            ["read_data conf.lmp\n"],
            ["read_restart restart.100\n", "atom_modify map yes\n"],
            [
                'if "${restart} > 0" then "read_restart restart.100" '
                'else "read_data conf.lmp"\n',
                "atom_modify map yes\n",
            ],
        ):
            with self.subTest(lines=lines):
                with self.assertRaisesRegex(ValueError, "atom_modify map"):
                    _validate_pt2_template_atom_map(lines)
        with self.assertRaisesRegex(ValueError, "read_data or read_restart"):
            _validate_pt2_template_atom_map(["atom_modify map yes\n"])

    def test_dpa_backend_and_compile_option_validation(self):
        _validate_dpa_training_config(
            {
                "train_backend": "pytorch",
                "model_format": "pt2",
                "default_training_param": {
                    "model": {
                        "type": "DPA4",
                        "use_compile": True,
                        "enable_tf32": True,
                    },
                    "training": {},
                },
            }
        )
        _validate_dpa_training_config(
            {
                "train_backend": "pytorch",
                "model_format": "pt2",
                "default_training_param": {"model": {"type": "SeZM"}},
            }
        )
        _validate_dpa_training_config(
            {
                "train_backend": "pt-expt",
                "model_format": "pt2",
                "default_training_param": {
                    "model": {"descriptor": {"type": "DPA4C"}},
                    "training": {
                        "enable_compile": True,
                        "enable_tf32": True,
                    },
                },
            }
        )
        with self.assertRaisesRegex(ValueError, "requires train_backend='pytorch'"):
            _validate_dpa_training_config(
                {
                    "train_backend": "pytorch-exportable",
                    "model_format": "pt2",
                    "default_training_param": {
                        "model": {"descriptor": {"type": "dpa4"}},
                        "training": {},
                    },
                }
            )
        with self.assertRaisesRegex(ValueError, "training.enable_compile"):
            _validate_dpa_training_config(
                {
                    "train_backend": "pytorch-exportable",
                    "model_format": "pt2",
                    "default_training_param": {
                        "model": {
                            "descriptor": {"type": "dpa4c"},
                            "use_compile": True,
                        },
                        "training": {},
                    },
                }
            )
        with self.assertRaisesRegex(ValueError, "cannot mix DPA4 and DPA4C"):
            _validate_dpa_training_config(
                {
                    "train_backend": "pytorch-exportable",
                    "model_format": "pt2",
                    "default_training_param": {
                        "model": {
                            "model_dict": {
                                "dpa4": {"descriptor": {"type": "dpa4"}},
                                "dpa4c": {"descriptor": {"type": "dpa4c"}},
                            }
                        },
                        "training": {},
                    },
                }
            )
        with self.assertRaisesRegex(ValueError, "only exports pt2 for DPA4/SeZM"):
            _validate_dpa_training_config(
                {"train_backend": "pytorch", "model_format": "pt2"}
            )

    def test_calypso_discovers_resolved_model_suffix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir)
            for name in ("graph.000.pb", "graph.000.pte", "graph.001.pte"):
                (path / name).touch()
            self.assertEqual(len(_find_models(path, ".pb")), 1)
            self.assertEqual(len(_find_models(path, ".pte")), 2)

    def test_calypso_optimizer_uses_resolved_model(self):
        command = _make_calypso_opt_command("python", "graph.000.pt2")
        self.assertIn("calypso_run_opt.py --model ../graph.000.pt2", command)

    def test_calypso_recovery_uses_resolved_model(self):
        command = _make_calypso_check_command("python", "graph.000.pt2")
        self.assertEqual(command, "python check_outcar.py --model ../graph.000.pt2")

    def test_export_command_uses_deployment_executable(self):
        self.assertEqual(
            _get_export_command(
                {"train_command": "/train/dp", "train_export_command": "/deploy/dp"},
                "--pt",
            ),
            "/deploy/dp --pt",
        )

    def test_export_command_is_accepted_in_machine_configuration(self):
        train = run_mdata_arginfo().sub_fields["train"]
        self.assertIn("export_command", train.sub_fields)


class TestRunTrainDeepmdBackend(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.old_cwd = os.getcwd()
        os.chdir(self.tmp.name)
        self.addCleanup(os.chdir, self.old_cwd)
        self.mdata = {
            "api_version": "1.0",
            "deepmd_version": "3.2.0",
            "train_command": "dp",
            "train_machine": {"name": "train"},
            "train_resources": {"queue": "train"},
            "model_devi_machine": {"name": "model-devi"},
            "model_devi_resources": {"queue": "model-devi"},
        }

    def _run(self, **updates):
        jdata = {"numb_models": 1, "one_h5": True}
        jdata.update(updates)
        with patch("dpgen.generator.run.make_submission") as make_submission:
            run_train_dp(0, jdata, self.mdata)
        calls = []
        for call in make_submission.call_args_list:
            details = dict(call.kwargs)
            details["machine"], details["resources"] = call.args[:2]
            calls.append(details)
        return calls[0] if len(calls) == 1 else calls

    def test_legacy_tensorflow_commands_are_preserved(self):
        call = self._run()
        self.assertIn("model.ckpt.index", call["commands"][0])
        self.assertEqual(call["commands"][1], "dp freeze")
        self.assertIn("frozen_model.pb", call["backward_files"])
        self.assertIn("model.ckpt.index", call["backward_files"])

    def test_prepare_training_input_preserves_member_specific_model(self):
        first = {"model": {"descriptor": {"type": "dpa2"}}, "training": {}}
        second = {"model": {"descriptor": {"type": "dpa3"}}, "training": {}}
        for item in (first, second):
            _prepare_training_input(
                item,
                "3.2.0",
                ["system"],
                [1],
                ["H"],
                0,
                None,
                0,
                None,
                "auto",
                None,
                None,
                None,
                None,
            )
        self.assertEqual(first["model"]["descriptor"]["type"], "dpa2")
        self.assertEqual(second["model"]["descriptor"]["type"], "dpa3")
        self.assertEqual(first["model"]["type_map"], ["H"])
        self.assertEqual(second["model"]["type_map"], ["H"])

    def test_prepare_training_input_sets_ele_temp_for_model_dict(self):
        for mode, enabled, removed in (
            (1, "numb_fparam", "numb_aparam"),
            (2, "numb_aparam", "numb_fparam"),
        ):
            item = {
                "model": {
                    "model_dict": {
                        "a": {"fitting_net": {}},
                        "b": {"fitting_net": {}},
                    }
                },
                "training": {},
            }
            _prepare_training_input(
                item,
                "3.2.0",
                ["system"],
                [1],
                ["H"],
                mode,
                None,
                0,
                None,
                "auto",
                None,
                None,
                None,
                None,
            )
            for branch in item["model"]["model_dict"].values():
                self.assertEqual(branch["fitting_net"][enabled], 1)
                self.assertNotIn(removed, branch["fitting_net"])

    def test_make_train_generates_cross_architecture_inputs(self):
        jdata = {
            "numb_models": 3,
            "init_data_prefix": "data",
            "init_data_sys": [],
            "sys_configs": [],
            "model_devi_jobs": [{}],
            "fp_task_min": 0,
            "type_map": ["H"],
            "default_training_param": [
                {"model": {"descriptor": {"type": "dpa2"}}, "training": {}},
                {"model": {"descriptor": {"type": "dpa3"}}, "training": {}},
                {"model": {"descriptor": {"type": "dpa4"}}, "training": {}},
            ],
        }
        mdata = {"deepmd_version": "3.2.0"}
        with patch("dpgen.generator.run.os.symlink"):
            make_train_dp(0, jdata, mdata)
        descriptors = []
        for index in range(3):
            with open(
                Path("iter.000000") / "00.train" / f"{index:03d}" / "input.json"
            ) as fp:
                descriptors.append(json.load(fp)["model"]["descriptor"]["type"])
        self.assertEqual(descriptors, ["dpa2", "dpa3", "dpa4"])

    def test_make_train_seeds_model_dict_branches(self):
        jdata = {
            "numb_models": 1,
            "init_data_prefix": "data",
            "init_data_sys": [],
            "sys_configs": [],
            "model_devi_jobs": [{}],
            "fp_task_min": 0,
            "type_map": ["H"],
            "default_training_param": {
                "model": {
                    "model_dict": {
                        "a": {
                            "descriptor": {"type": "dpa2"},
                            "fitting_net": {},
                        },
                        "b": {
                            "descriptor": {"type": "dpa3"},
                            "fitting_net": {},
                        },
                    }
                },
                "training": {},
            },
        }
        with patch("dpgen.generator.run.os.symlink"):
            make_train_dp(0, jdata, {"deepmd_version": "3.2.0"})
        with open(Path("iter.000000") / "00.train" / "000" / "input.json") as fp:
            model_dict = json.load(fp)["model"]["model_dict"]
        for branch in model_dict.values():
            self.assertIn("seed", branch["descriptor"])
            self.assertIn("seed", branch["fitting_net"])

    def test_make_train_seeds_shared_descriptor_branch_heads(self):
        jdata = {
            "numb_models": 2,
            "init_data_prefix": "data",
            "init_data_sys": [],
            "sys_configs": [],
            "model_devi_jobs": [{}],
            "fp_task_min": 0,
            "type_map": ["H"],
            "default_training_param": {
                "model": {
                    "shared_dict": {
                        "descriptor": {"type": "dpa2", "seed": 1},
                        "fitting": {"type": "ener", "seed": 1},
                    },
                    "model_dict": {
                        "a": {
                            "descriptor": "descriptor:0",
                            "fitting_net": "fitting",
                            "type_embedding": {},
                        }
                    },
                },
                "training": {},
            },
        }
        with patch("dpgen.generator.run.os.symlink"):
            make_train_dp(0, jdata, {"deepmd_version": "3.2.0"})
        shared_seeds = []
        for index in range(2):
            with open(
                Path("iter.000000") / "00.train" / f"{index:03d}" / "input.json"
            ) as fp:
                model = json.load(fp)["model"]
            branch = model["model_dict"]["a"]
            self.assertEqual(branch["descriptor"], "descriptor:0")
            self.assertEqual(branch["fitting_net"], "fitting")
            self.assertIn("seed", branch["type_embedding"])
            shared_seeds.append(
                (
                    model["shared_dict"]["descriptor"]["seed"],
                    model["shared_dict"]["fitting"]["seed"],
                )
            )
        self.assertNotEqual(shared_seeds[0], shared_seeds[1])

    def test_cross_architecture_submission_tracks_all_members(self):
        calls = self._run(
            numb_models=3,
            train_backend="pytorch",
            model_devi_engine="calypso",
            default_training_param=[
                {"model": {"descriptor": {"type": "dpa2"}}},
                {"model": {"descriptor": {"type": "dpa3"}}},
                {"model": {"descriptor": {"type": "dpa4"}}},
            ],
        )
        self.assertEqual(calls["run_tasks"], ["000", "001", "002"])

    def test_pytorch_dpa4_exports_pt2_on_model_devi_resources(self):
        train_call, export_call = self._run(
            train_backend="pytorch",
            model_format="pt2",
            default_training_param={"model": {"type": "dpa4"}},
        )
        self.assertEqual(train_call["machine"], self.mdata["train_machine"])
        self.assertEqual(len(train_call["commands"]), 1)
        self.assertIn("dp --pt train", train_call["commands"][0])
        self.assertIn(f"touch {training_complete_file}", train_call["commands"][0])
        self.assertIn(training_complete_file, train_call["backward_files"])
        self.assertIn("model.ckpt.pt", train_call["backward_files"])
        self.assertNotIn("frozen_model.pt2", train_call["backward_files"])
        self.assertEqual(export_call["machine"], self.mdata["model_devi_machine"])
        self.assertEqual(export_call["resources"], self.mdata["model_devi_resources"])
        self.assertEqual(
            export_call["commands"],
            ["dp --pt freeze -c model.ckpt.pt -o frozen_model"],
        )
        self.assertEqual(export_call["forward_files"], ["model.ckpt.pt"])
        self.assertIn("frozen_model.pt2", export_call["backward_files"])

    def test_pytorch_exportable_nlist_pt2_export(self):
        _, export_call = self._run(
            train_backend="pt-expt",
            model_format="pt2",
            default_training_param={"model": {"descriptor": {"type": "se_e2_a"}}},
        )
        self.assertEqual(
            export_call["commands"],
            [
                "dp --pt-expt freeze -c model.ckpt.pt "
                "-o frozen_model.pt2 --lower-kind nlist"
            ],
        )

    def test_pt2_export_uses_deployment_command(self):
        self.mdata["train_command"] = "/train/dp"
        self.mdata["train_export_command"] = "/deploy/dp"
        _, export_call = self._run(
            train_backend="pytorch",
            model_format="pt2",
            default_training_param={"model": {"type": "dpa4"}},
        )
        self.assertEqual(
            export_call["commands"],
            ["/deploy/dp --pt freeze -c model.ckpt.pt -o frozen_model"],
        )

    def test_pt2_export_reuses_completed_training_checkpoints(self):
        task = Path("iter.000000") / "00.train" / "000"
        task.mkdir(parents=True)
        (task / "model.ckpt.pt").touch()
        (task / training_complete_file).touch()
        with patch("dpgen.generator.run.make_submission") as make_submission:
            run_train_dp(
                0,
                {
                    "numb_models": 1,
                    "one_h5": True,
                    "train_backend": "pytorch",
                    "model_format": "pt2",
                    "default_training_param": {"model": {"type": "dpa4"}},
                },
                self.mdata,
            )
        self.assertEqual(make_submission.call_count, 1)
        self.assertEqual(
            make_submission.call_args.kwargs["commands"],
            ["dp --pt freeze -c model.ckpt.pt -o frozen_model"],
        )

    def test_pt2_export_does_not_skip_incomplete_training(self):
        task = Path("iter.000000") / "00.train" / "000"
        task.mkdir(parents=True)
        (task / "model.ckpt.pt").touch()
        with patch("dpgen.generator.run.make_submission") as make_submission:
            run_train_dp(
                0,
                {
                    "numb_models": 1,
                    "one_h5": True,
                    "train_backend": "pytorch",
                    "model_format": "pt2",
                    "default_training_param": {"model": {"type": "dpa4"}},
                },
                self.mdata,
            )
        self.assertEqual(make_submission.call_count, 2)

    def test_legacy_pytorch_commands_are_preserved(self):
        call = self._run(train_backend="pytorch")
        self.assertIn("dp --pt train", call["commands"][0])
        self.assertIn("model.ckpt.pt", call["commands"][0])
        self.assertEqual(call["commands"][1], "dp --pt freeze")
        self.assertIn("frozen_model.pth", call["backward_files"])
        self.assertIn("model.ckpt.pt", call["backward_files"])

    def test_pytorch_exportable_dense_pte(self):
        call = self._run(
            train_backend="pytorch-exportable", model_devi_engine="calypso"
        )
        self.assertIn("dp --pt-expt train", call["commands"][0])
        self.assertEqual(
            call["commands"][1],
            "dp --pt-expt freeze -o frozen_model.pte",
        )
        self.assertIn("frozen_model.pte", call["backward_files"])

    def test_pytorch_exportable_dpa4c_pt2_compression(self):
        train_call, export_call = self._run(
            train_backend="pt-expt",
            model_format="pt2",
            dp_compress=True,
            default_training_param={"model": {"descriptor": {"type": "dpa4c"}}},
        )
        self.assertIn("dp --pt-expt train", train_call["commands"][0])
        self.assertEqual(len(train_call["commands"]), 1)
        self.assertEqual(
            export_call["commands"],
            [
                "dp --pt-expt freeze -c model.ckpt.pt "
                "-o frozen_model.pt2 --lower-kind graph",
                "dp --pt-expt compress -i frozen_model.pt2 -o frozen_model_compressed.pt2",
            ],
        )
        self.assertEqual(export_call["forward_files"], ["model.ckpt.pt"])
        self.assertIn("frozen_model_compressed.pt2", export_call["backward_files"])

    def test_multiple_pt2_models_are_exported(self):
        train_call, export_call = self._run(
            numb_models=4,
            train_backend="pt-expt",
            model_format="pt2",
        )
        expected_tasks = [f"{index:03d}" for index in range(4)]
        self.assertEqual(train_call["run_tasks"], expected_tasks)
        self.assertEqual(export_call["run_tasks"], expected_tasks)

    def test_regular_pytorch_pt2_compression_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "cannot compress pt2"):
            self._run(
                train_backend="pytorch",
                model_format="pt2",
                dp_compress=True,
                default_training_param={"model": {"type": "dpa4"}},
            )

    def test_exportable_init_frozen_model_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "does not support"):
            self._run(
                train_backend="pt-expt",
                training_init_frozen_model=["model.pt2"],
            )

    def test_deepmd_31_is_rejected_for_pt2(self):
        self.mdata["deepmd_version"] = "3.1.0"
        with self.assertRaisesRegex(RuntimeError, "3.2 or later"):
            self._run(
                train_backend="pytorch",
                model_format="pt2",
                default_training_param={"model": {"type": "dpa4"}},
            )

    def test_finetune_keeps_source_model_suffix(self):
        train_call, _ = self._run(
            train_backend="pt-expt",
            model_format="pt2",
            training_finetune_model=["source.pt"],
        )
        self.assertIn("--finetune old/init.pt", train_call["commands"][0])
        self.assertIn(str(Path("old") / "init.pt"), train_call["forward_files"])

    def test_post_train_links_pt2_model(self):
        jdata = {
            "numb_models": 1,
            "train_backend": "pt-expt",
            "model_format": "pt2",
        }
        with patch("dpgen.generator.run.os.symlink") as symlink:
            post_train_dp(0, jdata, self.mdata)
        symlink.assert_called_once_with(
            str(Path("000") / "frozen_model.pt2"),
            str(Path("iter.000000") / "00.train" / "graph.000.pt2"),
        )

    def test_model_deviation_forwards_pt2_models(self):
        work_path = Path("iter.000000") / "01.model_devi"
        (work_path / "task.000.000000").mkdir(parents=True)
        (work_path / "graph.000.pt2").touch()
        (work_path / "cur_job.json").write_text(json.dumps({}), encoding="utf-8")
        jdata = {
            "train_backend": "pytorch",
            "model_format": "pt2",
            "model_devi_jobs": [{}],
        }
        mdata = {
            "api_version": "1.0",
            "model_devi_command": "lmp -k on g 1 -sf kk",
            "model_devi_group_size": 1,
            "model_devi_machine": {},
            "model_devi_resources": {},
        }
        with patch("dpgen.generator.run.make_submission") as make_submission:
            run_md_model_devi(0, jdata, mdata)
        call = make_submission.call_args.kwargs
        self.assertEqual(call["forward_common_files"], ["graph.000.pt2"])
        self.assertIn("lmp -k on g 1 -sf kk", call["commands"][0])

    def test_committee_training_to_model_deviation_handoff(self):
        train_jdata = {
            "numb_models": 3,
            "train_backend": "pytorch",
            "model_devi_engine": "calypso",
            "type_map": ["H", "O"],
            "model_devi_jobs": [{"times": [0], "PSTRESS": [0.0]}],
            "sys_configs": [],
            "default_training_param": [
                {
                    "model": {
                        "type_map": ["O", "H"],
                        "descriptor": {"type": "dpa2"},
                    }
                },
                {
                    "model": {
                        "type_map": ["O", "H"],
                        "descriptor": {"type": "dpa3"},
                    }
                },
                {
                    "model": {
                        "type_map": ["O", "H"],
                        "descriptor": {"type": "dpa4"},
                    }
                },
            ],
        }
        train_call = self._run(**train_jdata)
        self.assertEqual(train_call["run_tasks"], ["000", "001", "002"])

        train_path = Path("iter.000000") / "00.train"
        for index in range(3):
            task_path = train_path / f"{index:03d}"
            task_path.mkdir(parents=True)
            (task_path / "frozen_model.pth").write_text(
                f"committee member {index}", encoding="utf-8"
            )

        def copy_link(source, destination):
            destination = Path(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(destination.parent / source, destination)

        with patch("dpgen.generator.run.os.symlink", side_effect=copy_link):
            post_train_dp(0, train_jdata, self.mdata)
        train_models = sorted(train_path.glob("graph.*.pth"))
        self.assertEqual(
            [model.name for model in train_models],
            ["graph.000.pth", "graph.001.pth", "graph.002.pth"],
        )
        self.assertEqual(
            [model.read_text(encoding="utf-8") for model in train_models],
            ["committee member 0", "committee member 1", "committee member 2"],
        )

        with patch("dpgen.generator.run._make_model_devi_native_calypso"):
            with patch("dpgen.generator.run.os.symlink", side_effect=shutil.copyfile):
                self.assertTrue(make_model_devi(0, train_jdata, self.mdata))

        calypso_path = Path("iter.000000") / "01.model_devi"
        staged_models = sorted(
            (calypso_path / "gen_stru_analy.000").glob("graph.*.pth")
        )
        self.assertEqual(
            [model.name for model in staged_models],
            ["graph.000.pth", "graph.001.pth", "graph.002.pth"],
        )
        self.assertEqual(
            [model.read_text(encoding="utf-8") for model in staged_models],
            ["committee member 0", "committee member 1", "committee member 2"],
        )

        expected_models = " ".join(str(model.resolve()) for model in staged_models)
        record_path = (calypso_path / "record.calypso").resolve()
        record_path.write_text("3\n", encoding="utf-8")
        commands = []

        def run_calypso(command):
            commands.append(command)
            record_path.write_text("3\n4\n", encoding="utf-8")
            return 0

        with patch(
            "dpgen.generator.lib.run_calypso.os.system", side_effect=run_calypso
        ):
            run_model_devi(0, train_jdata, {"model_devi_deepmdkit_python": "python"})

        self.assertEqual(len(commands), 1)
        self.assertIn(
            f"--all_models {expected_models} --type_map H O --model_type_map O H",
            commands[0],
        )

    def test_gromacs_model_deviation_forwards_configured_script(self):
        work_path = Path("iter.000000") / "01.model_devi"
        (work_path / "task.000.000000").mkdir(parents=True)
        (work_path / "graph.000.pb").touch()
        (work_path / "cur_job.json").write_text(json.dumps({}), encoding="utf-8")
        jdata = {
            "model_devi_engine": "gromacs",
            "gromacs_settings": {"model_devi_script": "model_devi.py"},
            "model_devi_jobs": [{}],
        }
        mdata = {
            "api_version": "1.0",
            "model_devi_command": "gmx",
            "model_devi_group_size": 1,
            "model_devi_machine": {},
            "model_devi_resources": {},
        }
        with patch("dpgen.generator.run.make_submission") as make_submission:
            run_md_model_devi(0, jdata, mdata)

        forward_files = make_submission.call_args.kwargs["forward_files"]
        self.assertIn("model_devi.py", forward_files)


if __name__ == "__main__":
    unittest.main()

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dpgen.generator.lib.run_calypso import (
    _find_models,
    _make_calypso_opt_command,
)
from dpgen.generator.run import (
    _get_checkpoint_suffix,
    _get_input_model_suffix,
    _get_model_suffix,
    _get_train_backend_flag,
    _validate_dpa_training_config,
    _validate_pt2_template_atom_map,
    post_train_dp,
    run_md_model_devi,
    run_train_dp,
)


class TestDeepmdBackendConfig(unittest.TestCase):
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
                with self.assertRaisesRegex(ValueError, "atom_modify map yes"):
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

    def test_pytorch_dpa4_exports_pt2_on_model_devi_resources(self):
        train_call, export_call = self._run(
            train_backend="pytorch",
            model_format="pt2",
            default_training_param={"model": {"type": "dpa4"}},
        )
        self.assertEqual(train_call["machine"], self.mdata["train_machine"])
        self.assertEqual(len(train_call["commands"]), 1)
        self.assertIn("dp --pt train", train_call["commands"][0])
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
        )
        self.assertIn("dp --pt-expt train", train_call["commands"][0])
        self.assertEqual(len(train_call["commands"]), 1)
        self.assertEqual(
            export_call["commands"],
            [
                "dp --pt-expt freeze -c model.ckpt.pt -o frozen_model --lower-kind graph",
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


if __name__ == "__main__":
    unittest.main()

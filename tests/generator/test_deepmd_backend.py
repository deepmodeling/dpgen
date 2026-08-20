import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dpgen.generator.run import (
    _get_checkpoint_suffix,
    _get_input_model_suffix,
    _get_model_backend_flag,
    _get_model_suffix,
    _get_train_backend_flag,
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
                self.assertEqual(_get_model_suffix(jdata), ".pte")
                self.assertEqual(_get_checkpoint_suffix(jdata), ".pt")
                self.assertEqual(_get_train_backend_flag(jdata), "--pt-expt")

    def test_explicit_pt2_formats(self):
        cases = [
            {"train_backend": "pytorch", "model_format": "pt2"},
            {"train_backend": "pytorch-exportable", "model_format": "pt2"},
            {"train_backend": "pt-expt", "model_format": "pt2"},
            {
                "train_backend": "pytorch",
                "model_devi_backend": "pytorch-exportable",
                "model_format": "pt2",
            },
        ]
        for jdata in cases:
            with self.subTest(jdata=jdata):
                self.assertEqual(_get_model_suffix(jdata), ".pt2")
                self.assertEqual(_get_checkpoint_suffix(jdata), ".pt")

    def test_pytorch_checkpoint_can_use_exportable_deployment(self):
        jdata = {
            "train_backend": "pytorch",
            "model_devi_backend": "pt-expt",
            "model_format": "pt2",
        }
        self.assertEqual(_get_train_backend_flag(jdata), "--pt")
        self.assertEqual(_get_model_backend_flag(jdata), "--pt-expt")

    def test_rejects_other_cross_backend_exports(self):
        with self.assertRaisesRegex(ValueError, "Cannot export models"):
            _get_model_suffix(
                {
                    "train_backend": "tensorflow",
                    "model_devi_backend": "pytorch-exportable",
                }
            )

    def test_rejects_incompatible_model_format(self):
        with self.assertRaisesRegex(ValueError, "not available for backend"):
            _get_model_suffix({"train_backend": "tensorflow", "model_format": "pt2"})

    def test_input_models_must_share_suffix(self):
        self.assertEqual(_get_input_model_suffix(["a.pt", "b.pt"]), ".pt")
        with self.assertRaisesRegex(ValueError, "same non-empty file suffix"):
            _get_input_model_suffix(["a.pte", "b.pt2"])


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
            "train_machine": {},
            "train_resources": {},
        }

    def _run(self, **updates):
        jdata = {"numb_models": 1, "one_h5": True}
        jdata.update(updates)
        with patch("dpgen.generator.run.make_submission") as make_submission:
            run_train_dp(0, jdata, self.mdata)
        return make_submission.call_args.kwargs

    def test_legacy_tensorflow_commands_are_preserved(self):
        call = self._run()
        self.assertIn("model.ckpt.index", call["commands"][0])
        self.assertEqual(call["commands"][1], "dp freeze")
        self.assertIn("frozen_model.pb", call["backward_files"])
        self.assertIn("model.ckpt.index", call["backward_files"])

    def test_pytorch_dpa4_uses_exportable_pt2_deployment(self):
        call = self._run(
            train_backend="pytorch",
            model_devi_backend="pytorch-exportable",
            model_format="pt2",
        )
        self.assertIn("dp --pt train", call["commands"][0])
        self.assertIn("model.ckpt.pt", call["commands"][0])
        self.assertEqual(
            call["commands"][1],
            "dp --pt-expt freeze -o frozen_model.pt2 --lower-kind graph",
        )
        self.assertIn("frozen_model.pt2", call["backward_files"])
        self.assertIn("model.ckpt.pt", call["backward_files"])

    def test_legacy_pytorch_commands_are_preserved(self):
        call = self._run(train_backend="pytorch")
        self.assertIn("dp --pt train", call["commands"][0])
        self.assertIn("model.ckpt.pt", call["commands"][0])
        self.assertEqual(call["commands"][1], "dp --pt freeze")
        self.assertIn("frozen_model.pth", call["backward_files"])
        self.assertIn("model.ckpt.pt", call["backward_files"])

    def test_pytorch_exportable_dense_pte(self):
        call = self._run(train_backend="pytorch-exportable")
        self.assertIn("dp --pt-expt train", call["commands"][0])
        self.assertEqual(
            call["commands"][1],
            "dp --pt-expt freeze -o frozen_model.pte",
        )
        self.assertIn("frozen_model.pte", call["backward_files"])

    def test_pytorch_exportable_dpa4c_pt2_compression(self):
        call = self._run(train_backend="pt-expt", model_format="pt2", dp_compress=True)
        self.assertEqual(
            call["commands"][1],
            "dp --pt-expt freeze -o frozen_model.pt2 --lower-kind graph",
        )
        self.assertEqual(
            call["commands"][2],
            "dp --pt-expt compress -i frozen_model.pt2 -o frozen_model_compressed.pt2",
        )
        self.assertIn("frozen_model_compressed.pt2", call["backward_files"])

    def test_regular_pytorch_pt2_compression_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "cannot compress pt2"):
            self._run(train_backend="pytorch", model_format="pt2", dp_compress=True)

    def test_exportable_init_frozen_model_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "does not support"):
            self._run(
                train_backend="pt-expt",
                training_init_frozen_model=["model.pt2"],
            )

    def test_deepmd_31_is_rejected_for_pt2(self):
        self.mdata["deepmd_version"] = "3.1.0"
        with self.assertRaisesRegex(RuntimeError, "3.2 or later"):
            self._run(train_backend="pytorch", model_format="pt2")

    def test_finetune_keeps_source_model_suffix(self):
        call = self._run(
            train_backend="pt-expt",
            model_format="pt2",
            training_finetune_model=["source.pt"],
        )
        self.assertIn("--finetune old/init.pt", call["commands"][0])
        self.assertIn(str(Path("old") / "init.pt"), call["forward_files"])

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
            "model_devi_backend": "pytorch-exportable",
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

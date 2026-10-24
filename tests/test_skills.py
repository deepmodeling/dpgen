import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from dpgen.main import main, main_parser
from dpgen.skills import _read, _root


class TestSkills(unittest.TestCase):
    def invoke(self, argv):
        output = io.StringIO()
        with patch("sys.argv", ["dpgen", *argv]), contextlib.redirect_stdout(output):
            main()
        return output.getvalue()

    def test_list_and_read(self):
        names = self.invoke(["skills", "list"]).splitlines()
        self.assertEqual(names, ["dpgen-run", "dpgen-simplify"])
        for name in names:
            self.assertEqual(
                self.invoke(["skills", "read", name]),
                _root().joinpath(name, "SKILL.md").read_bytes().decode("utf-8"),
            )

    def test_references_and_templates(self):
        for name, path in [
            ("dpgen-run", "references/monitoring.md"),
            ("dpgen-simplify", "assets/param.template.json"),
        ]:
            self.assertEqual(
                self.invoke(["skills", "read", name, path]),
                _root().joinpath(name, path).read_bytes().decode("utf-8"),
            )

    def test_no_workflow_or_info(self):
        with (
            patch("dpgen.main.info") as info,
            patch("dpgen.main._run_command") as imports,
        ):
            self.invoke(["skills", "list"])
        info.assert_not_called()
        imports.assert_not_called()

    def test_bad_paths(self):
        for path in [
            "/etc/passwd",
            "../SKILL.md",
            "references/../../SKILL.md",
            "a\\b",
            "C:/a",
            "a\x00b",
            "./SKILL.md",
            "missing.md",
        ]:
            with self.subTest(path=path), self.assertRaises(ValueError):
                _read(_root(), "dpgen-run", path)
        for name in ["../dpgen-run", "/dpgen-run", "missing"]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                _read(_root(), name, "SKILL.md")

    def test_symlink_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "skills"
            skill = root / "example"
            skill.mkdir(parents=True)
            (skill / "SKILL.md").write_text("example", encoding="utf-8")
            secret = Path(directory) / "secret"
            secret.write_text("outside", encoding="utf-8")
            (skill / "escape.md").symlink_to(secret)
            with self.assertRaises(ValueError):
                _read(root, "example", "escape.md")

    def test_unknown_is_usage_error(self):
        with (
            contextlib.redirect_stderr(io.StringIO()),
            self.assertRaises(SystemExit) as error,
        ):
            self.invoke(["skills", "read", "missing"])
        self.assertEqual(error.exception.code, 2)

    def test_help(self):
        self.assertIn("skills", main_parser().format_help())
        for command in [
            [],
            ["run"],
            ["simplify"],
            ["init_surf"],
            ["init_bulk"],
            ["init_reaction"],
            ["autotest"],
            ["skills"],
            ["skills", "read"],
        ]:
            with self.subTest(command=command), self.assertRaises(SystemExit) as error:
                self.invoke([*command, "--help"])
            self.assertEqual(error.exception.code, 0)

    def test_lazy_dispatch_mapping(self):
        commands = [
            (["init_surf", "p"], "dpgen.data.surf", "gen_init_surf"),
            (["init_bulk", "p"], "dpgen.data.gen", "gen_init_bulk"),
            (["init_reaction", "p"], "dpgen.data.reaction", "gen_init_reaction"),
            (["auto_gen_param", "p"], "dpgen.tools.auto_gen_param", "auto_gen_param"),
            (["run", "p", "m"], "dpgen.generator.run", "gen_run"),
            (["run/report", "job"], "dpgen.tools.run_report", "run_report"),
            (["collect", "job", "out"], "dpgen.collect.collect", "gen_collect"),
            (["simplify", "p", "m"], "dpgen.simplify.simplify", "gen_simplify"),
            (["autotest", "make", "p"], "dpgen.auto_test.run", "gen_test"),
            (["db", "p"], "dpgen.database.run", "db_run"),
            (["gui"], "dpgen.gui", "start_dpgui"),
        ]
        for argv, module, function in commands:
            with self.subTest(argv=argv):
                args = main_parser().parse_args(argv)
                implementation = Mock()
                with patch(
                    "dpgen.main.importlib.import_module", return_value=implementation
                ) as imports:
                    args.func(args)
                imports.assert_called_once_with(module)
                getattr(implementation, function).assert_called_once_with(args)

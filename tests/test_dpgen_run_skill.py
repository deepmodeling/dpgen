import copy
import json
import re
import unittest
from pathlib import Path

from dargs.dargs import (
    ArgumentKeyError,
    ArgumentTypeError,
    ArgumentValueError,
)

from dpgen.generator.arginfo import run_jdata_arginfo
from dpgen.remote.decide_machine import convert_mdata
from dpgen.util import normalize

LINK_PATTERN = re.compile(r"\[[^]]+\]\(([^)]+)\)")
SKILL_RELATIVE_PATH = Path("skills") / "dpgen-run" / "SKILL.md"


def skill_markdown_files(repository_root):
    skill_root = repository_root / SKILL_RELATIVE_PATH.parent
    return sorted(skill_root.rglob("*.md"))


def linked_json_examples(repository_root, examples_root):
    examples = set()
    for markdown_path in skill_markdown_files(repository_root):
        for target in LINK_PATTERN.findall(markdown_path.read_text()):
            target_path = target.split("#", 1)[0]
            if not target_path.endswith(".json"):
                continue
            resolved_target = (markdown_path.parent / target_path).resolve()
            try:
                resolved_target.relative_to(examples_root)
            except ValueError:
                continue
            examples.add(resolved_target)
    return examples


class TestDPGenRunSkill(unittest.TestCase):
    def test_main_skill_is_small_progressive_router(self):
        repository_root = Path(__file__).resolve().parents[1]
        skill_path = repository_root / SKILL_RELATIVE_PATH
        skill_text = skill_path.read_text()

        self.assertLessEqual(len(skill_text.splitlines()), 50)
        reference_targets = {
            target
            for target in LINK_PATTERN.findall(skill_text)
            if target.startswith("references/")
        }
        self.assertEqual(len(reference_targets), 4)

    def test_all_repository_relative_links_exist(self):
        repository_root = Path(__file__).resolve().parents[1]
        missing = []
        for markdown_path in skill_markdown_files(repository_root):
            for target in LINK_PATTERN.findall(markdown_path.read_text()):
                if "://" in target or target.startswith("#"):
                    continue
                resolved_target = (
                    markdown_path.parent / target.split("#", 1)[0]
                ).resolve()
                try:
                    resolved_target.relative_to(repository_root)
                except ValueError:
                    missing.append(f"{markdown_path.name}: {target}")
                else:
                    if not resolved_target.exists():
                        missing.append(f"{markdown_path.name}: {target}")

        self.assertEqual(
            missing,
            [],
            "Missing repository-relative links in dpgen-run skill: "
            + ", ".join(missing),
        )

    def test_linked_parameter_examples_use_current_schema(self):
        repository_root = Path(__file__).resolve().parents[1]
        examples_root = repository_root / "examples" / "run"
        examples = linked_json_examples(repository_root, examples_root)
        self.assertTrue(examples, "No examples/run/*.json links found")

        for example in examples:
            with self.subTest(example=example):
                parameter_data = json.loads(example.read_text())
                normalized = normalize(
                    run_jdata_arginfo(),
                    parameter_data,
                    strict_check=False,
                )
                self.assertIsInstance(normalized, dict)
                self.assertTrue(normalized)

    def test_documented_parameter_schema_invariants(self):
        repository_root = Path(__file__).resolve().parents[1]
        example = (
            repository_root
            / "examples"
            / "run"
            / "dp2.x-lammps-cp2k"
            / "param_CH4_deepmd-kit-2.0.1.json"
        )
        parameter_data = json.loads(example.read_text())
        normalized = normalize(
            run_jdata_arginfo(),
            copy.deepcopy(parameter_data),
            strict_check=False,
        )

        self.assertEqual(normalized["model_devi_engine"], "lammps")
        self.assertEqual(normalized["train_backend"], "tensorflow")

        without_mass_map = copy.deepcopy(parameter_data)
        without_mass_map.pop("mass_map", None)
        normalized_without_mass_map = normalize(
            run_jdata_arginfo(),
            without_mass_map,
            strict_check=False,
        )
        self.assertEqual(normalized_without_mass_map["mass_map"], "auto")

        invalid_fp_style = copy.deepcopy(parameter_data)
        invalid_fp_style["fp_style"] = "none"
        with self.assertRaises(ArgumentValueError):
            normalize(
                run_jdata_arginfo(),
                invalid_fp_style,
                strict_check=False,
            )

        missing_model_devi_skip = copy.deepcopy(parameter_data)
        missing_model_devi_skip.pop("model_devi_skip")
        with self.assertRaises(ArgumentKeyError):
            normalize(
                run_jdata_arginfo(),
                missing_model_devi_skip,
                strict_check=False,
            )

        flat_sys_configs = copy.deepcopy(parameter_data)
        flat_sys_configs["sys_configs"] = ["POSCAR"]
        with self.assertRaises(ArgumentTypeError):
            normalize(
                run_jdata_arginfo(),
                flat_sys_configs,
                strict_check=False,
            )

        for backend in ("tensorflow", "pytorch"):
            with self.subTest(backend=backend):
                backend_data = copy.deepcopy(parameter_data)
                backend_data["train_backend"] = backend
                normalized_backend = normalize(
                    run_jdata_arginfo(),
                    backend_data,
                    strict_check=False,
                )
                self.assertEqual(normalized_backend["train_backend"], backend)

    def test_linked_machine_examples_use_current_schema(self):
        repository_root = Path(__file__).resolve().parents[1]
        examples_root = repository_root / "examples" / "machine"
        examples = linked_json_examples(repository_root, examples_root)
        self.assertTrue(examples, "No examples/machine/*.json links found")

        for example in examples:
            with self.subTest(example=example):
                machine_data = json.loads(example.read_text())
                stage_values = {}
                for stage in ("train", "model_devi", "fp"):
                    self.assertIn(stage, machine_data)
                    stage_data = machine_data[stage]
                    if isinstance(stage_data, list):
                        self.assertTrue(stage_data)
                        stage_data = stage_data[0]
                    self.assertIsInstance(stage_data, dict)
                    stage_values[stage] = stage_data

                    machine = stage_data.get("machine")
                    resources = stage_data.get("resources")
                    command = stage_data.get("command")
                    self.assertIsInstance(machine, dict)
                    self.assertTrue(machine)
                    self.assertIsInstance(resources, dict)
                    self.assertTrue(resources)
                    self.assertIsInstance(command, str)
                    self.assertTrue(command.strip())
                    for field in ("batch_type", "context_type", "local_root"):
                        self.assertIn(field, machine)
                        self.assertTrue(machine[field])
                    for field in (
                        "number_node",
                        "cpu_per_node",
                        "gpu_per_node",
                        "group_size",
                    ):
                        self.assertIn(field, resources)

                converted = convert_mdata(machine_data)
                for stage in ("train", "model_devi", "fp"):
                    self.assertIn(f"{stage}_machine", converted)
                    self.assertIn(f"{stage}_resources", converted)
                    self.assertIn(f"{stage}_command", converted)
                    converted_machine = converted[f"{stage}_machine"]
                    converted_resources = converted[f"{stage}_resources"]
                    converted_command = converted[f"{stage}_command"]
                    self.assertIsInstance(converted_machine, dict)
                    self.assertTrue(converted_machine)
                    self.assertIsInstance(converted_resources, dict)
                    self.assertTrue(converted_resources)
                    self.assertIsInstance(converted_command, str)
                    self.assertTrue(converted_command.strip())
                    self.assertEqual(
                        converted_machine,
                        stage_values[stage]["machine"],
                    )
                    self.assertEqual(
                        converted_resources,
                        stage_values[stage]["resources"],
                    )
                    self.assertEqual(
                        converted_command,
                        stage_values[stage]["command"],
                    )


if __name__ == "__main__":
    unittest.main()

import json
import re
import unittest
from pathlib import Path

from dpgen.generator.arginfo import run_jdata_arginfo
from dpgen.remote.decide_machine import convert_mdata
from dpgen.util import normalize

LINK_PATTERN = re.compile(r"\[[^]]+\]\(([^)]+)\)")


class TestDPGenRunSkill(unittest.TestCase):
    def test_all_repository_relative_links_exist(self):
        repository_root = Path(__file__).resolve().parents[1]
        skill_path = repository_root / "skills" / "dpgen-run" / "SKILL.md"

        missing = []
        for target in LINK_PATTERN.findall(skill_path.read_text()):
            if "://" in target or target.startswith("#"):
                continue
            resolved_target = (skill_path.parent / target.split("#", 1)[0]).resolve()
            try:
                resolved_target.relative_to(repository_root)
            except ValueError:
                missing.append(target)
            else:
                if not resolved_target.exists():
                    missing.append(target)

        self.assertEqual(
            missing,
            [],
            "Missing repository-relative links in dpgen-run skill: "
            + ", ".join(missing),
        )

    def test_linked_parameter_examples_use_current_schema(self):
        repository_root = Path(__file__).resolve().parents[1]
        skill_path = repository_root / "skills" / "dpgen-run" / "SKILL.md"
        parameter_examples_root = repository_root / "examples" / "run"
        examples = set()
        for target in LINK_PATTERN.findall(skill_path.read_text()):
            target_path = target.split("#", 1)[0]
            if not target_path.endswith(".json"):
                continue
            resolved_target = (skill_path.parent / target_path).resolve()
            try:
                resolved_target.relative_to(parameter_examples_root)
            except ValueError:
                continue
            examples.add(resolved_target)
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

    def test_linked_machine_examples_use_current_schema(self):
        repository_root = Path(__file__).resolve().parents[1]
        skill_path = repository_root / "skills" / "dpgen-run" / "SKILL.md"
        machine_examples_root = repository_root / "examples" / "machine"
        examples = set()
        for target in LINK_PATTERN.findall(skill_path.read_text()):
            target_path = target.split("#", 1)[0]
            if not target_path.endswith(".json"):
                continue
            resolved_target = (skill_path.parent / target_path).resolve()
            try:
                resolved_target.relative_to(machine_examples_root)
            except ValueError:
                continue
            examples.add(resolved_target)
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
                    self.assertEqual(converted_machine, stage_values[stage]["machine"])
                    self.assertEqual(
                        converted_resources, stage_values[stage]["resources"]
                    )
                    self.assertEqual(converted_command, stage_values[stage]["command"])


if __name__ == "__main__":
    unittest.main()

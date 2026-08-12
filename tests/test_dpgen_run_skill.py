import json
import re
import unittest
from pathlib import Path

from dpgen.remote.decide_machine import convert_mdata


class TestDPGenRunSkill(unittest.TestCase):
    def test_all_repository_relative_links_exist(self):
        repository_root = Path(__file__).resolve().parents[1]
        skill_path = repository_root / "skills" / "dpgen-run" / "SKILL.md"
        link_pattern = re.compile(r"\[[^]]+\]\(([^)]+)\)")

        missing = []
        for target in link_pattern.findall(skill_path.read_text()):
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

    def test_linked_machine_examples_use_current_schema(self):
        repository_root = Path(__file__).resolve().parents[1]
        examples = [
            repository_root
            / "examples"
            / "machine"
            / "DeePMD-kit-1.x"
            / "machine-local.json",
            repository_root
            / "examples"
            / "machine"
            / "DeePMD-kit-1.x"
            / "machine-lsf-slurm-cp2k.json",
        ]

        for example in examples:
            with self.subTest(example=example):
                machine_data = json.loads(example.read_text())
                for stage in ("train", "model_devi", "fp"):
                    self.assertIn(stage, machine_data)
                converted = convert_mdata(machine_data)
                for stage in ("train", "model_devi", "fp"):
                    self.assertIn(f"{stage}_machine", converted)
                    self.assertIn(f"{stage}_resources", converted)
                    self.assertIn(f"{stage}_command", converted)


if __name__ == "__main__":
    unittest.main()

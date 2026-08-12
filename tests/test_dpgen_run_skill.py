import re
import unittest
from pathlib import Path


class TestDPGenRunSkill(unittest.TestCase):
    def test_all_repository_relative_links_exist(self):
        repository_root = Path(__file__).resolve().parents[1]
        skill_path = repository_root / "skills" / "dpgen-run" / "SKILL.md"
        link_pattern = re.compile(r"\[[^]]+\]\(([^)]+)\)")

        missing = []
        for target in link_pattern.findall(skill_path.read_text()):
            if "://" in target or target.startswith("#"):
                continue
            target_path = target.split("#", 1)[0]
            if not (skill_path.parent / target_path).resolve().exists():
                missing.append(target)

        self.assertEqual(
            missing,
            [],
            "Missing repository-relative links in dpgen-run skill: "
            + ", ".join(missing),
        )


if __name__ == "__main__":
    unittest.main()

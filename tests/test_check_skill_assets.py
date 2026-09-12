"""Validate bundled skill assets against DP-GEN's public input schemas."""

import json
import unittest
from pathlib import Path

from dpgen.simplify.arginfo import simplify_jdata_arginfo, simplify_mdata_arginfo
from dpgen.util import normalize


class TestSimplifySkillAssets(unittest.TestCase):
    """Ensure every simplify skill JSON asset remains usable as a template."""

    assets = Path(__file__).parent.parent / "skills" / "dpgen-simplify" / "assets"

    def test_json_assets_match_simplify_schemas(self):
        """Normalize parameter and machine assets with their production schemas."""
        for path in sorted(self.assets.glob("*.json")):
            with self.subTest(path=path.name):
                with path.open() as file:
                    data = json.load(file)

                arginfo = (
                    simplify_jdata_arginfo()
                    if path.name.startswith("param.")
                    else simplify_mdata_arginfo()
                )
                normalize(arginfo, data)

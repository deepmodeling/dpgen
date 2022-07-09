"""This module ensures input in the examples directory
could pass the argument checking.
"""
import unittest
import json
from pathlib import Path

from dpgen.util import normalize
from dpgen.data.arginfo import (
    init_reaction_jdata_arginfo,
)
from dpgen.simplify.arginfo import (
    simplify_jdata_arginfo,
)

init_reaction_jdata = init_reaction_jdata_arginfo()
simplify_jdata = simplify_jdata_arginfo()

# directory of examples
p_examples = Path(__file__).parent.parent / "examples"

# input_files : tuple[tuple[Argument, Path]]
#   tuple of example list
input_files = (
    (init_reaction_jdata, p_examples / "init" / "reaction.json"),
    (simplify_jdata, p_examples / "simplify" / "qm7.json"),
)


class TestExamples(unittest.TestCase):
    def test_arguments(self):
        for arginfo, fn in input_files:
            fn = str(fn)
            with self.subTest(fn=fn):
                with open(fn) as f:
                    data = json.load(f)
                normalize(arginfo, data)

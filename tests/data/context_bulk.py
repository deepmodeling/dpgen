import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from dpgen.data.gen import *  # noqa: F403

param_file = "alloy.json"
abacus_param_file = "CuW.json"

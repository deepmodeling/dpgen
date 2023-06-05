import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from dpgen.data.gen import *

param_file = "alloy.json"
abacus_param_file = "CuW.json"

abacus_ref_Cu_coord = 3.76 * np.array(
    [[0.5, 0, 0.5], [0.5, 0, 1.5], [0.5, 1, 0.5], [0.5, 1, 1.5]]
)

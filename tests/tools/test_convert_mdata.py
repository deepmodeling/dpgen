import os,sys,json
import unittest

test_dir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.insert(0, os.path.join(test_dir, '..'))
__package__ = 'tools'
from dpgen.remote.decide_machine import convert_mdata
from .context import machine_file
from .context import setUpModule
class TestConvertMdata(unittest.TestCase):
    def test_convert_mdata (self):
        setUpModule()
        mdata = json.load(open(machine_file))
        mdata = convert_mdata(mdata, ["fp"])
        self.assertEqual(mdata["fp_command"], "vasp_std")
        self.assertEqual(mdata["fp_group_size"], 8)
        self.assertEqual(mdata["fp_machine"]["batch_type"], "PBS")
        self.assertEqual(mdata["fp_user_forward_files"], ["vdw_kernel.bindat"])

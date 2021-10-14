import os,sys,json,glob,shutil,textwrap
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import get_dumped_forces

class TestGetDumpForce(unittest.TestCase):
    def setUp(self):
        file_content = textwrap.dedent("""\
ITEM: TIMESTEP
40
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
-2.9180686220264818e-04 8.0855380329747089e+00 1.4011011277606830e-07
-2.9198257591541018e-04 8.0855378881632269e+00 3.3202396460852749e-08
-2.9180686326490957e-04 8.0855378891632768e+00 -1.7571268247505500e-07
ITEM: ATOMS id type x y z fx fy fz
1 1 2.09532 8.19528 2.00538 -0.00569269 -0.0200373 -0.0342394
2 1 -0.0727384 4.01773 4.05582 -0.0297083 0.0817184 0.0722508
""")
        with open('tmp.dump', 'w') as fp:
            fp.write(file_content)
        self.expected_f = [ -0.00569269, -0.0200373, -0.0342394, -0.0297083, 0.0817184, 0.0722508]

    def tearDown(self):
        if os.path.isfile('tmp.dump'):
            os.remove('tmp.dump')

    def test_read_dump(self):
        ff = get_dumped_forces('tmp.dump')
        self.assertEqual(ff.shape, (2, 3))
        ff = ff.reshape([-1])
        for ii in range(6):
            self.assertAlmostEqual(ff[ii], self.expected_f[ii])

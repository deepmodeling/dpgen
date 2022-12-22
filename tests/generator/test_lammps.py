import os,sys,json,glob,shutil,textwrap
import dpdata
import numpy as np
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'generator'
from .context import get_dumped_forces
from .context import get_all_dumped_forces

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

class TestGetDumpForce(unittest.TestCase):
    def setUp(self):
        file_content = textwrap.dedent("""\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
ITEM: ATOMS id type x y z fx fy fz
1 1 5.38154 4.06861 3.60573 0.000868817 -0.00100822 -0.000960258
2 2 3.9454 4.80321 4.38469 0.000503458 -0.000374043 -9.15676e-05
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
0.0000000000000000e+00 1.0000000000000000e+01 0.0000000000000000e+00
ITEM: ATOMS id type x y z fx fy fz
1 1 5.35629 3.93297 3.70556 -0.125424 0.0481604 -0.0833015
2 2 3.93654 4.79972 4.48179 0.134843 -0.0444238 -0.143111
""")
        with open('tmp.dump', 'w') as fp:
            fp.write(file_content)
        self.expected_f = [ 0.000868817 , -0.00100822 , -0.000960258 , 0.000503458 , -0.000374043 , -9.15676e-05 , -0.125424 , 0.0481604 , -0.0833015 , 0.134843 , -0.0444238 , -0.143111]
    def tearDown(self):
        if os.path.isfile('tmp.dump'):
            os.remove('tmp.dump')

    def test_read_all_dump(self):
        ff = get_all_dumped_forces('tmp.dump')
        ff = np.array(ff)
        self.assertEqual(ff.shape, (2,2,3))
        ff = ff.reshape([-1])
        for ii in range(12):
            self.assertAlmostEqual(ff[ii], self.expected_f[ii])

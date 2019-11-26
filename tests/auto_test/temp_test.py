import unittest
import filecmp

class TestTemp(unittest.TestCase):
    def test_temp(self):
        n,e,v,s=(1,1.00,2.0,None)
        #n,e,v = ((%d %10.3f %10.3f)%(n,e,v))
        self.assertEqual((n,e,v),(1,1.0,2.0))
    def test_file(self):
        result = './01.eos/Cu/std-fcc/deepmd/ref'
        ref = './01.eos/Cu/std-fcc/deepmd/ref1'
        self.assertTrue(filecmp.cmp(ref,result))

suit=unittest.TestSuite()
suit.addTest(TestTemp("test_temp"))
suit.addTest(TestTemp("test_file"))
runner = unittest.TextTestRunner()
runner.run(suit)

import os,json,glob,shutil,filecmp,uuid,time
import unittest

from context import FinRecord

class TestDispUtils(unittest.TestCase):
    def setUp(self):
        self.njobs = 10
        self.fr = FinRecord('.', self.njobs)

    def tearDown(self):
        if os.path.isfile('fin.record'):
            os.remove('fin.record')

    def test_all_false(self) :
        recd = self.fr.get_record()
        self.assertEqual(recd, [False]*self.njobs)

    def test_write_read(self) :
        recd = self.fr.get_record()
        recd[self.njobs//3] = True
        self.fr.write_record(recd)
        recd1 = self.fr.get_record()
        self.assertEqual(recd, recd1)

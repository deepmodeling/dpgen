import os,sys,json,glob,shutil,uuid,time
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'dispatcher'
# from .context import FinRecord
from .context import _split_tasks
from .context import setUpModule

# class TestFinRecord(unittest.TestCase):
#     def setUp(self):
#         self.njobs = 10
#         self.fr = FinRecord('.', self.njobs)

#     def tearDown(self):
#         if os.path.isfile('fin.record'):
#             os.remove('fin.record')

#     def test_all_false(self) :
#         recd = self.fr.get_record()
#         self.assertEqual(recd, [False]*self.njobs)

#     def test_write_read(self) :
#         recd = self.fr.get_record()
#         recd[self.njobs//3] = True
#         self.fr.write_record(recd)
#         recd1 = self.fr.get_record()
#         self.assertEqual(recd, recd1)

class TestDispatchSplit(unittest.TestCase):
    def test_split(self):
        tasks = [ii for ii in range(10)]
        chunks = _split_tasks(tasks, 5)
        self.assertEqual(chunks, [[0,2,4,6,8],[1,3,5,7,9]])

    def test_split_1(self):
        tasks = [ii for ii in range(13)]
        chunks = _split_tasks(tasks, 5)
        self.assertEqual(chunks, [[0,3,6,9,12],[1,4,7,10],[2,5,8,11]])

        

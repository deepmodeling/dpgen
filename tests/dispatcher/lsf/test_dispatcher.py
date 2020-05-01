import os,sys,json,glob,shutil,uuid,time
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'lsf'
from .context import LocalSession
from .context import LocalContext
from .context import LSF
from .context import JobStatus
from .context import Dispatcher
from .context import my_file_cmp
from .context import setUpModule

class TestDispatcher(unittest.TestCase) :
    def setUp(self) :
        os.makedirs('loc', exist_ok = True)
        os.makedirs('rmt', exist_ok = True)
        os.makedirs('loc/task0', exist_ok = True)
        os.makedirs('loc/task1', exist_ok = True)
        os.makedirs('loc/task2', exist_ok = True)
        for ii in ['loc/task0', 'loc/task1', 'loc/task2']:
            with open(os.path.join(ii, 'test0'),'w') as fp:
                fp.write('this is test0 from ' + ii + '\n')
        work_profile = {'work_path':'rmt'}
        self.disp = Dispatcher(work_profile, 'local', 'lsf')

    @unittest.skipIf(not shutil.which("bsub"), "requires LSF")
    def test_sub_success(self):
        tasks = ['task0', 'task1', 'task2']
        self.disp.run_jobs(None,
                           'cp test0 test1',
                           'loc',
                           tasks,
                           2,
                           [],
                           ['test0'],
                           ['test1', 'hereout.log', 'hereerr.log'],
                           outlog = 'hereout.log',
                           errlog = 'hereerr.log')
        for ii in tasks:            
            my_file_cmp(self,
                        os.path.join('loc', ii, 'test0'),
                        os.path.join('loc', ii, 'test1'))
            self.assertTrue(os.path.isfile(os.path.join('loc', ii, 'hereout.log')))
            self.assertTrue(os.path.isfile(os.path.join('loc', ii, 'hereerr.log')))
            

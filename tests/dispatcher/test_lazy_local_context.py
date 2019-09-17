import os,json,glob,shutil,uuid,time
import unittest
from pathlib import Path

from .context import LazyLocalContext
from .context import setUpModule

class TestLazyLocalContext(unittest.TestCase):
    def setUp(self) :
        os.makedirs('loc', exist_ok = True)
        os.makedirs('loc/task0', exist_ok = True)
        os.makedirs('loc/task1', exist_ok = True)
        for ii in ['loc/task0', 'loc/task1']:
            with open(os.path.join(ii, 'test0'),'w') as fp:
                fp.write(str(uuid.uuid4()))
            with open(os.path.join(ii, 'test1'),'w') as fp:
                fp.write(str(uuid.uuid4()))
            with open(os.path.join(ii, 'test2'),'w') as fp:
                fp.write(str(uuid.uuid4()))
            os.makedirs(os.path.join(ii, 'dir0'), exist_ok = True)

    def tearDown(self):
        shutil.rmtree('loc')

    def test_upload(self) :
        self.job  = LazyLocalContext('loc', None)
        self.job1 = LazyLocalContext('loc', None, job_uuid = self.job.job_uuid)
        tasks = ['task0', 'task1']
        files = ['test0', 'test1']
        self.job.upload(tasks, files)
        self.job1.upload(tasks, files)

    def test_download(self):        
        # upload files
        self.job  = LazyLocalContext('loc', None)
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        self.job.download(tasks, ['test0', 'dir0'])
                
    def test_block_call(self) :
        self.job  = LazyLocalContext('loc', None)
        tasks = ['task0', 'task1']
        files = ['test0', 'test1']
        self.job.upload(tasks, files)
        # ls
        code, stdin, stdout, stderr = self.job.block_call('ls')
        self.assertEqual(stdout.read().decode('utf-8'), 'task0\ntask1\n')
        self.assertEqual(stdout.readlines(), ['task0\n','task1\n'])
        self.assertEqual(code, 0)
        code, stdin, stdout, stderr = self.job.block_call('ls a')
        self.assertEqual(code, 2)
        # self.assertEqual(stderr.read().decode('utf-8'), 'ls: cannot access a: No such file or directory\n')
        err_msg = stderr.read().decode('utf-8')
        self.assertTrue('ls: cannot access' in err_msg)
        self.assertTrue('No such file or directory\n' in err_msg)

    def test_block_checkcall(self) :
        self.job  = LazyLocalContext('loc', None)
        tasks = ['task0', 'task1']
        files = ['test0', 'test1']
        self.job.upload(tasks, files)
        # ls
        stdin, stdout, stderr = self.job.block_checkcall('ls')
        self.assertEqual(stdout.read().decode('utf-8'), 'task0\ntask1\n')
        self.assertEqual(stdout.readlines(), ['task0\n','task1\n'])
        with self.assertRaises(RuntimeError):
            stdin, stdout, stderr = self.job.block_checkcall('ls a')
            
    def test_file(self) :
        self.job = LazyLocalContext('loc', None)
        self.assertFalse(self.job.check_file_exists('aaa'))
        tmp = str(uuid.uuid4())
        self.job.write_file('aaa', tmp)
        self.assertTrue(self.job.check_file_exists('aaa'))
        tmp1 = self.job.read_file('aaa')
        self.assertEqual(tmp, tmp1)
        

    def test_call(self) :
        self.job = LazyLocalContext('loc', None)
        proc = self.job.call('sleep 3')
        self.assertFalse(self.job.check_finish(proc))
        time.sleep(1)
        self.assertFalse(self.job.check_finish(proc))
        time.sleep(2.5)
        self.assertTrue(self.job.check_finish(proc))
        r,o,e=self.job.get_return(proc)
        self.assertEqual(r, 0)
        self.assertEqual(o.read(), b'')
        self.assertEqual(e.read(), b'')
        r,o,e=self.job.get_return(proc)
        self.assertEqual(r, 0)
        self.assertEqual(o, None)
        self.assertEqual(e, None)


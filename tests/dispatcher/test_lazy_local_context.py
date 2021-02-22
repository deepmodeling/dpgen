import os,sys,json,glob,shutil,uuid,time
import unittest
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'dispatcher'
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
                
    def test_download_check_mark(self):        
        # upload files
        self.job  = LazyLocalContext('loc', None)
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        record_uuid = []
        # generate extra donwload files
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6'):
                    continue
                with open(os.path.join('loc',ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        self.job.download(tasks, ['test6', 'test7', 'dir1'], check_exists = True, mark_failure = True)
        # check dlded
        cc = 0
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6') :
                    self.assertFalse(os.path.exists(os.path.join('loc', ii, jj)), 
                                     msg = 'found ' + os.path.join('loc', ii, jj))
                    self.assertTrue(os.path.exists(os.path.join('loc', ii, 'tag_failure_download_%s' % jj)), 
                                    msg = 'failed to find ' + os.path.join('loc', ii, 'tag_failure_download_%s' % jj))
                    continue
                with open(os.path.join('loc',ii,jj), 'r') as fp:
                    tmp = fp.read()
                    self.assertEqual(tmp, record_uuid[cc])
                    cc += 1
        for ii in tasks :
            for jj in ['dir1'] :
                self.assertFalse(os.path.exists(os.path.join('loc', ii, jj)))
                self.assertTrue(os.path.exists(os.path.join('loc', ii, 'tag_failure_download_%s' % jj)))


    def test_download_check_nomark(self):        
        # upload files
        self.job  = LazyLocalContext('loc', None)
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        record_uuid = []
        # generate extra donwload files
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6'):
                    continue
                with open(os.path.join('loc',ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        self.job.download(tasks, ['test6', 'test7', 'dir1'], check_exists = True, mark_failure = False)
        # check dlded
        cc = 0
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6') :
                    self.assertFalse(os.path.exists(os.path.join('loc', ii, jj)), 
                                     msg = 'found ' + os.path.join('loc', ii, jj))
                    self.assertFalse(os.path.exists(os.path.join('loc', ii, 'tag_failure_download_%s' % jj)), 
                                     msg = 'found ' + os.path.join('loc', ii, 'tag_failure_download_%s' % jj))
                    continue
                with open(os.path.join('loc',ii,jj), 'r') as fp:
                    tmp = fp.read()
                    self.assertEqual(tmp, record_uuid[cc])
                    cc += 1
        for ii in tasks :
            for jj in ['dir1'] :
                self.assertFalse(os.path.exists(os.path.join('loc', ii, jj)))
                self.assertFalse(os.path.exists(os.path.join('loc', ii, 'tag_failure_download_%s' % jj)))
        

                
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
        proc = self.job.call('sleep 1.5')
        self.assertFalse(self.job.check_finish(proc))
        time.sleep(1)
        self.assertFalse(self.job.check_finish(proc))
        time.sleep(2.5)
        self.assertTrue(self.job.check_finish(proc))
        r,o,e=self.job.get_return(proc)
        self.assertEqual(r, 0)
        self.assertEqual(o.read(), b'')
        self.assertEqual(e.read(), b'')
        # r,o,e=self.job.get_return(proc)
        # self.assertEqual(r, 0)
        # self.assertEqual(o, None)
        # self.assertEqual(e, None)


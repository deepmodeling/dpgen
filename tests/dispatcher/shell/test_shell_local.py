import os,sys,json,glob,shutil,uuid,time
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'shell'
from .context import LocalSession
from .context import LocalContext
from .context import Shell
from .context import JobStatus
from .context import my_file_cmp
from .context import setUpModule

class TestShell(unittest.TestCase) :
    def setUp(self) :
        os.makedirs('loc', exist_ok = True)
        os.makedirs('rmt', exist_ok = True)
        os.makedirs('loc/task0', exist_ok = True)
        os.makedirs('loc/task1', exist_ok = True)
        for ii in ['loc/task0', 'loc/task1']:
            with open(os.path.join(ii, 'test0'),'w') as fp:
                fp.write(str(uuid.uuid4()))
        work_profile = LocalSession({'work_path':'rmt'})
        self.ctx = LocalContext('loc', work_profile)
        self.shell = Shell(self.ctx)

    def tearDown(self):
        shutil.rmtree('loc')
        shutil.rmtree('rmt')
        if os.path.exists('dpgen.log'):
            os.remove('dpgen.log')
        if os.path.exists('run.sub'):
            os.remove('run.sub')
        if os.path.exists('run.sub.1'):
            os.remove('run.sub.1')

    def test_manual_cuda_devices(self):
        job_dirs = ['task0', 'task1']
        res = {'manual_cuda_devices': 3}
        ret = self.shell.sub_script(job_dirs, ['touch test1', 'touch test2'], res = res)
        with open('run.sub.gpu', 'w') as fp:
            fp.write(ret)        
            
    def test_manual_cuda_multiplicity(self):
        job_dirs = ['task0', 'task1', 'task2', 'task3']
        res = {'manual_cuda_devices': 2, 'manual_cuda_multiplicity': 2}
        ret = self.shell.sub_script(job_dirs, ['touch test1', 'touch test2'], res = res)
        with open('run.sub.gpu.multi', 'w') as fp:
            fp.write(ret)

    def test_gen_sub_script(self):
        job_dirs = ['task0', 'task1']
        self.shell.context.upload(job_dirs, ['test0'])
        ret = self.shell.sub_script(job_dirs, ['touch test1', 'touch test2'])
        with open('run.sub', 'w') as fp:
            fp.write(ret)
        ret1 = self.shell.sub_script(job_dirs, ['touch', 'touch'], args = [['test1 ', 'test1 '], ['test2 ', 'test2 ']])
        with open('run.sub.1', 'w') as fp:
            fp.write(ret1)
        time.sleep(1)
        my_file_cmp(self, 'run.sub.1', 'run.sub')
        # with open('run.sub', 'w') as fp:
        #     fp.write(ret)

    def test_sub_success(self) :
        job_dirs = ['task0', 'task1']
        self.shell.context.upload(job_dirs, ['test0'])
        self.shell.submit(job_dirs, ['touch test1', 'touch test2'])
        while True:
            ret = self.shell.check_status()
            if ret == JobStatus.finished  :
                break
            time.sleep(1)        
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_0_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_1_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_0_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_1_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, '%s_tag_finished' % self.shell.context.job_uuid)))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test2')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test2')))


    def test_sub_scancel(self) :
        job_dirs = ['task0', 'task1']
        self.shell.context.upload(job_dirs, ['test0'])
        # sub        
        self.shell.submit(job_dirs, ['touch test1', 'sleep 10'])
        while True:
            ret = self.shell.check_status()
            if ret == JobStatus.finished  :
                raise RuntimeError('should not finished')
            if ret == JobStatus.running :
                # wait for file writing
                time.sleep(2)
                # kill job
                self.shell.context.kill(self.shell.proc)
                break
            time.sleep(1)
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_0_finished')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_1_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_0_finished')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_1_finished')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, '%s_tag_finished' % self.shell.context.job_uuid)))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test1')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test2')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test2')))
        # sub restart
        self.shell.submit(job_dirs, ['rm test1', 'touch test2'], restart = True)
        while True:
            ret = self.shell.check_status()
            if ret == JobStatus.finished  :
                break
            time.sleep(1)
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_0_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_1_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_0_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_1_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, '%s_tag_finished' % self.shell.context.job_uuid)))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test2')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test2')))
        

import os,sys,json,glob,shutil,uuid,time
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'lsf'
from .context import LocalSession
from .context import LocalContext
from .context import LSF
from .context import JobStatus
from .context import setUpModule

class TestLSF(unittest.TestCase) :
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
        self.lsf = LSF(self.ctx)

    def tearDown(self):
        shutil.rmtree('loc')
        shutil.rmtree('rmt')
        if os.path.exists('dpgen.log'):
            os.remove('dpgen.log')

    def test_gen_sub_script(self):
        job_dirs = ['task0', 'task1']
        self.lsf.context.upload(job_dirs, ['test0'])
        ret = self.lsf.sub_script(job_dirs, ['touch test1', 'touch test2'])
        self.lsf.context.write_file('run.sub', ret)
        with open('run.sub', 'w') as fp:
            fp.write(ret)            

    @unittest.skipIf(not shutil.which("bsub"), "requires LSF")
    def test_sub_success(self) :
         job_dirs = ['task0', 'task1']
         self.lsf.context.upload(job_dirs, ['test0'])
         self.lsf.submit(job_dirs, ['touch test1', 'touch test2'])
         while True:
             ret = self.lsf.check_status()
             if ret == JobStatus.finished  :
                 break
             time.sleep(1)        
         self.assertTrue(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/tag_0_finished')))
         self.assertTrue(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/tag_1_finished')))
         self.assertTrue(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/tag_0_finished')))
         self.assertTrue(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/tag_1_finished')))
         self.assertTrue(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'tag_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/test1')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/test1')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/test2')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/test2')))

    @unittest.skipIf(not shutil.which("bsub"), "requires LSF")
    def test_sub_bkill(self) :
         job_dirs = ['task0', 'task1']
         self.lsf.context.upload(job_dirs, ['test0'])
         # sub        
         self.lsf.submit(job_dirs, ['touch test1', 'sleep 10'])
         while True:
             ret = self.lsf.check_status()
             if ret == JobStatus.finished  :
                 raise RuntimeError('should not finished')
             if ret == JobStatus.running :
                # wait for file writing
                 time.sleep(2)
                 job_id = self.lsf._get_job_id()
                 os.system('bkill ' + job_id)
                 break
             time.sleep(1)
         while True:
             ret = self.lsf.check_status()
             if ret == JobStatus.terminated  :
                 break
             time.sleep(1)
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/tag_0_finished')))
         self.assertFalse(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/tag_1_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/tag_0_finished')))
         self.assertFalse(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/tag_1_finished')))
         self.assertFalse(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'tag_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/test1')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/test1')))
         self.assertFalse(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/test2')))
         self.assertFalse(os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/test2')))
         # sub restart
         self.lsf.submit(job_dirs, ['rm test1', 'touch test2'], restart = True)
         while True:
             ret = self.lsf.check_status()
             if ret == JobStatus.finished  :
                 break
             time.sleep(1)
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/tag_0_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/tag_1_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/tag_0_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/tag_1_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'tag_finished')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/test1')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/test1')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task0/test2')))
         self.assertTrue (os.path.isfile(os.path.join('rmt', self.lsf.context.remote_root, 'task1/test2')))
       

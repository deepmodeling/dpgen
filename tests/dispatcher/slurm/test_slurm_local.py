import os,sys,json,glob,shutil,uuid,time
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'slurm'
from .context import LocalSession
from .context import LocalContext
from .context import Slurm
from .context import JobStatus
from .context import setUpModule

@unittest.skipIf(not shutil.which("sbatch"), "requires Slurm")
class TestSlurm(unittest.TestCase) :
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
        self.slurm = Slurm(self.ctx)

    def tearDown(self):
        shutil.rmtree('loc')
        shutil.rmtree('rmt')
        if os.path.exists('dpgen.log'):
            os.remove('dpgen.log')

    def test_gen_sub_script(self):
        job_dirs = ['task0', 'task1']
        self.slurm.context.upload(job_dirs, ['test0'])
        ret = self.slurm.sub_script(job_dirs, ['touch test1', 'touch test2'])
        self.slurm.context.write_file('run.sub', ret)
        with open('run.sub', 'w') as fp:
            fp.write(ret)            

    def test_sub_success(self) :
        job_dirs = ['task0', 'task1']
        self.slurm.context.upload(job_dirs, ['test0'])
        self.slurm.submit(job_dirs, ['touch test1', 'touch test2'])
        while True:
            ret = self.slurm.check_status()
            if ret == JobStatus.finished  :
                break
            time.sleep(1)        
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/tag_0_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/tag_1_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/tag_0_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/tag_1_finished')))
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, '%s_tag_finished' % self.slurm.context.job_uuid)))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/test2')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/test2')))

    def test_sub_scancel(self) :
        job_dirs = ['task0', 'task1']
        self.slurm.context.upload(job_dirs, ['test0'])
        # sub        
        self.slurm.submit(job_dirs, ['touch test1', 'sleep 10'])
        while True:
            ret = self.slurm.check_status()
            if ret == JobStatus.finished  :
                raise RuntimeError('should not finished')
            if ret == JobStatus.running :
                # wait for file writing
                time.sleep(2)
                job_id = self.slurm._get_job_id()
                os.system('scancel ' + job_id)
                break
            time.sleep(1)
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/tag_0_finished')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/tag_1_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/tag_0_finished')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/tag_1_finished')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, '%s_tag_finished' % self.slurm.context.job_uuid)))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/test1')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/test2')))
        self.assertFalse(os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/test2')))
        # sub restart
        self.slurm.submit(job_dirs, ['rm test1', 'touch test2'], restart = True)
        while True:
            ret = self.slurm.check_status()
            if ret == JobStatus.finished  :
                break
            time.sleep(1)
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/tag_0_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/tag_1_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/tag_0_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/tag_1_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, '%s_tag_finished' % self.slurm.context.job_uuid)))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task0/test2')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.slurm.context.remote_root, 'task1/test2')))
        

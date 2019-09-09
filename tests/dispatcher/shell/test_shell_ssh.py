import os,json,glob,shutil,filecmp,uuid,time,getpass
import unittest
from context import SSHSession
from context import SSHContext
from context import Shell
from context import JobStatus

class TestShell(unittest.TestCase) :
    def setUp(self) :
        os.makedirs('loc', exist_ok = True)
        os.makedirs('rmt', exist_ok = True)
        os.makedirs('loc/task0', exist_ok = True)
        os.makedirs('loc/task1', exist_ok = True)
        for ii in ['loc/task0', 'loc/task1']:
            with open(os.path.join(ii, 'test0'),'w') as fp:
                fp.write(str(uuid.uuid4()))
        ssh_session = SSHSession({'hostname' : 'localhost',
                                  'port': 5566,
                                  'username' : getpass.getuser(),
                                  'work_path' : os.path.join(os.getcwd(), 'rmt')})
        self.ctx = SSHContext('loc', ssh_session)
        self.shell = Shell(self.ctx)

    def tearDown(self):
        shutil.rmtree('loc')
        shutil.rmtree('rmt')
        if os.path.exists('dpgen.log'):
            os.remove('dpgen.log')

    def test_gen_sub_script(self):
        job_dirs = ['task0', 'task1']
        self.shell.context.upload(job_dirs, ['test0'])
        ret = self.shell.sub_script(job_dirs, ['touch test1', 'touch test2'])
        self.shell.context.write_file('run.sub', ret)
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
        self.assertTrue(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'tag_finished')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test1')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test2')))
        self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test2')))


    # def test_sub_scancel(self) :
    #     job_dirs = ['task0', 'task1']
    #     self.shell.context.upload(job_dirs, ['test0'])
    #     # sub        
    #     self.shell.submit(job_dirs, ['touch test1', 'sleep 10'])
    #     while True:
    #         ret = self.shell.check_status()
    #         if ret == JobStatus.finished  :
    #             raise RuntimeError('should not finished')
    #         if ret == JobStatus.running :
    #             # wait for file writing
    #             time.sleep(2)
    #             # kill job
    #             ##################################################
    #             # problematic killing remotly
    #             ##################################################
    #             self.shell.context.kill(self.shell.proc)
    #             break
    #         time.sleep(1)
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_0_finished')))
    #     self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_1_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_0_finished')))
    #     self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_1_finished')))
    #     self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'tag_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test1')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test1')))
    #     self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test2')))
    #     self.assertFalse(os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test2')))
    #     # sub restart
    #     self.shell.submit(job_dirs, ['rm test1', 'touch test2'], restart = True)
    #     while True:
    #         ret = self.shell.check_status()
    #         if ret == JobStatus.finished  :
    #             break
    #         time.sleep(1)
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_0_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/tag_1_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_0_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/tag_1_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'tag_finished')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test1')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test1')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task0/test2')))
    #     self.assertTrue (os.path.isfile(os.path.join('rmt', self.shell.context.remote_root, 'task1/test2')))
        

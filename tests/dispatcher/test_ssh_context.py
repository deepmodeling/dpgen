import os,sys,json,glob,shutil,uuid,getpass
import unittest
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'dispatcher'
from .context import SSHContext, SSHSession
from .context import setUpModule

class TestSSHContext(unittest.TestCase):
    def setUp(self) :
        os.makedirs('loc', exist_ok = True)
        os.makedirs('loc/task0', exist_ok = True)
        os.makedirs('loc/task1', exist_ok = True)
        for ii in ['loc/task0', 'loc/task1']:
            with open(os.path.join(ii, 'test0'),'w') as fp:
                fp.write(str(uuid.uuid4()))
            with open(os.path.join(ii, 'test1'),'w') as fp:
                fp.write(str(uuid.uuid4()))
            os.makedirs(os.path.join(ii, 'dir0'), exist_ok = True)
            with open(os.path.join(ii, 'dir0', 'test2'),'w') as fp:
                fp.write(str(uuid.uuid4()))
        os.makedirs('rmt', exist_ok = True)
        try :
            self.ssh_session = SSHSession({'hostname' : 'localhost',
                                           'port': 22,
                                           'username' : getpass.getuser(),
                                           'work_path' : os.path.join(os.getcwd(), 'rmt')})
        except Exception:
            # for tianhe-2
            try:
                self.ssh_session = SSHSession({'hostname' : 'localhost',
                                            'port': 5566,
                                            'username' : getpass.getuser(),
                                            'work_path' : os.path.join(os.getcwd(), 'rmt')})
            except Exception:
                self.skipTest("Network error")
        self.job = SSHContext('loc', self.ssh_session)
        self.job1 = SSHContext('loc', self.ssh_session, job_uuid = self.job.job_uuid)
    
    def tearDown(self):
        self.job.close()
        self.job1.close()
        shutil.rmtree('loc')
        shutil.rmtree('rmt')

    def test_upload(self) :
        tasks = ['task0', 'task1']
        files = ['test0', 'test1']
        self.job.upload(tasks, files)
        for ii in tasks :
            for jj in files :
                locf = os.path.join('loc', ii, jj)
                rmtf = os.path.join('rmt', self.job.job_uuid, ii, jj)
                with open(locf) as fp:
                    locs = fp.read()
                with open(rmtf) as fp:
                    rmts = fp.read()
                self.assertEqual(locs, rmts)
        self.job.upload(tasks, ['dir0'])
        for ii in tasks :
            for jj in ['dir0'] :
                for kk in ['test2'] :
                    locf = os.path.join('loc', ii, jj, kk)
                    rmtf = os.path.join('rmt', self.job.job_uuid, ii, jj, kk)
                    with open(locf) as fp:
                        locs = fp.read()
                    with open(rmtf) as fp:
                        rmts = fp.read()
                    self.assertEqual(locs, rmts)


    def test_donwload(self):
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        # generate extra donwload files
        record_uuid = []
        for ii in tasks :
            for jj in ['test4', 'test5'] :
                with open(os.path.join('rmt',self.job.job_uuid,ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        # generate extra donwload dirs and files
        for ii in tasks :
            for jj in ['dir1'] :
                os.makedirs(os.path.join('rmt',self.job.job_uuid,ii,jj), exist_ok=False)
                for kk in ['test6']:
                    with open(os.path.join('rmt',self.job.job_uuid,ii,jj,kk), 'w') as fp:
                        tmp = str(uuid.uuid4())
                        fp.write(tmp)
                        record_uuid.append(tmp)                        
        # donwload
        files = ['test4', 'test5', 'dir1']
        self.job.download(tasks, files)
        # check dlded
        cc = 0
        for ii in tasks :
            for jj in ['test4', 'test5'] :
                with open(os.path.join('loc',ii,jj), 'r') as fp:
                    tmp = fp.read()
                    self.assertEqual(tmp, record_uuid[cc])
                    cc += 1
        for ii in tasks :
            for jj in ['dir1'] :
                for kk in ['test6']:
                    with open(os.path.join('loc',ii,jj,kk), 'r') as fp:
                        tmp = fp.read()
                        self.assertEqual(tmp, record_uuid[cc])
                        cc += 1


    def test_donwload_check_mark(self):
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        # generate extra donwload files
        record_uuid = []
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6'):
                    continue
                with open(os.path.join('rmt',self.job.job_uuid,ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        # donwload
        files = ['test6', 'test7', 'dir1']
        self.job.download(tasks, files, check_exists = True, mark_failure = True)
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


    def test_donwload_check_nomark(self):
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        # generate extra donwload files
        record_uuid = []
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if ii == 'task1' and jj == 'test7' :
                    continue
                if ii == 'task0' and jj == 'test6' :
                    continue
                with open(os.path.join('rmt',self.job.job_uuid,ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        # donwload
        files = ['test6', 'test7', 'dir1']
        self.job.download(tasks, files, check_exists = True, mark_failure = False)
        # check dlded
        cc = 0
        for ii in tasks :
            for jj in ['test6', 'test7'] :
                if ii == 'task1' and jj == 'test7' :
                    self.assertFalse(os.path.exists(os.path.join('loc', ii, jj)), 
                                     msg = 'found ' + os.path.join('loc', ii, jj))
                    self.assertFalse(os.path.exists(os.path.join('loc', ii, 'tag_failure_download_%s' % jj)), 
                                     msg = 'found ' + os.path.join('loc', ii, 'tag_failure_download_%s' % jj))
                    continue
                if ii == 'task0' and jj == 'test6' :
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
        tasks = ['task0', 'task1']
        files = ['test0', 'test1']
        self.job.upload(tasks, files)        
        # ls
        code, stdin, stdout, stderr = self.job.block_call('ls')
        self.assertEqual(stdout.read(), b'task0\ntask1\n')
        self.assertEqual(code, 0)
        code, stdin, stdout, stderr = self.job.block_call('ls')
        self.assertEqual(stdout.readlines(), ['task0\n','task1\n'])
        code, stdin, stdout, stderr = self.job.block_call('ls a')
        self.assertEqual(code, 2)
        # self.assertEqual(stderr.read().decode('utf-8'), 'ls: cannot access a: No such file or directory\n')
        err_msg = stderr.read().decode('utf-8')
        self.assertTrue('ls: cannot access' in err_msg)
        self.assertTrue('No such file or directory\n' in err_msg)
                        
    def test_block_checkcall(self) :
        tasks = ['task0', 'task1']
        files = ['test0', 'test1']
        self.job.upload(tasks, files)        
        # ls
        stdin, stdout, stderr = self.job.block_checkcall('ls')
        self.assertEqual(stdout.read(), b'task0\ntask1\n')
        stdin, stdout, stderr = self.job.block_checkcall('ls')
        self.assertEqual(stdout.readlines(), ['task0\n','task1\n'])
        with self.assertRaises(RuntimeError):
            stdin, stdout, stderr = self.job.block_checkcall('ls a')

    def test_file(self) :
        self.assertFalse(self.job.check_file_exists('aaa'))
        tmp = str(uuid.uuid4())
        self.job.write_file('aaa', tmp)
        self.assertTrue(self.job.check_file_exists('aaa'))
        tmp1 = self.job.read_file('aaa')
        self.assertEqual(tmp, tmp1)



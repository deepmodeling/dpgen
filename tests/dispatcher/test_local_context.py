import os,sys,json,glob,shutil,uuid,time
import unittest
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
__package__ = 'dispatcher'
from .context import LocalContext, LocalSession
from .context import setUpModule
from .context import _identical_files

class TestIdFile(unittest.TestCase) :
    def test_id(self) :
        with open('f0', 'w') as fp:
            fp.write('foo')
        with open('f1', 'w') as fp:
            fp.write('foo')
        self.assertTrue(_identical_files('f0', 'f1'))
        os.remove('f0')
        os.remove('f1')

    def test_diff(self) :
        with open('f0', 'w') as fp:
            fp.write('foo')
        with open('f1', 'w') as fp:
            fp.write('bar')
        self.assertFalse(_identical_files('f0', 'f1'))
        os.remove('f0')
        os.remove('f1')


class TestLocalContext(unittest.TestCase):
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
            os.makedirs(os.path.join(ii, 'dir2'), exist_ok = True)
            with open(os.path.join(ii, 'dir2', 'dtest0'),'w') as fp:
                fp.write(str(uuid.uuid4()))
        os.makedirs('rmt', exist_ok = True)

    def tearDown(self):
        shutil.rmtree('loc')
        shutil.rmtree('rmt')

    def test_upload_non_exist(self) :
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        tasks = ['task0', 'task1']
        # test uploading non-existing file
        with self.assertRaises(OSError):
            self.job.upload(tasks, ['foo'])

    def test_upload(self) :
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        self.job1 = LocalContext('loc', work_profile, job_uuid = self.job.job_uuid)
        tasks = ['task0', 'task1']
        files = ['test0', 'test1', 'dir2/dtest0']
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
                locf = os.path.join('loc', ii, jj)
                rmtf = os.path.join('rmt', self.job.job_uuid, ii, jj)
                self.assertEqual(os.path.realpath(locf),
                                 os.path.realpath(rmtf))
        self.job1.upload(tasks, files)
        for ii in tasks :
            for jj in files :
                locf = os.path.join('loc', ii, jj)
                rmtf = os.path.join('rmt', self.job.job_uuid, ii, jj)
                with open(locf) as fp:
                    locs = fp.read()
                with open(rmtf) as fp:
                    rmts = fp.read()
                self.assertEqual(locs, rmts)

    def test_dl_f_f(self):
        # no local, no remote
        self.test_download_non_exist()
        
    def test_dl_t_f(self) :
        # has local, no remote
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        tasks = ['task0', 'task1']
        record_uuid = []        
        for ii in tasks :
            for jj in ['dir1'] :
                os.makedirs(os.path.join('loc',ii,jj), exist_ok=False)
                for kk in ['test6', 'test7']:
                    with open(os.path.join('loc',ii,jj,kk), 'w') as fp:
                        tmp = str(uuid.uuid4())
                        fp.write(tmp)
                        record_uuid.append(tmp)
        files = ['dir1']
        self.job.download(tasks, files)
        cc = 0
        for ii in tasks :
            for jj in ['dir1'] :
                for kk in ['test6', 'test7']:
                    with open(os.path.join('loc',ii,jj,kk), 'r') as fp:
                        tmp = fp.read()
                        self.assertEqual(tmp, record_uuid[cc])
                        cc += 1
        
    def test_dl_t_t(self) :
        # has local, has remote
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        tasks = ['task0', 'task1']
        for ii in tasks :
            for jj in ['dir1'] :
                os.makedirs(os.path.join('loc',ii,jj), exist_ok=False)
                for kk in ['test6', 'test7']:
                    with open(os.path.join('loc',ii,jj,kk), 'w') as fp:
                        tmp = str(uuid.uuid4())
                        fp.write(tmp)
        record_uuid = []        
        for ii in tasks :
            for jj in ['dir1'] :
                os.makedirs(os.path.join('rmt', self.job.job_uuid,ii,jj), exist_ok=False)
                for kk in ['test6', 'test7']:
                    with open(os.path.join('rmt', self.job.job_uuid,ii,jj,kk), 'w') as fp:
                        tmp = str(uuid.uuid4())
                        fp.write(tmp)
                        record_uuid.append(tmp)
        files = ['dir1']
        self.job.download(tasks, files)
        cc = 0
        for ii in tasks :
            for jj in ['dir1'] :
                for kk in ['test6', 'test7']:
                    with open(os.path.join('loc',ii,jj,kk), 'r') as fp:
                        tmp = fp.read()
                        self.assertEqual(tmp, record_uuid[cc])
                        cc += 1
        

    def test_download_non_exist(self):
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        tasks = ['task0', 'task1']
        # down load non-existing file
        with self.assertRaises(RuntimeError):
            self.job.download(tasks, ['foo'])        

    def test_download(self):        
        # upload files
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
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
        files = ['test0', 'dir0', 'test4', 'test5', 'dir1']
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
        # check links preserved
        for ii in tasks :
            for jj in ['test0'] :
                locf = os.path.join('loc', ii, jj)
                rmtf = os.path.join('rmt', self.job.job_uuid, ii, jj)
                self.assertEqual(os.path.realpath(locf),
                                 os.path.realpath(rmtf))
        for ii in tasks :
            for jj in ['dir0'] :
                for kk in ['test6'] :
                    locf = os.path.join('loc', ii, jj, kk)
                    rmtf = os.path.join('rmt', self.job.job_uuid, ii, jj, kk)
                    self.assertEqual(os.path.realpath(locf),
                                     os.path.realpath(rmtf))

    def test_download_check_mark(self):        
        # upload files
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        # generate extra donwload files
        record_uuid = []
        for ii in tasks :
            for jj in ['test7', 'test8'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6') :
                    continue
                with open(os.path.join('rmt',self.job.job_uuid,ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        # donwload
        files = ['test7', 'test8', 'dir1']
        self.job.download(tasks, files, check_exists = True, mark_failure = True)
        # check dlded
        cc = 0
        for ii in tasks :
            for jj in ['test7', 'test8'] :
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
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
        tasks = ['task0', 'task1']
        self.job.upload(tasks, ['test0', 'dir0'])
        # generate extra donwload files
        record_uuid = []
        for ii in tasks :
            for jj in ['test7', 'test8'] :
                if (ii == 'task1' and jj == 'test7') or \
                   (ii == 'task0' and jj == 'test6') :
                    continue
                with open(os.path.join('rmt',self.job.job_uuid,ii,jj), 'w') as fp:
                    tmp = str(uuid.uuid4())
                    fp.write(tmp)
                    record_uuid.append(tmp)
        # donwload
        files = ['test7', 'test8', 'dir1']
        self.job.download(tasks, files, check_exists = True, mark_failure = False)
        # check dlded
        cc = 0
        for ii in tasks :
            for jj in ['test7', 'test8'] :
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
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
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
        work_profile = LocalSession({'work_path':'rmt'})
        self.job  = LocalContext('loc', work_profile)
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
        work_profile = LocalSession({'work_path':'rmt'})
        self.job = LocalContext('loc', work_profile)
        self.assertFalse(self.job.check_file_exists('aaa'))
        tmp = str(uuid.uuid4())
        self.job.write_file('aaa', tmp)
        self.assertTrue(self.job.check_file_exists('aaa'))
        tmp1 = self.job.read_file('aaa')
        self.assertEqual(tmp, tmp1)
        

    def test_call(self) :
        work_profile = LocalSession({'work_path':'rmt'})
        self.job = LocalContext('loc', work_profile)
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


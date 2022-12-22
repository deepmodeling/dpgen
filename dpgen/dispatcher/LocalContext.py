import os,shutil,uuid,hashlib
import subprocess as sp
from glob import glob
from dpgen import dlog

class LocalSession (object) :
    def __init__ (self, jdata) :
        self.work_path = os.path.abspath(jdata['work_path'])
        assert(os.path.exists(self.work_path))

    def get_work_root(self) :
        return self.work_path

class SPRetObj(object) :
    def __init__ (self,
                  ret) :
        self.data = ret

    def read(self) :
        return self.data

    def readlines(self) :
        lines = self.data.decode('utf-8').splitlines()
        ret = []
        for aa in lines:
            ret.append(aa+'\n')
        return ret

def _check_file_path(fname) :
    dirname = os.path.dirname(fname)    
    if dirname != "":
        os.makedirs(dirname, exist_ok=True)

def _identical_files(fname0, fname1) :
    with open(fname0) as fp:
        code0 = hashlib.sha1(fp.read().encode('utf-8')).hexdigest()
    with open(fname1) as fp:
        code1 = hashlib.sha1(fp.read().encode('utf-8')).hexdigest()
    return code0 == code1


class LocalContext(object) :
    def __init__ (self,
                  local_root,
                  work_profile,
                  job_uuid = None) :
        """
        work_profile:
        local_root:
        """
        assert(type(local_root) == str)
        self.local_root = os.path.abspath(local_root)
        if job_uuid:
           self.job_uuid=job_uuid
        else:
           self.job_uuid = str(uuid.uuid4())

        self.remote_root = os.path.join(work_profile.get_work_root(), self.job_uuid)
        dlog.debug("local_root is %s"% local_root)
        dlog.debug("remote_root is %s"% self.remote_root)

        os.makedirs(self.remote_root, exist_ok = True)
        
    def get_job_root(self) :
        return self.remote_root

    def upload(self,
               job_dirs,
               local_up_files,
               dereference = True) :
        cwd = os.getcwd()
        for ii in job_dirs :
            local_job = os.path.join(self.local_root, ii)
            remote_job = os.path.join(self.remote_root, ii)
            os.makedirs(remote_job, exist_ok = True)
            os.chdir(remote_job)
            for jj in local_up_files :
                if not os.path.exists(os.path.join(local_job, jj)):
                    os.chdir(cwd)
                    raise OSError('cannot find upload file ' + os.path.join(local_job, jj))
                if os.path.exists(os.path.join(remote_job, jj)) :
                    os.remove(os.path.join(remote_job, jj))
                _check_file_path(jj)
                os.symlink(os.path.join(local_job, jj),
                           os.path.join(remote_job, jj))
        os.chdir(cwd)

    def download(self, 
                 job_dirs,
                 remote_down_files,
                 check_exists = False,
                 mark_failure = True,
                 back_error=False) :
        cwd = os.getcwd()
        for ii in job_dirs :
            local_job = os.path.join(self.local_root, ii)
            remote_job = os.path.join(self.remote_root, ii)
            flist = remote_down_files
            if back_error :
                os.chdir(remote_job)
                flist += glob('error*')                        
                os.chdir(cwd)
            for jj in flist :
                rfile = os.path.join(remote_job, jj)
                lfile = os.path.join(local_job, jj)
                if not os.path.realpath(rfile) == os.path.realpath(lfile) :
                    if (not os.path.exists(rfile)) and (not os.path.exists(lfile)):
                        if check_exists :
                            if mark_failure:
                                with open(os.path.join(self.local_root, ii, 'tag_failure_download_%s' % jj), 'w') as fp: pass
                            else :
                                pass
                        else :
                            raise RuntimeError('do not find download file ' + rfile)
                    elif (not os.path.exists(rfile)) and (os.path.exists(lfile)) :
                        # already downloaded
                        pass
                    elif (os.path.exists(rfile)) and (not os.path.exists(lfile)) :
                        # trivial case, download happily
                        # If the file to be downloaded is a softlink, `cp` should be performed instead of `mv`.
                        # Otherwise, `lfile` is still a file linked to some original file,
                        # and when this file's removed, `lfile` will be invalid.
                        if os.path.islink(rfile):
                            shutil.copyfile(rfile,lfile)
                        else:
                            shutil.move(rfile, lfile)
                    elif (os.path.exists(rfile)) and (os.path.exists(lfile)) :
                        # both exists, replace!
                        dlog.info('find existing %s, replacing by %s' % (lfile, rfile))
                        if os.path.isdir(lfile):
                            shutil.rmtree(lfile, ignore_errors=True)
                        elif os.path.isfile(lfile) or os.path.islink(lfile):
                            os.remove(lfile)
                        shutil.move(rfile, lfile)
                    else :
                        raise RuntimeError('should not reach here!')
                else :
                    # no nothing in the case of linked files
                    pass
        os.chdir(cwd)

    def block_checkcall(self,
                        cmd) :
        cwd = os.getcwd()
        os.chdir(self.remote_root)
        proc = sp.Popen(cmd, shell=True, stdout = sp.PIPE, stderr = sp.PIPE)
        o, e = proc.communicate()
        stdout = SPRetObj(o)
        stderr = SPRetObj(e)
        code = proc.returncode
        if code != 0:
            os.chdir(cwd)        
            raise RuntimeError("Get error code %d in locally calling %s with job: %s ", (code, cmd, self.job_uuid))
        os.chdir(cwd)        
        return None, stdout, stderr
        
    def block_call(self, cmd) :
        cwd = os.getcwd()
        os.chdir(self.remote_root)
        proc = sp.Popen(cmd, shell=True, stdout = sp.PIPE, stderr = sp.PIPE)
        o, e = proc.communicate()
        stdout = SPRetObj(o)
        stderr = SPRetObj(e)
        code = proc.returncode
        os.chdir(cwd)        
        return code, None, stdout, stderr

    def clean(self) :
        shutil.rmtree(self.remote_root, ignore_errors=True)

    def write_file(self, fname, write_str):
        with open(os.path.join(self.remote_root, fname), 'w') as fp :
            fp.write(write_str)

    def read_file(self, fname):
        with open(os.path.join(self.remote_root, fname), 'r') as fp:
            ret = fp.read()
        return ret

    def check_file_exists(self, fname):
        return os.path.isfile(os.path.join(self.remote_root, fname))
        
    def call(self, cmd) :
        cwd = os.getcwd()
        os.chdir(self.remote_root)
        proc = sp.Popen(cmd, shell=True, stdout = sp.PIPE, stderr = sp.PIPE)
        os.chdir(cwd)        
        return proc

    def kill(self, proc):
        proc.kill()

    def check_finish(self, proc):
        return (proc.poll() != None)

    def get_return(self, proc):
        ret = proc.poll()
        if ret is None:
            return None, None, None
        else :
            try:
                o, e = proc.communicate()
                stdout = SPRetObj(o)
                stderr = SPRetObj(e)
            except ValueError:
                stdout = None
                stderr = None
        return ret, stdout, stderr
    
    

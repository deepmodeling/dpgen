import os,shutil,uuid
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
                    raise RuntimeError('cannot file upload file ' + os.path.join(local_job, jj))
                if os.path.exists(os.path.join(remote_job, jj)) :
                    os.remove(os.path.join(remote_job, jj))
                _check_file_path(jj)
                os.symlink(os.path.join(local_job, jj),
                           os.path.join(remote_job, jj))
        os.chdir(cwd)

    def download(self, 
                 job_dirs,
                 remote_down_files,
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
                if not os.path.exists(rfile) :
                    os.chdir(cwd)
                    raise RuntimeError('cannot file download file ' + rfile)
                if not os.path.realpath(rfile) == os.path.realpath(lfile) :
                    shutil.move(rfile, lfile)
                else :
                    # print('skip ' + rfile)
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
        shutil.rmtree(self.remote_root)

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
            except:
                stdout = None
                stderr = None
        return ret, stdout, stderr
    
    

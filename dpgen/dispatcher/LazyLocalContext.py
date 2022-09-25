import os,shutil,uuid
import subprocess as sp
from glob import glob
from dpgen import dlog

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

class LazyLocalContext(object) :
    def __init__ (self,
                  local_root,
                  work_profile = None,
                  job_uuid = None) :
        """
        work_profile:
        local_root:
        """
        assert(type(local_root) == str)
        self.local_root = os.path.abspath(local_root)
        self.remote_root = self.local_root
        if job_uuid:
           self.job_uuid=job_uuid
        else:
           self.job_uuid = str(uuid.uuid4())
        
    def get_job_root(self) :
        return self.local_root

    def upload(self,
               job_dirs,
               local_up_files,
               dereference = True) :
        pass

    def download(self, 
                 job_dirs,
                 remote_down_files,
                 check_exists = False,
                 mark_failure = True,
                 back_error=False) :
        for ii in job_dirs :
            for jj in remote_down_files :
                fname = os.path.join(self.local_root, ii, jj)
                exists = os.path.exists(fname)
                if not exists:
                    if check_exists:
                        if mark_failure:
                            with open(os.path.join(self.local_root, ii, 'tag_failure_download_%s' % jj), 'w') as fp: pass
                        else:
                            pass
                    else:
                        raise OSError('do not find download file ' + fname)


    def block_checkcall(self,
                        cmd) :
        cwd = os.getcwd()
        os.chdir(self.local_root)
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
        os.chdir(self.local_root)
        proc = sp.Popen(cmd, shell=True, stdout = sp.PIPE, stderr = sp.PIPE)
        o, e = proc.communicate()
        stdout = SPRetObj(o)
        stderr = SPRetObj(e)
        code = proc.returncode
        os.chdir(cwd)        
        return code, None, stdout, stderr

    def clean(self) :
        pass

    def write_file(self, fname, write_str):
        with open(os.path.join(self.local_root, fname), 'w') as fp :
            fp.write(write_str)

    def read_file(self, fname):
        with open(os.path.join(self.local_root, fname), 'r') as fp:
            ret = fp.read()
        return ret

    def check_file_exists(self, fname):
        return os.path.isfile(os.path.join(self.local_root, fname))
        
    def call(self, cmd) :
        cwd = os.getcwd()
        os.chdir(self.local_root)
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
    
    

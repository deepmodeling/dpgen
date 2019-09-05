#!/usr/bin/env python
# coding: utf-8

import os, sys, paramiko, json, uuid, tarfile, time, stat, shutil
from glob import glob
from dpgen import dlog

class SSHSession (object) :
    def __init__ (self, jdata) :
        self.remote_profile = jdata
        # with open(remote_profile) as fp :
        #     self.remote_profile = json.load(fp)
        self.remote_host = self.remote_profile['hostname']
        self.remote_port = self.remote_profile['port']
        self.remote_uname = self.remote_profile['username']
        self.remote_password = None
        if 'password' in self.remote_profile :
            self.remote_password = self.remote_profile['password']
        self.remote_workpath = self.remote_profile['work_path']
        self.ssh = None
        self._setup_ssh(self.remote_host,
                        self.remote_port,
                        username=self.remote_uname,
                        password=self.remote_password)

    def ensure_alive(self,
                     max_check = 10,
                     sleep_time = 10):
        count = 1
        while not self._check_alive():
            if count == max_check:
                raise RuntimeError('cannot connect ssh after %d failures at interval %d s' %
                                   (max_check, sleep_time))
            self._setup_ssh(self.remote_host,
                            self.remote_port,
                            username=self.remote_uname,
                            password=self.remote_password)
            count += 1
            time.sleep(sleep)

    def _check_alive(self):
        if self.ssh == None:
            return False
        try :
            transport = self.ssh.get_transport()
            transport.send_ignore()
            return True
        except EOFError:
            return False        

    def _setup_ssh(self,
                   hostname,
                   port, 
                   username = None,
                   password = None):
        self.ssh = paramiko.SSHClient()        
        # ssh_client.load_system_host_keys()        
        self.ssh.set_missing_host_key_policy(paramiko.WarningPolicy)
        self.ssh.connect(hostname, port=port, username=username, password=password)
        assert(self.ssh.get_transport().is_active())
        transport = self.ssh.get_transport()
        transport.set_keepalive(60)

    def get_ssh_client(self) :
        return self.ssh

    def get_session_root(self) :
        return self.remote_workpath

    def close(self) :
        self.ssh.close()


class RemoteContext (object):
    def __init__ (self,
                  ssh_session,
                  local_root,
                  job_uuid=None,
    ) :
        self.local_root = os.path.abspath(local_root)
        if job_uuid:
           self.job_uuid=job_uuid
        else:
           self.job_uuid = str(uuid.uuid4())
        self.remote_root = os.path.join(ssh_session.get_session_root(), self.job_uuid)
        self.ssh_session = ssh_session
        self.ssh = self.ssh_session.get_ssh_client()        
        self.ssh_session.ensure_alive()
        try:
           sftp = self.ssh.open_sftp()        
           sftp.mkdir(self.remote_root)
           sftp.close()
        except: 
           pass

    def get_job_root(self) :
        return self.remote_root
        
    def upload(self,
               job_dirs,
               local_up_files,
               dereference = True) :
        self.ssh_session.ensure_alive()
        cwd = os.getcwd()
        os.chdir(self.local_root) 
        file_list = []
        for ii in job_dirs :
            for jj in local_up_files :
                file_list.append(os.path.join(ii,jj))        
        self._put_files(file_list, dereference = dereference)
        os.chdir(cwd)

    def download(self, 
                 job_dirs,
                 remote_down_files,
                 back_error=False) :
        self.ssh_session.ensure_alive()
        cwd = os.getcwd()
        os.chdir(self.local_root) 
        file_list = []
        for ii in job_dirs :
            for jj in remote_down_files :
                file_list.append(os.path.join(ii,jj))
            if back_error:
               errors=glob(os.path.join(ii,'error*'))
               file_list.extend(errors)
        self._get_files(file_list)
        os.chdir(cwd)
        
    def block_checkcall(self, 
                        cmd) :
        self.ssh_session.ensure_alive()
        stdin, stdout, stderr = self.ssh.exec_command(('cd %s ;' % self.remote_root) + cmd)
        exit_status = stdout.channel.recv_exit_status() 
        if exit_status != 0:
            raise RuntimeError("Get error code %d in calling %s through ssh with job: %s . message:",
                               (exit_status, cmd, self.job_uuid, stderr.read().decode('utf-8')))
        return stdin, stdout, stderr    

    def block_call(self, 
                   cmd) :
        self.ssh_session.ensure_alive()
        stdin, stdout, stderr = self.ssh.exec_command(('cd %s ;' % self.remote_root) + cmd)
        exit_status = stdout.channel.recv_exit_status() 
        return exit_status, stdin, stdout, stderr

    def clean(self) :        
        self.ssh_session.ensure_alive()
        sftp = self.ssh.open_sftp()        
        self._rmtree(sftp, self.remote_root)
        sftp.close()

    def write_file(self, fname, write_str):
        self.ssh_session.ensure_alive()
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, fname), 'w') as fp :
            fp.write(write_str)
        sftp.close()

    def read_file(self, fname):
        self.ssh_session.ensure_alive()
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, fname), 'r') as fp:
            ret = fp.read().decode('utf-8')
        sftp.close()
        return ret

    def check_file_exists(self, fname):
        self.ssh_session.ensure_alive()
        sftp = self.ssh.open_sftp()
        try:
            sftp.stat(os.path.join(self.remote_root, fname)) 
            ret = True
        except IOError:
            ret = False
        sftp.close()
        return ret        
        
    def call(self, cmd):
        stdin, stdout, stderr = self.ssh.exec_command(cmd)
        # stdin, stdout, stderr = self.ssh.exec_command('echo $$; exec ' + cmd)
        # pid = stdout.readline().strip()
        # print(pid)
        return {'stdin':stdin, 'stdout':stdout, 'stderr':stderr}
    
    def check_finish(self, cmd_pipes):
        return cmd_pipes['stdout'].channel.exit_status_ready()
        
    def get_return(self, cmd_pipes):
        if not self.check_finish(cmd_pipes):
            return None, None, None
        else :
            retcode = cmd_pipes['stdout'].channel.recv_exit_status()
            return retcode, cmd_pipes['stdout'], cmd_pipes['stderr']

    def kill(self, cmd_pipes) :
        raise RuntimeError('dose not work! we do not know how to kill proc through paramiko.SSHClient')
        self.block_checkcall('kill -15 %s' % cmd_pipes['pid'])


    def _rmtree(self, sftp, remotepath, level=0, verbose = False):
        for f in sftp.listdir_attr(remotepath):
            rpath = os.path.join(remotepath, f.filename)
            if stat.S_ISDIR(f.st_mode):
                self._rmtree(sftp, rpath, level=(level + 1))
            else:
                rpath = os.path.join(remotepath, f.filename)
                if verbose: dlog.info('removing %s%s' % ('    ' * level, rpath))
                sftp.remove(rpath)
        if verbose: dlog.info('removing %s%s' % ('    ' * level, remotepath))
        sftp.rmdir(remotepath)

    def _put_files(self,
                   files,
                   dereference = True) :
        of = self.job_uuid + '.tgz'
        # local tar
        cwd = os.getcwd()
        os.chdir(self.local_root)
        if os.path.isfile(of) :
            os.remove(of)
        with tarfile.open(of, "w:gz", dereference = dereference) as tar:
            for ii in files :
                tar.add(ii)
        os.chdir(cwd)
        # trans
        from_f = os.path.join(self.local_root, of)
        to_f = os.path.join(self.remote_root, of)
        sftp = self.ssh.open_sftp()
        sftp.put(from_f, to_f)
        # remote extract
        self.block_checkcall('tar xf %s' % of)
        # clean up
        os.remove(from_f)
        sftp.remove(to_f)
        sftp.close()

    def _get_files(self, 
                   files) :
        of = self.job_uuid + '.tgz'
        flist = ""
        for ii in files :
            flist += " " + ii
        # remote tar
        self.block_checkcall('tar czf %s %s' % (of, flist))
        # trans
        from_f = os.path.join(self.remote_root, of)
        to_f = os.path.join(self.local_root, of)
        if os.path.isfile(to_f) :
            os.remove(to_f)
        sftp = self.ssh.open_sftp()
        sftp.get(from_f, to_f)
        # extract
        cwd = os.getcwd()
        os.chdir(self.local_root)
        with tarfile.open(of, "r:gz") as tar:
            tar.extractall()
        os.chdir(cwd)        
        # cleanup
        os.remove(to_f)
        sftp.remove(from_f)

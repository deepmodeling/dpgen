#!/usr/bin/env python3

import os, sys, paramiko, json, uuid, tarfile, time, stat
from enum import Enum

class JobStatus (Enum) :
    unsubmitted = 1
    waiting = 2
    running = 3
    terminated = 4
    finished = 5
    unknow = 100

class SSHSession (object) :
    def __init__ (self, remote_profile) :
        with open(remote_profile) as fp :
            self.remote_profile = json.load(fp)
        self.remote_host = self.remote_profile['hostname']
        self.remote_port = self.remote_profile['port']
        self.remote_uname = self.remote_profile['username']
        self.remote_workpath = self.remote_profile['work_path']
        self.ssh = self._setup_ssh(self.remote_host, self.remote_port, username = self.remote_uname)
        
    def _setup_ssh(self,
                   hostname,
                   port, 
                   username = None,
                   password = None):
        ssh_client = paramiko.SSHClient()        
        ssh_client.load_system_host_keys()
        ssh_client.set_missing_host_key_policy(paramiko.WarningPolicy)
        ssh_client.connect(hostname, port=port, username=username, password=password)
        return ssh_client

    def get_ssh_client(self) :
        return self.ssh

    def get_session_root(self) :
        return self.remote_workpath

    def close(self) :
        self.ssh.close()


class RemoteJob (object):
    def __init__ (self,
                  ssh_session,
                  job_dir
    ) :
        self.local_dir = os.path.abspath(job_dir)
        self.job_uuid = str(uuid.uuid4())
        # self.job_uuid = 'a21d0017-c9f1-4d29-9a03-97df06965cef'
        self.remote_workpath = os.path.join(ssh_session.get_session_root(), self.job_uuid)
        self.ssh = ssh_session.get_ssh_client()        
        self.sftp = self.ssh.open_sftp()        
        self.sftp.mkdir(self.remote_workpath)
        open('job_uuid', 'w').write(self.job_uuid)

    def upload(self,
               local_up_files) :
        self._put_files(local_up_files)

    def download(self, 
                 remote_down_files) :
        self._get_files(remote_down_files)
        
    def block_call(self, 
                   cmd) :
        stdin, stdout, stderr = self.ssh.exec_command(('cd %s ;' % self.remote_workpath) + cmd)
        exit_status = stdout.channel.recv_exit_status() 
        if exit_status != 0:
            raise RuntimeError("Get error code %d in calling through ssh with job: %s ", (exit_status, self.job_uuid))
        return stdin, stdout, stderr    

    def clean(self) :        
        self._rmtree(self.remote_workpath)

    def _rmtree(self, remotepath, level=0, verbose = False):
        sftp = self.sftp
        for f in sftp.listdir_attr(remotepath):
            rpath = os.path.join(remotepath, f.filename)
            if stat.S_ISDIR(f.st_mode):
                self._rmtree(rpath, level=(level + 1))
            else:
                rpath = os.path.join(remotepath, f.filename)
                if verbose: print('removing %s%s' % ('    ' * level, rpath))
                sftp.remove(rpath)
        if verbose: print('removing %s%s' % ('    ' * level, remotepath))
        sftp.rmdir(remotepath)

    def _put_files(self,
                   files) :
        of = self.job_uuid + '.tgz'
        # local tar
        cwd = os.getcwd()
        os.chdir(self.local_dir)
        if os.path.isfile(of) :
            os.remove(of)
        with tarfile.open(of, "w:gz") as tar:
            for ii in files :
                tar.add(ii)
        os.chdir(cwd)
        # trans
        from_f = os.path.join(self.local_dir, of)
        to_f = os.path.join(self.remote_workpath, of)
        self.sftp.put(from_f, to_f)
        # remote extract
        self.block_call('tar xf %s' % of)
        # clean up
        os.remove(from_f)
        self.sftp.remove(to_f)

    def _get_files(self, 
                   files) :
        of = self.job_uuid + '.tgz'
        flist = ""
        for ii in files :
            flist += " " + ii
        # remote tar
        self.block_call('tar czf %s %s' % (of, flist))
        # trans
        from_f = os.path.join(self.remote_workpath, of)
        to_f = os.path.join(self.local_dir, of)
        if os.path.isfile(to_f) :
            os.remove(to_f)
        self.sftp.get(from_f, to_f)
        # extract
        cwd = os.getcwd()
        os.chdir(self.local_dir)
        with tarfile.open(of, "r:gz") as tar:
            tar.extractall()
        os.chdir(cwd)        
        # cleanup
        os.remove(to_f)
        self.sftp.remove(from_f)

class CloudMachineJob (RemoteJob) :
    def submit(self, 
               cmd) :
        print(self.remote_workpath)
        self.stdin, self.stdout, self.stderr = self.ssh.exec_command(('cd %s ;' % self.remote_workpath) + cmd)

    def check_status(self) :
        if not self._check_finish(self.stdout) :
            return JobStatus.running
        elif self._get_exit_status(self.stdout) == 0 :
            return JobStatus.finished
        else :
            return JobStatus.terminated

    def _check_finish(self, stdout) :
        return stdout.channel.exit_status_ready()
    
    def _get_exit_status(self, stdout) :
        return stdout.channel.recv_exit_status() 
    



ssh_session = SSHSession('localhost.json')        
rjob = CloudMachineJob(ssh_session, '.')
# can upload dirs and normal files
rjob.upload(['batch_exec.py', 'test'])
rjob.block_call('touch a')
rjob.submit('sleep 2')
while rjob.check_status() == JobStatus.running :
    time.sleep(2)
print(rjob.check_status())
# can download dirs and normal files
rjob.download(['a'])
rjob.clean()

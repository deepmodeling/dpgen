#!/usr/bin/env python
# coding: utf-8

import os, sys, paramiko, json, uuid, tarfile, time, stat, shutil
from glob import glob
from enum import Enum
from dpgen import dlog


class JobStatus (Enum) :
    unsubmitted = 1
    waiting = 2
    running = 3
    terminated = 4
    finished = 5
    unknown = 100

class awsMachineJob(object):
    def __init__ (self,
                  remote_root,
                  work_path,
                  job_uuid=None,
    ) :
        self.remote_root=os.path.join(remote_root,work_path)
        self.local_root = os.path.abspath(work_path)
        if job_uuid:
            self.job_uuid=job_uuid
        else:
            self.job_uuid = str(uuid.uuid4())
        
        dlog.info("local_root is %s"% self.local_root)
        dlog.info("remote_root is %s"% self.remote_root)
        
    def upload(self,
               job_dir,
               local_up_files,
               dereference = True) :
        cwd = os.getcwd()
        print('cwd=',cwd)
        os.chdir(self.local_root)       
        for ii in local_up_files :
            print('self.local_root=',self.local_root,'remote_root=',self.remote_root,'job_dir=',job_dir,'ii=',ii)
            if os.path.isfile(os.path.join(job_dir,ii)):
                if not os.path.exists(os.path.join(self.remote_root,job_dir)):
                    os.makedirs(os.path.join(self.remote_root,job_dir))
                shutil.copyfile(os.path.join(job_dir,ii),os.path.join(self.remote_root,job_dir,ii))
            elif os.path.isdir(os.path.join(job_dir,ii)):
                shutil.copytree(os.path.join(job_dir,ii),os.path.join(self.remote_root,job_dir,ii))
            else:
                print('unknownfile','local_root=',self.local_root,'job_dir=',job_dir,'filename=',ii)
        os.chdir(cwd)
    def download(self,
               job_dir,
               remote_down_files,
               dereference = True) :
        for ii in remote_down_files:
         #   print('self.local_root=',self.local_root,'remote_root=',self.remote_root,'job_dir=',job_dir,'ii=',ii)
            file_succ_copy_flag=False
            while not file_succ_copy_flag:
                if os.path.isfile(os.path.join(self.remote_root,job_dir,ii)):
                    shutil.copyfile(os.path.join(self.remote_root,job_dir,ii),os.path.join(self.local_root,job_dir,ii))
                    file_succ_copy_flag=True
                elif os.path.isdir(os.path.join(self.remote_root,job_dir,ii)):
                    try:
                        os.rmdir(os.path.join(self.local_root,job_dir,ii))
                    except Exception:
                        print('dir is not empty   '+str(os.path.join(self.local_root,job_dir,ii)))
                    else:
                        shutil.copytree(os.path.join(self.remote_root,job_dir,ii),os.path.join(self.local_root,job_dir,ii))
                        file_succ_copy_flag=True
                else:
                    print('unknownfile,maybe need for waiting for a while','local_root=',self.local_root,'job_dir=',job_dir,'filename=',ii)
                    time.sleep(5)

def _default_item(resources, key, value) :
    if key not in resources :
        resources[key] = value

def _set_default_resource(res) :
    if res == None :
        res = {}
    _default_item(res, 'numb_node', 1)
    _default_item(res, 'task_per_node', 1)
    _default_item(res, 'numb_gpu', 0)
    _default_item(res, 'time_limit', '1:0:0')
    _default_item(res, 'mem_limit', -1)
    _default_item(res, 'partition', '')
    _default_item(res, 'account', '')
    _default_item(res, 'qos', '')
    _default_item(res, 'constraint_list', [])
    _default_item(res, 'license_list', [])
    _default_item(res, 'exclude_list', [])
    _default_item(res, 'module_unload_list', [])
    _default_item(res, 'module_list', [])
    _default_item(res, 'source_list', [])
    _default_item(res, 'envs', None)
    _default_item(res, 'with_mpi', False)


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
        self.local_key_filename = None
        if 'key_filename' in self.remote_profile:
            self.local_key_filename = self.remote_profile['key_filename']
        self.remote_timeout = None
        if 'timeout' in self.remote_profile:
            self.remote_timeout = self.remote_profile['timeout']
        self.local_key_passphrase = None
        if 'passphrase' in self.remote_profile:
            self.local_key_passphrase = self.remote_profile['passphrase']
        self.remote_workpath = self.remote_profile['work_path']
        self.ssh = self._setup_ssh(hostname=self.remote_host,
                                   port=self.remote_port,
                                   username=self.remote_uname,
                                   password=self.remote_password,
                                   key_filename=self.local_key_filename, 
                                   timeout=self.remote_timeout,
                                   passphrase=self.local_key_passphrase)
        
    def _setup_ssh(self,
                   hostname,
                   port=22,
                   username=None,
                   password=None,
                   key_filename=None,
                   timeout=None,
                   passphrase=None
                   ):
        ssh_client = paramiko.SSHClient()
        ssh_client.load_system_host_keys()
        ssh_client.set_missing_host_key_policy(paramiko.WarningPolicy)
        ssh_client.connect(hostname, port, username, password,
                           key_filename, timeout, passphrase)
        assert(ssh_client.get_transport().is_active())
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
                  local_root,
                  job_uuid=None,
    ) :
        self.local_root = os.path.abspath(local_root)
        if job_uuid:
           self.job_uuid=job_uuid
        else:
           self.job_uuid = str(uuid.uuid4())

        self.remote_root = os.path.join(ssh_session.get_session_root(), self.job_uuid)
        dlog.info("local_root is %s"% local_root)
        dlog.info("remote_root is %s"% self.remote_root)
        self.ssh = ssh_session.get_ssh_client()        
        # keep ssh alive
        transport = self.ssh.get_transport()
        transport.set_keepalive(60)
        try:
           sftp = self.ssh.open_sftp()        
           sftp.mkdir(self.remote_root)
           sftp.close()
        except Exception: 
           pass
        # open('job_uuid', 'w').write(self.job_uuid)

    def get_job_root(self) :
        return self.remote_root
        
    def upload(self,
               job_dirs,
               local_up_files,
               dereference = True) :
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
        stdin, stdout, stderr = self.ssh.exec_command(('cd %s ;' % self.remote_root) + cmd)
        exit_status = stdout.channel.recv_exit_status() 
        if exit_status != 0:
            dlog.info("Error info: %s "%(stderr.readlines()[0]))
            raise RuntimeError("Get error code %d in calling %s through ssh with job: %s "% (exit_status, cmd, self.job_uuid))
        return stdin, stdout, stderr    

    def block_call(self, 
                   cmd) :
        stdin, stdout, stderr = self.ssh.exec_command(('cd %s ;' % self.remote_root) + cmd)
        exit_status = stdout.channel.recv_exit_status() 
        return exit_status, stdin, stdout, stderr

    def clean(self) :        
        sftp = self.ssh.open_sftp()        
        self._rmtree(sftp, self.remote_root)
        sftp.close()

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
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar)
        os.chdir(cwd)        
        # cleanup
        os.remove(to_f)
        sftp.remove(from_f)

class CloudMachineJob (RemoteJob) :
    def submit(self, 
               job_dirs,
               cmd, 
               args = None, 
               resources = None) :
        
        #dlog.info("Current path is",os.getcwd())

        #for ii in job_dirs :
        #    if not os.path.isdir(ii) :
        #        raise RuntimeError("cannot find dir %s" % ii)
        # dlog.info(self.remote_root)
        script_name = self._make_script(job_dirs, cmd, args, resources)
        self.stdin, self.stdout, self.stderr = self.ssh.exec_command(('cd %s; bash %s' % (self.remote_root, script_name)))
        # dlog.info(self.stderr.read().decode('utf-8'))
        # dlog.info(self.stdout.read().decode('utf-8'))

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
    
    def _make_script(self, 
                     job_dirs,
                     cmd, 
                     args = None, 
                     resources = None) :
        _set_default_resource(resources)
        envs = resources['envs']
        module_list = resources['module_list']
        module_unload_list = resources['module_unload_list']
        task_per_node = resources['task_per_node']

        script_name = 'run.sh'
        if args == None :
            args = []
            for ii in job_dirs:
                args.append('')
        script = os.path.join(self.remote_root, script_name)
        sftp = self.ssh.open_sftp()
        with sftp.open(script, 'w') as fp :
            fp.write('#!/bin/bash\n\n')
            # fp.write('set -euo pipefail\n')
            if envs != None :
                for key in envs.keys() :
                    fp.write('export %s=%s\n' % (key, envs[key]))
                fp.write('\n')
            if module_unload_list is not None :
                for ii in module_unload_list :
                    fp.write('module unload %s\n' % ii)
                fp.write('\n')
            if module_list is not None :
                for ii in module_list :
                    fp.write('module load %s\n' % ii)
                fp.write('\n')
            for ii,jj in zip(job_dirs, args) :
                fp.write('cd %s\n' % ii)                
                fp.write('test $? -ne 0 && exit\n')
                if resources['with_mpi'] == True :
                    fp.write('mpirun -n %d %s %s\n' 
                             % (task_per_node, cmd, jj))
                else :
                    fp.write('%s %s\n' % (cmd, jj))
                if 'allow_failure' not in resources or resources['allow_failure'] is False:
                    fp.write('test $? -ne 0 && exit\n')
                fp.write('cd %s\n' % self.remote_root)         
                fp.write('test $? -ne 0 && exit\n')  
            fp.write('\ntouch tag_finished\n')
        sftp.close()
        return script_name


class SlurmJob (RemoteJob) :
    def submit(self, 
               job_dirs,
               cmd,
               args = None, 
               resources = None,
               restart=False) :

        def _submit():
           script_name = self._make_script(job_dirs, cmd, args, res = resources)
           stdin, stdout, stderr = self.block_checkcall(('cd %s; sbatch %s' % (self.remote_root, script_name)))
           subret = (stdout.readlines())
           job_id = subret[0].split()[-1]
           sftp = self.ssh.open_sftp()
        
           with sftp.open(os.path.join(self.remote_root, 'job_id'), 'w') as fp:
             fp.write(job_id)
           sftp.close()

        dlog.debug(restart)
        if restart:
           try:
               status = self.check_status()
               dlog.debug(status)
               if status in [  JobStatus.unsubmitted, JobStatus.unknown, JobStatus.terminated ]:
                  dlog.debug('task restart point !!!')
                  _submit()    
               elif status==JobStatus.waiting:
                  dlog.debug('task is waiting')
               elif status==JobStatus.running:
                  dlog.debug('task is running')
               else:
                  dlog.debug('task is finished')
                 
           except Exception:
               dlog.debug('no job_id file')
               dlog.debug('task restart point !!!')
               _submit()
        else:
           dlog.debug('new task!!!')
           _submit()
   
    def check_status(self) :
        job_id = self._get_job_id()
        if job_id == "" :
            raise RuntimeError("job %s has not been submitted" % self.remote_root)
        ret, stdin, stdout, stderr\
            = self.block_call ("squeue --job " + job_id)
        err_str = stderr.read().decode('utf-8')
        if (ret != 0) :
            if str("Invalid job id specified") in err_str :
                if self._check_finish_tag() :
                    return JobStatus.finished
                else :
                    return JobStatus.terminated
            else :
                raise RuntimeError\
                    ("status command squeue fails to execute\nerror message:%s\nreturn code %d\n" % (err_str, ret))
        status_line = stdout.read().decode('utf-8').split ('\n')[-2]
        status_word = status_line.split ()[-4]
        if      status_word in ["PD","CF","S"] :
            return JobStatus.waiting
        elif    status_word in ["R","CG"] :
            return JobStatus.running
        elif    status_word in ["C","E","K","BF","CA","CD","F","NF","PR","SE","ST","TO"] :
            if self._check_finish_tag() :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown
    
    def _get_job_id(self) :
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, 'job_id'), 'r') as fp:            
            ret = fp.read().decode('utf-8')
        sftp.close()
        return ret

    def _check_finish_tag(self) :
        sftp = self.ssh.open_sftp()
        try:
            sftp.stat(os.path.join(self.remote_root, 'tag_finished')) 
            ret =  True
        except IOError:
            ret = False
        sftp.close()
        return ret

    def _make_squeue(self,mdata1, res):
        ret = ''
        ret += 'squeue -u %s ' % mdata1['username']
        ret += '-p %s ' % res['partition']
        ret += '| grep PD'
        return ret

    def _make_script(self, 
                     job_dirs,
                     cmd,
                     args = None, 
                     res = None) :
        _set_default_resource(res)
        ret = ''
        ret += "#!/bin/bash -l\n"
        ret += "#SBATCH -N %d\n" % res['numb_node']
        ret += "#SBATCH --ntasks-per-node %d\n" % res['task_per_node']
        ret += "#SBATCH -t %s\n" % res['time_limit']
        if res['mem_limit'] > 0 :
            ret += "#SBATCH --mem %dG \n" % res['mem_limit']
        if len(res['account']) > 0 :
            ret += "#SBATCH --account %s \n" % res['account']
        if len(res['partition']) > 0 :
            ret += "#SBATCH --partition %s \n" % res['partition']
        if len(res['qos']) > 0 :
            ret += "#SBATCH --qos %s \n" % res['qos']
        if res['numb_gpu'] > 0 :
            ret += "#SBATCH --gres=gpu:%d\n" % res['numb_gpu']
        for ii in res['constraint_list'] :
            ret += '#SBATCH -C %s \n' % ii
        for ii in res['license_list'] :
            ret += '#SBATCH -L %s \n' % ii
        if len(res['exclude_list']) >0:
            temp_exclude = ""
            for ii in res['exclude_list'] :
                temp_exclude += ii
                temp_exclude += ","
            temp_exclude = temp_exclude[:-1]
            ret += '#SBATCH --exclude %s \n' % temp_exclude
        ret += "\n"
        # ret += 'set -euo pipefail\n\n'
        for ii in res['module_unload_list'] :
            ret += "module unload %s\n" % ii
        for ii in res['module_list'] :
            ret += "module load %s\n" % ii
        ret += "\n"
        for ii in res['source_list'] :
            ret += "source %s\n" %ii
        ret += "\n"
        envs = res['envs']
        if envs != None :
            for key in envs.keys() :
                ret += 'export %s=%s\n' % (key, envs[key])
            ret += '\n'

        if args == None :
            args = []
            for ii in job_dirs:
                args.append('')

        try:
           cvasp=res['cvasp']
           try:
              fp_max_errors = res['fp_max_errors']
           except Exception:
              fp_max_errors = 3
        except Exception:
           cvasp=False

        for ii,jj in zip(job_dirs, args) :
            ret += 'cd %s\n' % ii
            ret += 'test $? -ne 0 && exit\n\n'

            if cvasp:
                cmd=cmd.split('1>')[0].strip()
                if res['with_mpi'] :
                    ret += 'if [ -f tag_finished ] ;then\n'
                    ret += '  echo gogogo \n'
                    ret += 'else\n'
                    ret += '  python ../cvasp.py "srun %s" %s %s 1>log 2>log\n' % (cmd, fp_max_errors, jj)
                    ret += '  if test $? -ne 0 \n'
                    ret += '  then\n'
                    ret += '     exit\n'
                    ret += '  else\n'
                    ret += '     touch tag_finished\n'
                    ret += '  fi\n'
                    ret += 'fi\n\n'
                else :
                    ret += 'if [ -f tag_finished ] ;then\n'
                    ret += '  echo gogogo \n'
                    ret += 'else\n'
                    ret += '  python ../cvasp.py "%s" %s %s 1>log 2>log\n' % (cmd, fp_max_errors, jj)
                    ret += '  if test $? -ne 0 \n'
                    ret += '  then\n'
                    ret += '     exit\n'
                    ret += '  else\n'
                    ret += '     touch tag_finished\n'
                    ret += '  fi\n'
                    ret += 'fi\n\n'
            else:
                if res['with_mpi'] :
                    ret += 'if [ -f tag_finished ] ;then\n'
                    ret += '  echo gogogo \n'
                    ret += 'else\n'
                    ret += '  srun %s %s\n' % (cmd, jj)
                    ret += '  if test $? -ne 0 \n'
                    ret += '  then\n'
                    ret += '     exit\n'
                    ret += '  else\n'
                    ret += '     touch tag_finished\n'
                    ret += '  fi\n'
                    ret += 'fi\n\n'
                else :
                    ret += 'if [ -f tag_finished ] ;then\n'
                    ret += '  echo gogogo \n'
                    ret += 'else\n'
                    ret += '  %s %s\n' % (cmd, jj)
                    ret += '  if test $? -ne 0 \n'
                    ret += '  then\n'
                    ret += '     exit\n'
                    ret += '  else\n'
                    ret += '     touch tag_finished\n'
                    ret += '  fi\n'
                    ret += 'fi\n\n'
            if 'allow_failure' not in res or res['allow_failure'] is False:
                ret += 'test $? -ne 0 && exit\n'
            ret += 'cd %s\n' % self.remote_root
            ret += 'test $? -ne 0 && exit\n'
        ret += '\ntouch tag_finished\n'

        script_name = 'run.sub'
        script = os.path.join(self.remote_root, script_name)
        sftp = self.ssh.open_sftp()
        with sftp.open(script, 'w') as fp :
            fp.write(ret)
        sftp.close()

        return script_name


class PBSJob (RemoteJob) :
    def submit(self, 
               job_dirs,
               cmd,
               args = None, 
               resources = None) :
        script_name = self._make_script(job_dirs, cmd, args, res = resources)
        stdin, stdout, stderr = self.block_checkcall(('cd %s; qsub %s' % (self.remote_root, script_name)))
        subret = (stdout.readlines())
        job_id = subret[0].split()[0]
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, 'job_id'), 'w') as fp:
            fp.write(job_id)
        sftp.close()

    def check_status(self) :
        job_id = self._get_job_id()
        if job_id == "" :
            raise RuntimeError("job %s is has not been submitted" % self.remote_root)
        ret, stdin, stdout, stderr\
            = self.block_call ("qstat " + job_id)
        err_str = stderr.read().decode('utf-8')
        if (ret != 0) :
            if str("qstat: Unknown Job Id") in err_str :
                if self._check_finish_tag() :
                    return JobStatus.finished
                else :
                    return JobStatus.terminated
            else :
                raise RuntimeError ("status command qstat fails to execute. erro info: %s return code %d"
                                    % (err_str, ret))
        status_line = stdout.read().decode('utf-8').split ('\n')[-2]
        status_word = status_line.split ()[-2]        
        # dlog.info (status_word)
        if      status_word in ["Q","H"] :
            return JobStatus.waiting
        elif    status_word in ["R"] :
            return JobStatus.running
        elif    status_word in ["C","E","K"] :
            if self._check_finish_tag() :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown
    
    def _get_job_id(self) :
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, 'job_id'), 'r') as fp:
            ret = fp.read().decode('utf-8')
        sftp.close()
        return ret

    def _check_finish_tag(self) :
        sftp = self.ssh.open_sftp()
        try:
            sftp.stat(os.path.join(self.remote_root, 'tag_finished')) 
            ret =  True
        except IOError:
            ret = False
        sftp.close()
        return ret

    def _make_script(self, 
                     job_dirs,
                     cmd,
                     args = None, 
                     res = None) :
        _set_default_resource(res)
        ret = ''
        ret += "#!/bin/bash -l\n"
        if res['numb_gpu'] == 0:
            ret += '#PBS -l nodes=%d:ppn=%d\n' % (res['numb_node'], res['task_per_node'])
        else :
            ret += '#PBS -l nodes=%d:ppn=%d:gpus=%d\n' % (res['numb_node'], res['task_per_node'], res['numb_gpu'])
        ret += '#PBS -l walltime=%s\n' % (res['time_limit'])
        if res['mem_limit'] > 0 :
            ret += "#PBS -l mem=%dG \n" % res['mem_limit']
        ret += '#PBS -j oe\n'
        if len(res['partition']) > 0 :
            ret += '#PBS -q %s\n' % res['partition']
        ret += "\n"
        for ii in res['module_unload_list'] :
            ret += "module unload %s\n" % ii
        for ii in res['module_list'] :
            ret += "module load %s\n" % ii
        ret += "\n"
        for ii in res['source_list'] :
            ret += "source %s\n" %ii
        ret += "\n"
        envs = res['envs']
        if envs != None :
            for key in envs.keys() :
                ret += 'export %s=%s\n' % (key, envs[key])
            ret += '\n'
        ret += 'cd $PBS_O_WORKDIR\n\n'

        if args == None :
            args = []
            for ii in job_dirs:
                args.append('')
        for ii,jj in zip(job_dirs, args) :
            ret += 'cd %s\n' % ii
            ret += 'test $? -ne 0 && exit\n'
            if res['with_mpi'] :
                ret += 'mpirun -machinefile $PBS_NODEFILE -n %d %s %s\n' % (res['numb_node'] * res['task_per_node'], cmd, jj)
            else :
                ret += '%s %s\n' % (cmd, jj)                
            if 'allow_failure' not in res or res['allow_failure'] is False:
                ret += 'test $? -ne 0 && exit\n'
            ret += 'cd %s\n' % self.remote_root
            ret += 'test $? -ne 0 && exit\n'
        ret += '\ntouch tag_finished\n'

        script_name = 'run.sub'
        script = os.path.join(self.remote_root, script_name)
        sftp = self.ssh.open_sftp()
        with sftp.open(script, 'w') as fp :
            fp.write(ret)
        sftp.close()

        return script_name


# ssh_session = SSHSession('localhost.json')        
# rjob = CloudMachineJob(ssh_session, '.')
# # can upload dirs and normal files
# rjob.upload(['job0', 'job1'], ['batch_exec.py', 'test'])
# rjob.submit(['job0', 'job1'], 'touch a; sleep 2')
# while rjob.check_status() == JobStatus.running :
#     dlog.info('checked')
#     time.sleep(2)
# dlog.info(rjob.check_status())
# # can download dirs and normal files
# rjob.download(['job0', 'job1'], ['a'])
# # rjob.clean()


class LSFJob (RemoteJob) :
    def submit(self,
               job_dirs,
               cmd,
               args = None,
               resources = None,
               restart = False):
        dlog.debug(restart)
        if restart:
            status = self.check_status()
            if status in [  JobStatus.unsubmitted, JobStatus.unknown, JobStatus.terminated ]:
                dlog.debug('task restart point !!!')
                if 'task_max' in resources and resources['task_max'] > 0:
                    while self.check_limit(task_max=resources['task_max']):
                        time.sleep(60)
                self._submit(job_dirs, cmd, args, resources)
            elif status==JobStatus.waiting:
                dlog.debug('task is waiting')
            elif status==JobStatus.running:
                dlog.debug('task is running')
            else:
                dlog.debug('task is finished')
            #except Exception:
               #dlog.debug('no job_id file')
               #dlog.debug('task restart point !!!')
               #self._submit(job_dirs, cmd, args, resources)
        else:
            dlog.debug('new task!!!')
            if 'task_max' in resources and resources['task_max'] > 0:
                while self.check_limit(task_max=resources['task_max']):
                    time.sleep(60)
            self._submit(job_dirs, cmd, args, resources)
        if resources.get('wait_time', False):
            time.sleep(resources['wait_time']) # For preventing the crash of the tasks while submitting.

    def _submit(self, 
               job_dirs,
               cmd,
               args = None, 
               resources = None) :
        script_name = self._make_script(job_dirs, cmd, args, res = resources)
        stdin, stdout, stderr = self.block_checkcall(('cd %s; bsub < %s' % (self.remote_root, script_name)))
        subret = (stdout.readlines())
        job_id = subret[0].split()[1][1:-1]
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, 'job_id'), 'w') as fp:
            fp.write(job_id)
        sftp.close()

    def check_limit(self, task_max):
        stdin_run, stdout_run, stderr_run = self.block_checkcall("bjobs | grep RUN | wc -l")
        njobs_run = int(stdout_run.read().decode('utf-8').split ('\n')[0])
        stdin_pend, stdout_pend, stderr_pend = self.block_checkcall("bjobs | grep PEND | wc -l")
        njobs_pend = int(stdout_pend.read().decode('utf-8').split ('\n')[0])
        if (njobs_pend + njobs_run) < task_max:
            return False
        else:
            return True

    def check_status(self) :
        try:
            job_id = self._get_job_id()
        except Exception:
            return JobStatus.terminated
        if job_id == "" :
            raise RuntimeError("job %s is has not been submitted" % self.remote_root)
        ret, stdin, stdout, stderr\
            = self.block_call ("bjobs " + job_id)
        err_str = stderr.read().decode('utf-8')
        if ("Job <%s> is not found" % job_id) in err_str :
            if self._check_finish_tag() :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        elif ret != 0 :
            raise RuntimeError ("status command bjobs fails to execute. erro info: %s return code %d"
                                    % (err_str, ret))
        status_out = stdout.read().decode('utf-8').split('\n')
        if len(status_out) < 2:
            return JobStatus.unknown
        else:
            status_line = status_out[1]
            status_word = status_line.split()[2]

        # ref: https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bjobs.1.html
        if status_word in ["PEND", "WAIT", "PSUSP"] :
            return JobStatus.waiting
        elif status_word in ["RUN", "USUSP"] :
            return JobStatus.running
        elif status_word in ["DONE","EXIT"] :
            if self._check_finish_tag() :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown
    
    def _get_job_id(self) :
        sftp = self.ssh.open_sftp()
        with sftp.open(os.path.join(self.remote_root, 'job_id'), 'r') as fp:
            ret = fp.read().decode('utf-8')
        sftp.close()
        return ret

    def _check_finish_tag(self) :
        sftp = self.ssh.open_sftp()
        try:
            sftp.stat(os.path.join(self.remote_root, 'tag_finished')) 
            ret =  True
        except IOError:
            ret = False
        sftp.close()
        return ret

    def _make_script(self, 
                     job_dirs,
                     cmd,
                     args = None, 
                     res = None) :
        _set_default_resource(res)
        ret = ''
        ret += "#!/bin/bash -l\n#BSUB -e %J.err\n#BSUB -o %J.out\n"
        if res['numb_gpu'] == 0:
            ret += '#BSUB -n %d\n#BSUB -R span[ptile=%d]\n' % (
                res['numb_node'] * res['task_per_node'], res['node_cpu'])
        else:
            if res['node_cpu']:
                ret += '#BSUB -R span[ptile=%d]\n' % res['node_cpu']
            if 'new_lsf_gpu' in res and res['new_lsf_gpu'] == True:
                # supportted in LSF >= 10.1.0 SP6
                # ref: https://www.ibm.com/support/knowledgecenter/en/SSWRJV_10.1.0/lsf_resource_sharing/use_gpu_res_reqs.html
                ret += '#BSUB -n %d\n#BSUB -gpu "num=%d:mode=shared:j_exclusive=yes"\n' % (
                    res['numb_gpu'], res['task_per_node'])
            else:
                ret += '#BSUB -n %d\n#BSUB -R "select[ngpus >0] rusage[ngpus_excl_p=%d]"\n' % (
                    res['numb_gpu'], res['task_per_node'])
        if res['time_limit']:
            ret += '#BSUB -W %s\n' % (res['time_limit'].split(':')[
                0] + ':' + res['time_limit'].split(':')[1])
        if res['mem_limit'] > 0 :
            ret += "#BSUB -M %d \n" % (res['mem_limit'])
        ret += '#BSUB -J %s\n' % (res['job_name'] if 'job_name' in res else 'dpgen')
        if len(res['partition']) > 0 :
            ret += '#BSUB -q %s\n' % res['partition']
        ret += "\n"
        for ii in res['module_unload_list'] :
            ret += "module unload %s\n" % ii
        for ii in res['module_list'] :
            ret += "module load %s\n" % ii
        ret += "\n"
        for ii in res['source_list'] :
            ret += "source %s\n" %ii
        ret += "\n"
        envs = res['envs']
        if envs != None :
            for key in envs.keys() :
                ret += 'export %s=%s\n' % (key, envs[key])
            ret += '\n'

        if args == None :
            args = []
            for ii in job_dirs:
                args.append('')
        for ii,jj in zip(job_dirs, args) :
            ret += 'cd %s\n' % ii
            ret += 'test $? -ne 0 && exit\n'
            if res['with_mpi']:
                ret += 'mpirun -machinefile $LSB_DJOB_HOSTFILE -n %d %s %s\n' % (
                    res['numb_node'] * res['task_per_node'], cmd, jj)
            else :
                ret += '%s %s\n' % (cmd, jj)                
            if 'allow_failure' not in res or res['allow_failure'] is False:
                ret += 'test $? -ne 0 && exit\n'
            ret += 'cd %s\n' % self.remote_root
            ret += 'test $? -ne 0 && exit\n'
        ret += '\ntouch tag_finished\n'

        script_name = 'run.sub'
        script = os.path.join(self.remote_root, script_name)
        sftp = self.ssh.open_sftp()
        with sftp.open(script, 'w') as fp :
            fp.write(ret)
        sftp.close()

        return script_name

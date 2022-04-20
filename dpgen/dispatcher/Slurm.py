import os,getpass,time
from dpgen.dispatcher.Batch import Batch
from dpgen.dispatcher.JobStatus import JobStatus

def _default_item(resources, key, value) :
    if key not in resources :
        resources[key] = value

class Slurm(Batch) :

    def check_status(self) :
        """
        check the status of a job
        """
        job_id = self._get_job_id()
        if job_id == '' :
            return JobStatus.unsubmitted
        while True:
            stat = self._check_status_inner(job_id)
            if stat != JobStatus.completing:
                return stat
            else:
                time.sleep(5)

    def do_submit(self, 
                  job_dirs,
                  cmd,
                  args = None, 
                  res = None,
                  outlog = 'log',
                  errlog = 'err'):
        if res == None:
            res = self.default_resources(res)
        if 'task_max' in res and res['task_max'] > 0:
            while self._check_sub_limit(task_max=res['task_max']):
                time.sleep(60)
        script_str = self.sub_script(job_dirs, cmd, args=args, res=res, outlog=outlog, errlog=errlog)
        self.context.write_file(self.sub_script_name, script_str)
        stdin, stdout, stderr = self.context.block_checkcall('cd %s && %s %s' % (self.context.remote_root, 'sbatch', self.sub_script_name))
        subret = (stdout.readlines())
        job_id = subret[0].split()[-1]
        self.context.write_file(self.job_id_name, job_id)        
                
    def default_resources(self, res_) :
        """
        set default value if a key in res_ is not fhound
        """
        if res_ == None :
            res = {}
        else:
            res = res_
        _default_item(res, 'numb_node', 1)
        _default_item(res, 'task_per_node', 1)
        _default_item(res, 'cpus_per_task', -1)
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
        _default_item(res, 'cuda_multi_tasks', False)
        _default_item(res, 'allow_failure', False)
        _default_item(res, 'cvasp', False)
        return res

    def sub_script_head(self, res):
        ret = ''
        ret += "#!/bin/bash -l\n"
        ret += "#SBATCH -N %d\n" % res['numb_node']
        ret += "#SBATCH --ntasks-per-node=%d\n" % res['task_per_node']
        if res['cpus_per_task'] > 0 :            
            ret += "#SBATCH --cpus-per-task=%d\n" % res['cpus_per_task']
        ret += "#SBATCH -t %s\n" % res['time_limit']
        if res['mem_limit'] > 0 :
            ret += "#SBATCH --mem=%dG \n" % res['mem_limit']
        if 'job_name' in res:
            if len(res['job_name']) > 0:
                ret += '#SBATCH --job-name=%s\n' % res['job_name']
        if len(res['account']) > 0 :
            ret += "#SBATCH --account=%s \n" % res['account']
        if len(res['partition']) > 0 :
            ret += "#SBATCH --partition=%s \n" % res['partition']
        if len(res['qos']) > 0 :
            ret += "#SBATCH --qos=%s \n" % res['qos']
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
            ret += '#SBATCH --exclude=%s \n' % temp_exclude
        for flag in res.get('custom_flags', []):
            ret += '#SBATCH %s \n' % flag
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
        return ret

    def sub_script_cmd(self,
                       cmd,
                       arg,
                       res) :
        try:
            cvasp=res['cvasp']
            fp_max_errors = 3
            try:
                fp_max_errors = res['fp_max_errors']
            except Exception:
                pass
        except Exception:
            cvasp=False

        _cmd = cmd.split('1>')[0].strip()
        if cvasp :
            if res['with_mpi']:
                _cmd = 'python cvasp.py "srun %s %s" %s' % (_cmd, arg, fp_max_errors)
            else :
                _cmd = 'python cvasp.py "%s %s" %s' % (_cmd, arg, fp_max_errors)
        else :
            if res['with_mpi']:
                _cmd = 'srun %s %s' % (_cmd, arg)
            else :
                _cmd = '%s %s' % (_cmd, arg)        
        return _cmd

    def _get_job_id(self) :
        if self.context.check_file_exists(self.job_id_name) :
            return self.context.read_file(self.job_id_name)
        else:
            return ""

    def _check_status_inner(self, job_id, retry=0):
        ret, stdin, stdout, stderr\
            = self.context.block_call ('squeue -o "%.18i %.2t" -j ' + job_id)
        if (ret != 0) :
            err_str = stderr.read().decode('utf-8')
            if str("Invalid job id specified") in err_str :
                if self.check_finish_tag() :
                    return JobStatus.finished
                else :
                    return JobStatus.terminated
            else :
                # retry 3 times
                if retry < 3:
                    # rest 60s
                    time.sleep(60)
                    return self._check_status_inner(job_id, retry=retry+1)
                raise RuntimeError\
                    ("status command squeue fails to execute\nerror message:%s\nreturn code %d\n" % (err_str, ret))
        status_line = stdout.read().decode('utf-8').split ('\n')[-2]
        status_word = status_line.split ()[-1]
        if not (len(status_line.split()) == 2 and status_word.isupper()): 
            raise RuntimeError("Error in getting job status, " +
                              f"status_line = {status_line}, " + 
                              f"parsed status_word = {status_word}")
        if status_word in ["PD","CF","S"] :
            return JobStatus.waiting
        elif status_word in ["R"] :
            return JobStatus.running
        elif status_word in ["CG"] :
            return JobStatus.completing
        elif status_word in ["C","E","K","BF","CA","CD","F","NF","PR","SE","ST","TO"] :
            if self.check_finish_tag() :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown                    


    def _check_sub_limit(self, task_max, **kwarg) :
        if task_max <= 0:
            return True
        username = getpass.getuser()
        stdin, stdout, stderr = self.context.block_checkcall('squeue -u %s -h' % username)
        nj = len(stdout.readlines())
        return nj >= task_max

    def _make_squeue(self,mdata1, res):
        ret = ''
        ret += 'squeue -u %s ' % mdata1['username']
        ret += '-p %s ' % res['partition']
        ret += '| grep PD'
        return ret

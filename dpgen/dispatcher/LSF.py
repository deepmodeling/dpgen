import os,getpass,time
from dpgen.dispatcher.Batch import Batch
from dpgen.dispatcher.JobStatus import JobStatus

def _default_item(resources, key, value) :
    if key not in resources :
        resources[key] = value

class LSF(Batch) :
    
    def check_status(self):
        try:
            job_id = self._get_job_id()
        except:
            return JobStatus.terminated
        if job_id == "" :
            raise RuntimeError("job %s has not been submitted" % self.context.remote_root)
        ret, stdin, stdout, stderr\
            = self.context.block_call ("bjobs " + job_id)
        err_str = stderr.read().decode('utf-8')
        if ("Job <%s> is not found" % job_id) in err_str :
            if self.check_finish_tag() :
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
        if      status_word in ["PEND", "WAIT", "PSUSP"] :
            return JobStatus.waiting
        elif    status_word in ["RUN", "USUSP"] :
            return JobStatus.running
        elif    status_word in ["DONE","EXIT"] :
            if self.check_finish_tag() :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown


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
        stdin, stdout, stderr = self.context.block_checkcall('cd %s && %s < %s' % (self.context.remote_root, 'bsub', self.sub_script_name))
        subret = (stdout.readlines())
        job_id = subret[0].split()[1][1:-1]
        self.context.write_file(self.job_id_name, job_id)        


    def default_resources(self, res_) :
        """
        set default value if a key in res_ is not fhound
        """
        if res_ == None :
            res = {}
        else:
            res = res_
        _default_item(res, 'node_cpu', 1)
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
        ret += "#!/bin/bash -l\n#BSUB -e %J.err\n#BSUB -o %J.out\n"
        if res['numb_gpu'] == 0:
            ret += '#BSUB -n %d\n#BSUB -R span[ptile=%d]\n' % (
                res['numb_node'] * res['task_per_node'], res['node_cpu'])
        else:
            if res['node_cpu']:
                ret += '#BSUB -R span[ptile=%d]\n' % res['node_cpu']
            if res.get('new_lsf_gpu', False):
                # supported in LSF >= 10.1.0.3
                # ref: https://www.ibm.com/support/knowledgecenter/en/SSWRJV_10.1.0
                # /lsf_resource_sharing/use_gpu_res_reqs.html
                if res.get('exclusive', False):
                    j_exclusive = "no"
                else:
                    j_exclusive = "yes"
                ret += '#BSUB -n %d\n#BSUB -gpu "num=%d:mode=shared:j_exclusive=%s"\n' % (
                    res['task_per_node'], res['numb_gpu'], j_exclusive)
            else:
                ret += '#BSUB -n %d\n#BSUB -R "select[ngpus >0] rusage[ngpus_excl_p=%d]"\n' % (
                    res['task_per_node'], res['numb_gpu'])
        if res['time_limit']:
            ret += '#BSUB -W %s\n' % (res['time_limit'].split(':')[
                0] + ':' + res['time_limit'].split(':')[1])
        if res['mem_limit'] > 0 :
            ret += "#BSUB -M %d \n" % (res['mem_limit'])
        ret += '#BSUB -J %s\n' % (res['job_name'] if 'job_name' in res else 'dpgen')
        if len(res['partition']) > 0 :
            ret += '#BSUB -q %s\n' % res['partition']
        if len(res['exclude_list']) > 0:
            ret += '#BSUB -R "select['
            temp_exclude = []
            for ii in res['exclude_list']:
                temp_exclude.append('hname != %s' % ii)
            ret += ' && '.join(temp_exclude)
            ret += ']"\n'
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
        if res['with_mpi']:
            ret = 'mpirun -machinefile $LSB_DJOB_HOSTFILE -n %d %s %s' % (
                    res['numb_node'] * res['task_per_node'], cmd, arg)
        else :
            ret = '%s %s' % (cmd, arg)
        return ret


    def _get_job_id(self) :
        if self.context.check_file_exists(self.job_id_name) :
            return self.context.read_file(self.job_id_name)
        else:
            return ""


    def _check_sub_limit(self, task_max, **kwarg) :
        stdin_run, stdout_run, stderr_run = self.context.block_checkcall("bjobs | grep RUN | wc -l")
        njobs_run = int(stdout_run.read().decode('utf-8').split ('\n')[0])
        stdin_pend, stdout_pend, stderr_pend = self.context.block_checkcall("bjobs | grep PEND | wc -l")
        njobs_pend = int(stdout_pend.read().decode('utf-8').split ('\n')[0])
        if (njobs_pend + njobs_run) < task_max:
            return False
        else:
            return True


    def _make_squeue(self, mdata1, res):
        ret = ''
        ret += 'bjobs -u %s ' % mdata1['username']
        ret += '-q %s ' % res['partition']
        ret += '| grep PEND '
        return ret

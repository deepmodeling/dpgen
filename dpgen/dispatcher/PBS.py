import os,getpass,time
from dpgen.dispatcher.Batch import Batch
from dpgen.dispatcher.JobStatus import JobStatus

def _default_item(resources, key, value) :
    if key not in resources :
        resources[key] = value

class PBS(Batch) :

    def check_status(self) :
        job_id = self._get_job_id()
        if job_id == "" :
            return JobStatus.unsubmitted
        ret, stdin, stdout, stderr\
            = self.context.block_call ("qstat " + job_id)
        err_str = stderr.read().decode('utf-8')
        if (ret != 0) :
            if str("qstat: Unknown Job Id") in err_str or str("Job has finished") in err_str:
                if self.check_finish_tag() :
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
        # if 'task_max' in res and res['task_max'] > 0:
        #     while self._check_sub_limit(task_max=res['task_max']):
        #         time.sleep(60)
        script_str = self.sub_script(job_dirs, cmd, args=args, res=res, outlog=outlog, errlog=errlog)
        self.context.write_file(self.sub_script_name, script_str)
        stdin, stdout, stderr = self.context.block_checkcall('cd %s && %s %s' % (self.context.remote_root, 'qsub', self.sub_script_name))
        subret = (stdout.readlines())
        job_id = subret[0].split()[0]
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
        _default_item(res, 'allow_failure', True)
        _default_item(res, 'cvasp', False)
        return res

    def sub_script_head(self, res):
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
        return ret

    def sub_script_cmd(self,
                       cmd,
                       arg,
                       res) :
        if res['with_mpi']:
            ret = 'mpirun -machinefile $PBS_NODEFILE -n %d %s %s' % (
                    res['numb_node'] * res['task_per_node'], cmd, arg)
        else :
            ret = '%s %s' % (cmd, arg)
        return ret        

    def _get_job_id(self) :
        if self.context.check_file_exists(self.job_id_name) :
            return self.context.read_file(self.job_id_name)
        else:
            return ""


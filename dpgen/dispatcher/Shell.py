import os,getpass,time
from dpgen.dispatcher.Batch import Batch
from dpgen.dispatcher.JobStatus import JobStatus

def _default_item(resources, key, value) :
    if key not in resources :
        resources[key] = value

def _set_default_resource(res_) :
    if res_ == None :
        res = {}
    else:
        res = res_
    _default_item(res, 'task_per_node', 1)
    _default_item(res, 'module_list', [])
    _default_item(res, 'module_unload_list', [])
    _default_item(res, 'source_list', [])
    _default_item(res, 'envs', {})
    _default_item(res, 'with_mpi', False)
    _default_item(res, 'cuda_multi_tasks', False)
    _default_item(res, 'allow_failure', False)
    _default_item(res, 'cvasp', False)
    return res


class Shell(Batch) :

    def check_status(self) :
        if not self.context.check_finish(self.proc) :
            return JobStatus.running
        elif (self.context.get_return(self.proc))[0] == 0 :
            return JobStatus.finished
        else :
            return JobStatus.terminated

    def do_submit(self, 
                  job_dirs,
                  cmd,
                  args = None, 
                  res = None):
        if res == None:
            res = {}
        script_str = self._sub_script(job_dirs, cmd, args=args, res=res)
        self.context.write_file('run.sub', script_str)
        self.proc = self.context.call('cd %s && exec bash %s' % (self.context.remote_root, 'run.sub'))


    def _script_head(self, resources) :
        envs = resources['envs']
        module_list = resources['module_list']
        module_unload_list = resources['module_unload_list']
        task_per_node = resources['task_per_node']
        source_list = resources['source_list']
        
        ret = ''
        ret += ('#!/bin/bash\n\n')
        # fp.write('set -euo pipefail\n')
        for key in envs.keys() :
            ret += ('export %s=%s\n' % (key, envs[key]))
        ret += ('\n')
        for ii in module_unload_list :
            ret += ('module unload %s\n' % ii)
        ret += ('\n')
        for ii in module_list :
            ret += ('module load %s\n' % ii)
        ret += ('\n')
        for ii in source_list :
            ret += 'source %s\n'
        ret += ('\n')
        return ret


    def _script_cmd(self,
                    res,
                    job_dirs,
                    cmd,
                    args,
                    idx,
                    outlog = 'log',
                    errlog = 'err',
                    cvasp = False,
                    fp_max_errors = 3) :
        ret = ""
        for ii,jj in zip(job_dirs, args[idx]) :
            ret += 'cd %s\n' % ii
            ret += 'test $? -ne 0 && exit\n\n'
            _cmd = cmd[idx].split('1>')[0].strip()
            if cvasp :
                if res['with_mpi']:
                    _cmd = 'python ../cvasp.py "mpirun -n %d %s" %s' % (res['task_per_node'], _cmd, fp_max_errors)
                else :
                    _cmd = 'python ../cvasp.py "%s" %s' % (_cmd, fp_max_errors)
            else :
                if res['with_mpi']:
                    _cmd = 'mpirun -n %d %s ' % (res['task_per_node'],  _cmd)
                else :
                    _cmd = '%s ' % (_cmd)
            _cmd += ' %s 1> %s 2> %s ' % (jj, outlog, errlog)
            ret += 'if [ ! -f tag_%d_finished ] ;then\n' % idx
            ret += '  %s\n' % (_cmd)
            if res['allow_failure'] is False:
                ret += '  if test $? -ne 0; then exit; else touch tag_%d_finished; fi \n' % idx
            ret += 'fi\n\n'
            ret += 'cd %s\n' % self.context.remote_root
            ret += 'test $? -ne 0 && exit\n'
        return ret        
        

    def _sub_script(self, 
                    job_dirs,
                    cmd, 
                    args = None,
                    res  = None,
                    outlog = 'log',
                    errlog = 'err') :
        res = _set_default_resource(res)
        ret = self._script_head(res)
        if not isinstance(cmd, list):
            cmd = [cmd]
        if args == None :
            args = []
            for ii in cmd:
                _args = []
                for jj in job_dirs:
                    _args.append('')
                args.append(_args)
        try:
            cvasp=res['cvasp']
            fp_max_errors = 3
            try:
                fp_max_errors = res['fp_max_errors']
            except:
                pass
        except:
            cvasp=False
        # loop over commands 
        for ii in range(len(cmd)):            
            # for one command
            ret += self._script_cmd(res,
                                    job_dirs,
                                    cmd,
                                    args,
                                    ii,
                                    outlog=outlog,
                                    errlog=errlog,
                                    cvasp=cvasp,
                                    fp_max_errors=fp_max_errors)
        # append finish tag
        ret += '\ntouch tag_finished\n'
        return ret

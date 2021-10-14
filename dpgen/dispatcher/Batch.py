import os,sys,time

from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog


class Batch(object) :
    def __init__ (self,
                  context, 
                  uuid_names = True) :
        self.context = context
        self.uuid_names = uuid_names
        if uuid_names:
            self.upload_tag_name = '%s_tag_upload' % self.context.job_uuid
            self.finish_tag_name = '%s_tag_finished' % self.context.job_uuid
            self.sub_script_name = '%s.sub' % self.context.job_uuid
            self.job_id_name = '%s_job_id' % self.context.job_uuid
        else:
            self.upload_tag_name = 'tag_upload'
            self.finish_tag_name = 'tag_finished'
            self.sub_script_name = 'run.sub'
            self.job_id_name = 'job_id'

    def check_status(self) :
        raise RuntimeError('abstract method check_status should be implemented by derived class')        
        
    def default_resources(self, res) :
        raise RuntimeError('abstract method sub_script_head should be implemented by derived class')        

    def sub_script_head(self, res) :
        raise RuntimeError('abstract method sub_script_head should be implemented by derived class')        

    def sub_script_cmd(self, cmd, res):
        raise RuntimeError('abstract method sub_script_cmd should be implemented by derived class')        

    def do_submit(self,
                  job_dirs,
                  cmd,
                  args = None, 
                  res = None,
                  outlog = 'log',
                  errlog = 'err'):
        '''
        submit a single job, assuming that no job is running there.
        '''
        raise RuntimeError('abstract method check_status should be implemented by derived class')        

    def sub_script(self,
                   job_dirs,
                   cmd,
                   args = None,
                   res  = None,
                   outlog = 'log',
                   errlog = 'err') :
        """
        make submit script

        job_dirs(list):         directories of jobs. size: n_job
        cmd(list):              commands to be executed. size: n_cmd
        args(list of list):     args of commands. size of n_cmd x n_job
                                can be None
        res(dict):              resources available
        outlog(str):            file name for output
        errlog(str):            file name for error
        """
        res = self.default_resources(res)
        ret = self.sub_script_head(res)
        if not isinstance(cmd, list):
            cmd = [cmd]
        if args == None :
            args = []
            for ii in cmd:
                _args = []
                for jj in job_dirs:
                    _args.append('')
                args.append(_args)
        # loop over commands 
        self.cmd_cnt = 0
        try:
            self.manual_cuda_devices = res['manual_cuda_devices']
        except KeyError:
            self.manual_cuda_devices = 0
        try:
            self.manual_cuda_multiplicity = res['manual_cuda_multiplicity']
        except KeyError:
            self.manual_cuda_multiplicity = 1
        for ii in range(len(cmd)):            
            # for one command
            ret += self._sub_script_inner(job_dirs,
                                          cmd[ii],
                                          args[ii],
                                          ii,
                                          res,
                                          outlog=outlog,
                                          errlog=errlog)
        ret += '\ntouch %s\n' % self.finish_tag_name
        return ret

    def submit(self,
               job_dirs,
               cmd,
               args = None,
               res = None,
               restart = False,
               outlog = 'log',
               errlog = 'err'):
        if restart:
            dlog.debug('restart task')
            status = self.check_status()
            if status in [  JobStatus.unsubmitted, JobStatus.unknown, JobStatus.terminated ]:
                dlog.debug('task restart point !!!')
                self.do_submit(job_dirs, cmd, args, res, outlog=outlog, errlog=errlog)
            elif status==JobStatus.waiting:
                dlog.debug('task is waiting')
            elif status==JobStatus.running:
                dlog.debug('task is running')
            elif status==JobStatus.finished:
                dlog.debug('task is finished')
            else:
                raise RuntimeError('unknow job status, must be wrong')
        else:
            dlog.debug('new task')
            self.do_submit(job_dirs, cmd, args, res, outlog=outlog, errlog=errlog)
        if res is None:
            sleep = 0
        else:
            sleep = res.get('submit_wait_time', 0)
        time.sleep(sleep) # For preventing the crash of the tasks while submitting

    def check_finish_tag(self) :
        return self.context.check_file_exists(self.finish_tag_name)

    def _sub_script_inner(self, 
                          job_dirs,
                          cmd,
                          args,
                          idx,
                          res,
                          outlog = 'log',
                          errlog = 'err') :
        ret = ""
        allow_failure = res.get('allow_failure', False)
        for ii,jj in zip(job_dirs, args) :
            ret += 'cd %s\n' % ii
            ret += 'test $? -ne 0 && exit 1\n\n'
            if self.manual_cuda_devices > 0:
                # set CUDA_VISIBLE_DEVICES
                ret += 'export CUDA_VISIBLE_DEVICES=%d\n' % (self.cmd_cnt % self.manual_cuda_devices)
            ret += '{ if [ ! -f tag_%d_finished ] ;then\n' % idx
            ret += '  %s 1>> %s 2>> %s \n' % (self.sub_script_cmd(cmd, jj, res), outlog, errlog)
            if res['allow_failure'] is False:
                ret += '  if test $? -ne 0; then exit 1; else touch tag_%d_finished; fi \n' % idx
            else :
                ret += '  if test $? -ne 0; then touch tag_failure_%d; fi \n' % idx
                ret += '  touch tag_%d_finished \n' % idx
            ret += 'fi }'
            if self.manual_cuda_devices > 0:
                ret += '&'
                self.cmd_cnt += 1
            ret += '\n\n'
            ret += 'cd %s\n' % self.context.remote_root
            ret += 'test $? -ne 0 && exit 1\n'
            if self.manual_cuda_devices > 0 and self.cmd_cnt % (self.manual_cuda_devices * self.manual_cuda_multiplicity) == 0:
                ret += '\nwait\n\n'
        ret += '\nwait\n\n'
        return ret

import os,sys,time

from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog


class Batch(object) :
    def __init__ (self,
                  context) :
        self.context = context

    def check_status(self) :
        raise RuntimeError('abstract method check_status should be implemented by derived class')        
        
    def default_resources(self, res) :
        raise RuntimeError('abstract method sub_script_head should be implemented by derived class')        

    def sub_script_head(self, res) :
        raise RuntimeError('abstract method sub_script_head should be implemented by derived class')        

    def sub_script_cmd(self, cmd, res, errlog, outlog):
        raise RuntimeError('abstract method sub_script_cmd should be implemented by derived class')        

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
        for ii in range(len(cmd)):            
            # for one command
            ret += self._sub_script_inner(job_dirs,
                                          cmd[ii],
                                          args[ii],
                                          ii,
                                          res,
                                          outlog=outlog,
                                          errlog=errlog)
        ret += '\ntouch tag_finished\n'
        return ret


    def do_submit(self) :
        '''
        submit a single job, assuming that no job is running there.
        '''
        raise RuntimeError('abstract method check_status should be implemented by derived class')        

    def submit(self,
               job_dirs,
               cmd,
               args = None,
               res = None,
               restart = False,
               sleep = 0):
        if restart:
            dlog.debug('restart task')
            status = self.check_status()
            if status in [  JobStatus.unsubmitted, JobStatus.unknown, JobStatus.terminated ]:
                dlog.debug('task restart point !!!')
                self.do_submit(job_dirs, cmd, args, res)
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
            self.do_submit(job_dirs, cmd, args, res)
        time.sleep(sleep) # For preventing the crash of the tasks while submitting        


    def _sub_script_inner(self, 
                          job_dirs,
                          cmd,
                          args,
                          idx,
                          res,
                          outlog = 'log',
                          errlog = 'err') :
        ret = ""
        try:
            allow_failure = res['allow_failure']
        except:
            allow_failure = False
        for ii,jj in zip(job_dirs, args) :
            ret += 'cd %s\n' % ii
            ret += 'test $? -ne 0 && exit\n\n'
            ret += 'if [ ! -f tag_%d_finished ] ;then\n' % idx
            ret += '  %s 1> %s 2> %s \n' % (self.sub_script_cmd(cmd, jj, res), outlog, errlog)
            if res['allow_failure'] is False:
                ret += '  if test $? -ne 0; then exit; else touch tag_%d_finished; fi \n' % idx
            ret += 'fi\n\n'
            ret += 'cd %s\n' % self.context.remote_root
            ret += 'test $? -ne 0 && exit\n'
        return ret
        

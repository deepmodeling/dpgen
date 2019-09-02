import os,sys,time

from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog


class Batch(object) :
    def __init__ (self,
                  context) :
        self.context = context

    def sub_cmd(self) :
        raise RuntimeError('abstract method sub_cmd should be implemented by derived class')

    def sub_script(self, job_dirs, cmd, args = None, res = None) :
        raise RuntimeError('abstract method sub_script should be implemented by derived class')

    def check_sub_limit(self, task_max, **kwarg) :
        raise RuntimeError('abstract method check_sub_limit should be implemented by derived class')

    def check_status(self) :
        raise RuntimeError('abstract method check_status should be implemented by derived class')        

    def get_job_id(self) :
        if self.context.check_file_exists('job_id') :
            return self.context.read_file('job_id')
        else:
            return ""

    def check_finish_tag(self) :
        return self.context.check_file_exists('tag_finished') 

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
                self._do_submit(job_dirs, cmd, args, res)
            elif status==JobStatus.waiting:
                dlog.debug('task is waiting')
            elif status==JobStatus.running:
                dlog.debug('task is running')
            else:
                dlog.debug('task is finished')
        else:
            dlog.debug('new task')
            self._do_submit(job_dirs, cmd, args, res)
        time.sleep(sleep) # For preventing the crash of the tasks while submitting        


    def _do_submit(self, 
                   job_dirs,
                   cmd,
                   args = None, 
                   res = None) :
        if res == None:
            res = {}
        if 'task_max' in res and res['task_max'] > 0:
            while self.check_sub_limit(task_max=res['task_max']):
                time.sleep(60)
        script_str = self.sub_script(job_dirs, cmd, args=args, res=res)
        self.context.write_file('run.sub', script_str)
        stdin, stdout, stderr = self.context.block_checkcall('cd %s; %s %s' % (self.context.remote_root, self.sub_cmd(), 'run.sub'))
        subret = (stdout.readlines())
        job_id = subret[0].split()[-1]
        self.context.write_file('job_id', job_id)

       

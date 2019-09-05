import os,sys,time

from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog


class Batch(object) :
    def __init__ (self,
                  context) :
        self.context = context

    def check_status(self) :
        raise RuntimeError('abstract method check_status should be implemented by derived class')        
        
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
            else:
                dlog.debug('task is finished')
        else:
            dlog.debug('new task')
            self.do_submit(job_dirs, cmd, args, res)
        time.sleep(sleep) # For preventing the crash of the tasks while submitting        


       

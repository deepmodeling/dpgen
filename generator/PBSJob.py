#!/usr/bin/env python3

import os
import sys
from enum               import Enum
from subprocess         import Popen, PIPE

from BatchJob import JobStatus
from BatchJob import BatchJob
        
class PBSJob (BatchJob) :
    def submit_command (self):
        return "qsub"
    def check_status (self):
        job_id = self.get_job_id ()
        if len(job_id) == 0 :
            return JobStatus.unsubmitted
        ret = Popen (["qstat", job_id],  stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if (ret.returncode != 0) :
            if str("qstat: Unknown Job Id") in str(stderr, encoding='ascii') :
                if os.path.exists (self.job_finish_tag) :
                    return JobStatus.finished
                else :
                    return JobStatus.terminated
            else :
                raise RuntimeError ("status command " + "qstat" + " fails to execute\n" + 
                                    "erro info: " + str(stderr, encoding='ascii') + "\n" +
                                    "return code: " + str(ret.returncode))
        status_line = str(stdout, encoding='ascii').split ('\n')[-2]
        status_word = status_line.split ()[-2]        
#        print (status_word)
        if      status_word in ["Q","H"] :
            return JobStatus.waiting
        elif    status_word in ["R"] :
            return JobStatus.running
        elif    status_word in ["C","E","K"] :
            if os.path.exists (self.job_finish_tag) :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown        

if __name__ == "__main__" :
    cwd = os.getcwd()
    job = PBSJob (cwd + "/test", "cu01.sleep")
    job.submit ()
    print ("submit done")
    stat = job.check_status ()
    print ("check done")
    print (stat)
    

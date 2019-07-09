#!/usr/bin/env python3

import os
import sys
from enum               import Enum
from subprocess         import Popen, PIPE
from dpgen.auto_test.lib.BatchJob import BatchJob
from dpgen.auto_test.lib.BatchJob import JobStatus
        
class SlurmJob (BatchJob) :
    def submit_command (self):
        return "sbatch"
    def check_status (self):
        job_id = self.get_job_id ()
        if len(job_id) == 0 :
            return JobStatus.unsubmitted
        ret = Popen (["squeue --job " + job_id], shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if (ret.returncode != 0) :
            if str("Invalid job id specified") in str(stderr, encoding='ascii') :
                if os.path.exists (self.job_finish_tag) :
                    return JobStatus.finished
                else :
                    return JobStatus.terminated
            else :
                Logger.error ("status command " + "squeue" + " fails to execute")
                Logger.error ("erro info: " + str(stderr, encoding='ascii'))
                Logger.error ("return code: " + str(ret.returncode))
                sys.exit ()
        status_line = str(stdout, encoding='ascii').split ('\n')[-2]
        status_word = status_line.split ()[4]
#        status_word = status_line.split ()[-4]
#        print ("status line: " + status_line)
#        print ("status word: " + status_word)
#        print (status_word)
        if      status_word in ["PD","CF","S"] :
            return JobStatus.waiting
        elif    status_word in ["R","CG"] :
            return JobStatus.running
        elif    status_word in ["C","E","K","BF","CA","CD","F","NF","PR","SE","ST","TO"] :
            if os.path.exists (self.job_finish_tag) :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown        

if __name__ == "__main__" :
    job = SlurmJob ("/home/han.wang/data/test/string/test", "cu01.sleep")
    job.submit ()
    print ("submit done")
    stat = job.check_status ()
    print ("check done")
    print (stat)
    

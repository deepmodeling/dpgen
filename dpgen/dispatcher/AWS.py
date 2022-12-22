import os,getpass,time
from datetime import datetime
from itertools import zip_longest
from dpgen.dispatcher.Batch import Batch 
from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog


class AWS(Batch):
    _query_time_interval = 30
    _job_id_map_status = {}
    _jobQueue = ""
    _query_next_allow_time = datetime.now().timestamp()

    def __init__(self, context, uuid_names=True):
        import boto3
        self.batch_client = boto3.client('batch')
        super().__init__(context, uuid_names)

    @staticmethod
    def map_aws_status_to_dpgen_status(aws_status):
        map_dict = {'SUBMITTED': JobStatus.waiting,
                'PENDING': JobStatus.waiting,
                'RUNNABLE': JobStatus.waiting,
                'STARTING': JobStatus.waiting,
                'RUNNING': JobStatus.running,
                'SUCCEEDED': JobStatus.finished,
                'FAILED': JobStatus.terminated,
                'UNKNOWN': JobStatus.unknown}
        return map_dict.get(aws_status, JobStatus.unknown)

    @classmethod
    def AWS_check_status(cls, job_id=""):
        """
        to aviod query jobStatus too often, set a time interval
        query_dict example:
        {job_id: JobStatus}

        {'40fb24b2-d0ca-4443-8e3a-c0906ea03622': <JobStatus.running: 3>,
         '41bda50c-0a23-4372-806c-87d16a680d85': <JobStatus.waiting: 2>}
           
        """
        query_dict ={}
        if datetime.now().timestamp() > cls._query_next_allow_time:
            cls._query_next_allow_time=datetime.now().timestamp()+cls._query_time_interval
            for status in ['SUBMITTED', 'PENDING', 'RUNNABLE', 'STARTING', 'RUNNING','SUCCEEDED', 'FAILED']:
                nextToken = ''
                while nextToken is not None:
                    status_response = self.batch_client.list_jobs(jobQueue=cls._jobQueue, jobStatus=status, maxResults=100, nextToken=nextToken)
                    status_list=status_response.get('jobSummaryList')
                    nextToken = status_response.get('nextToken', None)
                    for job_dict in status_list:
                        cls._job_id_map_status.update({job_dict['jobId']: cls.map_aws_status_to_dpgen_status(job_dict['status'])})
            dlog.debug('20000:_map: %s' %(cls._job_id_map_status))
        dlog.debug('62000:job_id:%s, _query: %s, _map: %s' %(job_id, query_dict, cls._job_id_map_status))
        if job_id:
            return cls._job_id_map_status.get(job_id)
                    
        return cls._job_id_map_status

    @property
    def job_id(self):
        try:
            self._job_id
        except AttributeError:
            if self.context.check_file_exists(self.job_id_name):
                self._job_id = self.context.read_file(self.job_id_name)
                response_list = self.batch_client.describe_jobs(jobs=[self._job_id]).get('jobs')
                try:
                    response = response_list[0]
                    jobQueue = response['jobQueue']
                except IndexError:
                    pass
                else:
                    self.job_id = (response, jobQueue)
                    return self._job_id
            dlog.debug("50000, self._job_id:%s,_Queue:%s,_map:%s,"%(self._job_id, self.__class__._jobQueue, self.__class__._job_id_map_status  ))       
            return ""
        return self._job_id

    @job_id.setter
    def job_id(self, values):
        response, jobQueue = values
        self._job_id = response['jobId']
        self._job_name = response['jobName']
        self.__class__._jobQueue = jobQueue
        self.__class__._job_id_map_status[self._job_id] = self.map_aws_status_to_dpgen_status(response.get('status', 'SUBMITTED'))
        self.context.write_file(self.job_id_name, self._job_id)
        dlog.debug("15000, _job_id:%s, _job_name:%s, _map:%s, _Queue:%s" % (self._job_id, self._job_name, self.__class__._job_id_map_status, self.__class__._jobQueue))

    def check_status(self):
        return self.__class__.AWS_check_status(job_id=self.job_id)
    
    def sub_script(self, job_dirs, cmd, args, res, outlog, errlog):
        if args is None:
            args=[]
        multi_command = ""
        for job_dir in job_dirs:
            for idx,t in enumerate(zip_longest(cmd, args, fillvalue='')):
                c_str =  f"cd {self.context.remote_root}/{job_dir} && ( test -f tag_{idx}_finished || ( ({t[0]} {t[1]} && touch tag_{idx}_finished  2>>{errlog} || exit 52 ) | tee -a {outlog}) ) || exit 51;"
                multi_command += c_str
        multi_command +="exit 0;"
        dlog.debug("10000, %s" % multi_command)
        return multi_command
        
    def default_resources(self, res):
        if res == None:
            res = {}
        else:
            # res.setdefault(jobDefinition)
            res.setdefault('cpu_num', 32)
            res.setdefault('memory_size', 120000)
            res.setdefault('jobQueue', 'deepmd_m5_v1_7')
        return res
    
    def do_submit(self,
            job_dirs,
            cmd,
            args = None,
            res = None,
            outlog = 'log',
            errlog = 'err'):

        res = self.default_resources(res)
        dlog.debug("2000, params=(%s, %s, %s, %s, %s, %s, )" % (job_dirs, cmd, args, res,  outlog, errlog ))
        dlog.debug('2200, self.context.remote_root: %s , self.context.local_root: %s' % (self.context.remote_root, self.context.local_root))
        # concreate_command = 
        script_str = self.sub_script(job_dirs, cmd, args=args, res=res, outlog=outlog, errlog=errlog)
        dlog.debug('2300, script_str: %s, self.sub_script_name: %s' % (script_str, self.sub_script_name))
        """
        jobName example:
        home-ec2-user-Ag_init-run_gen-iter_000000-01_model_devi-task_000_000048
        """
        jobName = os.path.join(self.context.remote_root,job_dirs.pop())[1:].replace('/','-').replace('.','_')
        jobName += ("_" + str(self.context.job_uuid))
        response = self.batch_client.submit_job(jobName=jobName, 
                jobQueue=res['jobQueue'], 
                jobDefinition=res['jobDefinition'],
                parameters={'task_command':script_str},
                containerOverrides={'vcpus':res['cpu_num'], 'memory':res['memory_size']})
        dlog.debug('4000, response:%s' % response)
        self.job_id = (response, res['jobQueue'])

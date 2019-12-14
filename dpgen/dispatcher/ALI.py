from aliyunsdkecs.request.v20140526.DescribeInstancesRequest import DescribeInstancesRequest
from aliyunsdkcore.client import AcsClient
from aliyunsdkcore.acs_exception.exceptions import ClientException
from aliyunsdkcore.acs_exception.exceptions import ServerException
from aliyunsdkecs.request.v20140526.RunInstancesRequest import RunInstancesRequest
from aliyunsdkecs.request.v20140526.DeleteInstancesRequest import DeleteInstancesRequest
import time, json, os, glob
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks
from os import listdir
from os.path import isfile, join

determine_machine = {
    "gpu": {
            1: "ecs.gn5-c8g1.2xlarge",
    },
    "cpu": {
            1: "ecs.c6.large",
            4: "ecs.c6.2xlarge",
            8: "ecs.c6.4xlarge"
    }
}

class ALI():
    def __init__(self, adata, mdata, nchunks, work_path, cwd):
        self.ip_list = None
        self.instance_list = None
        self.dispatchers = None
        self.job_handlers = None
        self.adata = adata
        self.mdata = mdata
        self.nchunks = nchunks
        self.work_path = work_path
        self.cwd = cwd
        self.regionID = 'cn-hangzhou'
        
    def init(self):
        if self.check_restart():
            pass
        else:
            self.create_machine()
            self.dispatchers = self.make_dispatchers()

    def check_restart(self):
        os.chdir(self.work_path)
        dispatchers = []
        instance_list = []
        if len(glob.glob('jr.*.json')) == self.nchunks:
            for ii in range(self.nchunks):
                with open('jr.%.06d.json' % ii) as fp:
                    job_record = json.load(fp)
                    key = list(job_record.keys())[0]
                    ip, instance_id = job_record[key]['context'][-2], job_record[key]['context'][-1]
                    instance_list.append(instance_id)
                    profile = {
                        'type': 'ALI',
                        'hostname': ip,
                        'instance_id': instance_id,
                        'port': 22,
                        'username': 'root',
                        'password': self.adata['password'],
                        'work_path': '/root/dpgen_work'
                    }
                    disp = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
                    max_check = 10
                    cnt = 0
                    while not disp.session._check_alive():
                        cnt += 1
                        if cnt == max_check:
                            break
                    if cnt != max_check:
                        dispatchers.append(disp)
        restart = False
        if len(dispatchers) == self.nchunks:
            restart = True
            self.dispatchers = dispatchers
            self.instance_list = instance_list
        os.chdir(self.cwd)
        return restart

    def run_jobs(self,
                 resources,
                 command,
                 work_path,
                 tasks,
                 group_size,
                 forward_common_files,
                 forward_task_files,
                 backward_task_files,
                 forward_task_deference = True,
                 outlog = 'log',
                 errlog = 'err'):
        task_chunks = _split_tasks(tasks, group_size)
        job_handlers = []
        for ii in range(self.nchunks):
            job_handler = self.dispatchers[ii].submit_jobs(resources,
                                                           command,
                                                           work_path,
                                                           task_chunks[ii],
                                                           group_size,
                                                           forward_common_files,
                                                           forward_task_files,
                                                           backward_task_files,
                                                           forward_task_deference,
                                                           outlog,
                                                           errlog)
            job_handlers.append(job_handler)
        while True:
            cnt = 0
            for ii in range(self.nchunks):
                if self.dispatchers[ii].all_finished(job_handlers[ii]):
                    cnt += 1
            if cnt == self.nchunks:
                break
            else:
                time.sleep(10)
        self.delete_machine()

    def make_dispatchers(self):
        dispatchers = []
        for ii in range(self.nchunks):
            remote_profile = {
                'type': 'ALI',
                'hostname': self.ip_list[ii],
                'instance_id': self.instance_list[ii],
                'port': 22,
                'username': 'root',
                'password': self.adata['password'],
                'work_path': '/root/dpgen_work'
            }
            disp = Dispatcher(remote_profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
            dispatchers.append(disp)
        return dispatchers

    def create_machine(self):
        AccessKey_ID = self.adata["AccessKey_ID"]
        AccessKey_Secret = self.adata["AccessKey_Secret"]
        strategy = self.adata["pay_strategy"]
        pwd = self.adata["password"]
        regionID = self.regionID
        instance_type = determine_machine[self.mdata['partition']][self.mdata['numb_gpu']]
        if True:
            client = AcsClient(AccessKey_ID,AccessKey_Secret, 'cn-hangzhou')
            request = RunInstancesRequest()
            request.set_accept_format('json')
            request.set_UniqueSuffix(True)
            request.set_Password(pwd)
            request.set_Amount(self.nchunks)
            request.set_LaunchTemplateName(instance_type + '_cn-hangzhou_i')
            response = client.do_action_with_exception(request)
            response = json.loads(response)
            self.instance_list = response["InstanceIdSets"]["InstanceIdSet"]
            time.sleep(50)
            request = DescribeInstancesRequest()
            request.set_accept_format('json')
            request.set_InstanceIds(self.instance_list)
            response = client.do_action_with_exception(request)
            response = json.loads(response)
            ip = []
            for i in range(len(response["Instances"]["Instance"])):
                ip.append(response["Instances"]["Instance"][i]["PublicIpAddress"]['IpAddress'][0])
            self.ip_list = ip
        else:
            return "create failed"

    def delete_machine(self):
        AccessKey_ID = self.adata["AccessKey_ID"]
        AccessKey_Secret = self.adata["AccessKey_Secret"]
        regionID = self.regionID
        client = AcsClient(AccessKey_ID,AccessKey_Secret, regionID)
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        request.set_InstanceIds(self.instance_list)
        request.set_Force(True)
        response = client.do_action_with_exception(request)


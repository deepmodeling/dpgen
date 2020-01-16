from aliyunsdkecs.request.v20140526.DescribeInstancesRequest import DescribeInstancesRequest
from aliyunsdkcore.client import AcsClient
from aliyunsdkcore.acs_exception.exceptions import ClientException
from aliyunsdkcore.acs_exception.exceptions import ServerException
from aliyunsdkecs.request.v20140526.RunInstancesRequest import RunInstancesRequest
from aliyunsdkecs.request.v20140526.DeleteInstancesRequest import DeleteInstancesRequest
import time, json, os, glob
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks
from os.path import join

class ALI():
    def __init__(self, adata, mdata_resources, mdata_machine, nchunks, work_path):
        self.ip_list = None
        self.instance_list = None
        self.dispatchers = None
        self.adata = adata
        self.mdata_resources = mdata_resources
        self.mdata_machine = mdata_machine
        self.nchunks = nchunks
        self.work_path = work_path
        
    def init(self):
        if self.check_restart():
            pass
        else:
            self.create_machine()
            self.dispatchers = self.make_dispatchers()

    def check_restart(self):
        dispatchers = []
        instance_list = []
        if os.path.exists('machine_record.json'):
            with open('machine_record.json', 'r') as fp:
                machine_record = json.load(fp)
                for ii in range(self.nchunks):
                    ip, instance_id = machine_record['ip'][ii], machine_record['instance_id'][ii]
                    instance_list.append(instance_id)
                    profile = self.mdata_machine.copy()
                    profile['hostname'] = ip
                    profile['instance_id'] = instance_id
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
            return restart
        else:
            return False

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
            profile = self.mdata_machine.copy()
            profile['hostname'] = self.ip_list[ii]
            profile['instance_id'] = self.instance_list[ii]
            disp = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
            dispatchers.append(disp)
        return dispatchers

    def create_machine(self):
        AccessKey_ID = self.adata["AccessKey_ID"]
        AccessKey_Secret = self.adata["AccessKey_Secret"]
        strategy = self.adata["pay_strategy"]
        pwd = self.adata["password"]
        regionID = self.mdata_machine['regionID']
        template_name = '%s_%s_%s' % (self.mdata_resources['partition'], self.mdata_resources['numb_gpu'], strategy)
        instance_name = self.adata["instance_name"]
        client = AcsClient(AccessKey_ID, AccessKey_Secret, regionID)
        request = RunInstancesRequest()
        request.set_accept_format('json')
        request.set_UniqueSuffix(True)
        request.set_Password(pwd)
        request.set_InstanceName(instance_name)
        request.set_Amount(self.nchunks)
        request.set_LaunchTemplateName(template_name)
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
        with open('machine_record.json', 'w') as fp:
            json.dump({'ip': self.ip_list, 'instance_id': self.instance_list}, fp, indent=4)
<<<<<<< HEAD

=======
        
>>>>>>> 919c46a94f950c32eff6e7c308bafb3f3d14233f
    def delete_machine(self):
        AccessKey_ID = self.adata["AccessKey_ID"]
        AccessKey_Secret = self.adata["AccessKey_Secret"]
        regionID = self.mdata_machine['regionID']
        client = AcsClient(AccessKey_ID,AccessKey_Secret, regionID)
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        request.set_InstanceIds(self.instance_list)
        request.set_Force(True)
        response = client.do_action_with_exception(request)
        os.remove('machine_record.json')

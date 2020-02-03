from aliyunsdkecs.request.v20140526.DescribeInstancesRequest import DescribeInstancesRequest
from aliyunsdkcore.client import AcsClient
from aliyunsdkcore.acs_exception.exceptions import ClientException
from aliyunsdkcore.acs_exception.exceptions import ServerException
from aliyunsdkecs.request.v20140526.RunInstancesRequest import RunInstancesRequest
from aliyunsdkecs.request.v20140526.DeleteInstancesRequest import DeleteInstancesRequest
import time, json, os, glob
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks
from os.path import join
from dpgen  import dlog

def manual_delete(regionID):
    with open('machine-ali.json') as fp1:
        mdata = json.load(fp1)
        AccessKey_ID = mdata['train'][0]['machine']['ali_auth']['AccessKey_ID']
        AccessKey_Secret = mdata['train'][0]['machine']['ali_auth']['AccessKey_Secret']
        with open('machine_record.json', 'r') as fp2:
            machine_record = json.load(fp2)
            instance_list = machine_record['instance_id']
            client = AcsClient(AccessKey_ID, AccessKey_Secret, regionID)
            request = DeleteInstancesRequest()
            request.set_accept_format('json')
            if len(instance_list) <= 100:
                request.set_InstanceIds(instance_list)
                request.set_Force(True)
                response = client.do_action_with_exception(request)
            else:
                iteration = len(instance_list) // 100
                for i in range(iteration):
                    request.set_InstanceIds(instance_list[i*100:(i+1)*100])
                    request.set_Force(True)
                    response = client.do_action_with_exception(request)
                if len(instance_list) - iteration * 100 != 0:
                    request.set_InstanceIds(instance_list[iteration*100:])
                    request.set_Force(True)
                    response = client.do_action_with_exception(request)
            os.remove('machine_record.json')

class ALI():
    def __init__(self, adata, mdata_resources, mdata_machine, nchunks):
        self.ip_list = None
        self.instance_list = None
        self.regionID = mdata_machine["regionID"]
        self.dispatchers = None
        self.adata = adata
        self.mdata_resources = mdata_resources
        self.mdata_machine = mdata_machine
        self.nchunks = nchunks

        
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
        tmp = [1 for ii in range(self.nchunks)]
        while True:
            for ii in range(self.nchunks):
                if tmp[ii] and self.dispatchers[ii].all_finished(job_handlers[ii]):
                    self.delete(self.instance_list[ii])
                    tmp[ii] = 0
            if sum(tmp) == 0:
                os.remove('machine_record.json')
                break
            else:
                time.sleep(10)
        #self.delete_machine()

    def delete(self, instance_list):
        AccessKey_ID = self.adata["AccessKey_ID"]
        AccessKey_Secret = self.adata["AccessKey_Secret"]
        regionID = self.regionID
        client = AcsClient(AccessKey_ID, AccessKey_Secret, regionID)
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        request.set_InstanceIds([instance_list])
        request.set_Force(True)
        response = client.do_action_with_exception(request)

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
        regionID = self.regionID
        template_name = '%s_%s_%s' % (self.mdata_resources['partition'], self.mdata_resources['numb_gpu'], strategy)
        instance_name = self.adata["instance_name"]
        client = AcsClient(AccessKey_ID, AccessKey_Secret, regionID)
        self.instance_list = []
        self.ip_list = []
        request = RunInstancesRequest()
        request.set_accept_format('json')
        request.set_UniqueSuffix(True)
        request.set_Password(pwd)
        request.set_InstanceName(instance_name)
        request.set_LaunchTemplateName(template_name)
        if self.nchunks <= 100:
            request.set_Amount(self.nchunks)
            try:
                response = client.do_action_with_exception(request)
                response = json.loads(response)
                for instanceID in response["InstanceIdSets"]["InstanceIdSet"]:
                    self.instance_list.append(instanceID)
            except:
                dlog.debug("Create failed, please check the console.")
                if len(self.instance_list) > 0:
                    self.delete_machine()
        else:
            iteration = self.nchunks // 100 
            try:
                for i in range(iteration):
                    request.set_Amount(100)
                    response = client.do_action_with_exception(request)
                    response = json.loads(response)
                    for instanceID in response["InstanceIdSets"]["InstanceIdSet"]:
                        self.instance_list.append(instanceID)
            except:
                dlog.debug("Create failed, please check the console.")
                if len(self.instance_list) > 0:
                    self.delete_machine()
            if self.nchunks - iteration * 100 != 0:
                try:
                    request.set_Amount(self.nchunks - iteration * 100)
                    response = client.do_action_with_exception(request)
                    response = json.loads(response)
                    for instanceID in response["InstanceIdSets"]["InstanceIdSet"]:
                        self.instance_list.append(instanceID)
                except:
                    dlog.debug("Create failed, please check the console.")
                    if len(self.instance_list) > 0:
                        self.delete_machine()
        time.sleep(60)
        request = DescribeInstancesRequest()
        request.set_accept_format('json')
        if len(self.instance_list) <= 10:
            request.set_InstanceIds(self.instance_list)
            response = client.do_action_with_exception(request)
            response = json.loads(response)
            for i in range(len(response["Instances"]["Instance"])):
                self.ip_list.append(response["Instances"]["Instance"][i]["PublicIpAddress"]['IpAddress'][0])
        else:
            iteration = len(self.instance_list) // 10
            for i in range(iteration):
                request.set_InstanceIds(self.instance_list[i*10:(i+1)*10])
                response = client.do_action_with_exception(request)
                response = json.loads(response)
                for j in range(len(response["Instances"]["Instance"])):
                    self.ip_list.append(response["Instances"]["Instance"][j]["PublicIpAddress"]['IpAddress'][0])
            if len(self.instance_list) - iteration * 10 != 0:
                request.set_InstanceIds(self.instance_list[iteration*10:])
                response = client.do_action_with_exception(request)
                response = json.loads(response)
                for j in range(len(response["Instances"]["Instance"])):
                    self.ip_list.append(response["Instances"]["Instance"][j]["PublicIpAddress"]['IpAddress'][0])
        with open('machine_record.json', 'w') as fp:
            json.dump({'ip': self.ip_list, 'instance_id': self.instance_list}, fp, indent=4)

    def delete_machine(self):
        AccessKey_ID = self.adata["AccessKey_ID"]
        AccessKey_Secret = self.adata["AccessKey_Secret"]
        regionID = self.mdata_machine['regionID']
        client = AcsClient(AccessKey_ID,AccessKey_Secret, regionID)
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        if len(self.instance_list) <= 100:
            request.set_InstanceIds(self.instance_list)
            request.set_Force(True)
            response = client.do_action_with_exception(request)
        else:
            iteration = len(self.instance_list) // 100
            for i in range(iteration):
                request.set_InstanceIds(self.instance_list[i*100:(i+1)*100])
                request.set_Force(True)
                response = client.do_action_with_exception(request)
            if len(self.instance_list) - iteration * 100 != 0:
                request.set_InstanceIds(self.instance_list[iteration*100:])
                request.set_Force(True)
                response = client.do_action_with_exception(request)
        self.instance_list = []
        os.remove('machine_record.json')
        dlog.debug("Successfully free the machine!")

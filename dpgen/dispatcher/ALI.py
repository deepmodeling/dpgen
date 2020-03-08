from aliyunsdkecs.request.v20140526.DescribeInstancesRequest import DescribeInstancesRequest
from aliyunsdkcore.client import AcsClient
from aliyunsdkcore.acs_exception.exceptions import ClientException
from aliyunsdkcore.acs_exception.exceptions import ServerException
from aliyunsdkecs.request.v20140526.RunInstancesRequest import RunInstancesRequest
from aliyunsdkecs.request.v20140526.DeleteInstancesRequest import DeleteInstancesRequest
import time, json, os, glob
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, JobRecord
from os.path import join
from dpgen  import dlog
from hashlib import sha1

def manual_delete(stage):
    with open('machine-ali.json') as fp1:
        mdata = json.load(fp1)
        adata = mdata[stage][0]['machine']['ali_auth']
        mdata_resources = mdata[stage][0]['resources']
        mdata_machine = mdata[stage][0]['machine']
        ali = ALI(adata, mdata_resources, mdata_machine, 0)
        with open('machine_record.json', 'r') as fp2:
            machine_record = json.load(fp2)
            ali.instance_list = machine_record['instance_id']
            ali.delete_machine()

def manual_create(stage, machine_number):
    with open('machine-ali.json') as fp:
        mdata = json.load(fp)
        adata = mdata[stage][0]['machine']['ali_auth']
        mdata_resources = mdata[stage][0]['resources']
        mdata_machine = mdata[stage][0]['machine']
        ali = ALI(adata, mdata_resources, mdata_machine, machine_number)
        ali.alloc_machine()
        print(ali.ip_list)

class ALI():
    def __init__(self, adata, mdata_resources, mdata_machine, nchunks):
        self.ip_list = []
        self.instance_list = []
        self.dispatchers = None
        self.job_handlers = None
        self.task_chunks = None
        self.adata = adata
        self.regionID = adata["regionID"]
        self.client = AcsClient(adata["AccessKey_ID"], adata["AccessKey_Secret"], self.regionID)
        self.mdata_resources = mdata_resources
        self.mdata_machine = mdata_machine
        self.nchunks = nchunks

    def init(self, work_path, tasks, group_size):
        if self.check_restart(work_path, tasks, group_size):
            pass
        else:
            self.alloc_machine()
            self.dispatchers = self.make_dispatchers()

    def check_restart(self, work_path, tasks, group_size):
        if os.path.exists('machine_record.json'):
            dispatchers = []
            instance_list = []
            ip_list = []
            self.task_chunks = []
            task_chunks = _split_tasks(tasks, group_size)
            task_chunks_str = ['+'.join(ii) for ii in task_chunks]
            task_hashes = [sha1(ii.encode('utf-8')).hexdigest() for ii in task_chunks_str]
            nchunks = len(task_chunks)
            for ii in range(nchunks):
                job_record = JobRecord(work_path, task_chunks, fname = 'jr.%.06d.json' % ii)
                cur_chunk = task_chunks[ii]
                cur_hash = task_hashes[ii]
                if not job_record.check_finished(cur_hash): 
                    self.task_chunks.append(cur_chunk)
                    with open(os.path.join(work_path, 'jr.%.06d.json' % ii)) as fp:
                        jr = json.load(fp)
                        ip = jr[cur_hash]['context'][3]
                        instance_id = jr[cur_hash]['context'][4]
                        ip_list.append(ip)
                        instance_list.append(instance_id)
                        profile = self.mdata_machine.copy()
                        profile['hostname'] = ip
                        profile['instance_id'] = instance_id
                        disp = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
                        dispatchers.append(disp)
            self.dispatchers = dispatchers
            self.instance_list = instance_list
            self.ip_list = ip_list
            return True
        else:
            return False

    def alloc_machine(self):
        for district, type_num in self.adata["avail_resources"].items():
            for machine_type, number in type_num.items():
                if self.nchunks > number:
                    self.nchunks -= number
                    self.create_machine(machine_type, number, district)
                elif self.nchunks > 0:
                    self.create_machine(machine_type, self.nchunks, district)
        time.sleep(60)
        self.get_ip()

    def create_machine(self, machine_type, number, district):
        request = RunInstancesRequest()
        request.set_accept_format('json')
        request.set_UniqueSuffix(True)
        request.set_Password(self.mdata_machine["password"])
        request.set_InstanceName(self.adata['instance_name'])

        if self.mdata_resources['partition'] == 'gpu':
            template_name = 'gpu_%s_%s_%s_%s' % (self.mdata_resources['numb_gpu'], self.adata["pay_strategy"], district, machine_type)
        elif self.mdata_resources['partition'] == 'cpu':
            template_name = 'cpu_%s_%s_%s_%s' % (self.mdata_resources['task_per_node'], self.adata["pay_strategy"], district, machine_type)
        request.set_LaunchTemplateName(template_name)

        if number <= 100 and number > 0:
            request.set_Amount(number)
            response = self.client.do_action_with_exception(request)
            response = json.loads(response)
            for instanceID in response["InstanceIdSets"]["InstanceIdSet"]:
                self.instance_list.append(instanceID)
        else:
            iteration = number // 100
            for i in range(iteration):
                request.set_Amount(100)
                response = self.client.do_action_with_exception(request)
                response = json.loads(response)
                for instanceID in response["InstanceIdSets"]["InstanceIdSet"]:
                    self.instance_list.append(instanceID)
            if number - iteration * 100 != 0:
                request.set_Amount(number - iteration * 100)
                response = self.client.do_action_with_exception(request)
                response = json.loads(response)
                for instanceID in response["InstanceIdSets"]["InstanceIdSet"]:
                    self.instance_list.append(instanceID)

    def get_ip(self):
        request = DescribeInstancesRequest()
        request.set_accept_format('json')
        if len(self.instance_list) <= 10:
            for i in range(len(self.instance_list)):
                request.set_InstanceIds([self.instance_list[i]])
                response = self.client.do_action_with_exception(request)
                response = json.loads(response)
                self.ip_list.append(response["Instances"]["Instance"][0]["PublicIpAddress"]['IpAddress'][0])
        else:
            iteration = len(self.instance_list) // 10
            for i in range(iteration):
                for j in range(10):
                    request.set_InstanceIds([self.instance_list[i*10+j]])
                    response = self.client.do_action_with_exception(request)
                    response = json.loads(response)
                    self.ip_list.append(response["Instances"]["Instance"][0]["PublicIpAddress"]['IpAddress'][0])
            if len(self.instance_list) - iteration * 10 != 0:
                for j in range(len(self.instance_list) - iteration * 10):
                    request.set_InstanceIds([self.instance_list[iteration*10+j]])
                    response = self.client.do_action_with_exception(request)
                    response = json.loads(response)
                    self.ip_list.append(response["Instances"]["Instance"][0]["PublicIpAddress"]['IpAddress'][0])
        dlog.info('create machine successfully, following are the ip addresses')
        for ip in self.ip_list:
            dlog.info(ip)
        with open('machine_record.json', 'w') as fp:
            json.dump({'ip': self.ip_list, 'instance_id': self.instance_list}, fp, indent=4)

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
        if not self.task_chunks:
            self.task_chunks = _split_tasks(tasks, group_size)
        self.job_handlers = []
        for ii in range(len(self.dispatchers)):
            job_handler = self.dispatchers[ii].submit_jobs(resources,
                                                           command,
                                                           work_path,
                                                           self.task_chunks[ii],
                                                           group_size,
                                                           forward_common_files,
                                                           forward_task_files,
                                                           backward_task_files,
                                                           forward_task_deference,
                                                           outlog,
                                                           errlog)
            self.job_handlers.append(job_handler)
        while True:
            for ii in range(len(self.dispatchers)):
                if self.dispatchers[ii].all_finished(self.job_handlers[ii]):
                    self.delete(ii)
                    break
            if len(self.dispatchers) == 0:
                break
            else:
                time.sleep(10)

    def delete(self, ii):
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        request.set_InstanceIds([self.instance_list[ii]])
        request.set_Force(True)
        response = self.client.do_action_with_exception(request)
        self.nchunks -= 1
        self.instance_list.pop(ii)
        self.ip_list.pop(ii)
        self.dispatchers.pop(ii)
        self.job_handlers.pop(ii)
        if len(self.ip_list) == 0:
            os.remove('machine_record.json')
        else:
            with open('machine_record.json', 'w') as fp:
                json.dump({'ip': self.ip_list, 'instance_id': self.instance_list}, fp, indent=4)

    def make_dispatchers(self):
        dispatchers = []
        for ii in range(self.nchunks):
            profile = self.mdata_machine.copy()
            profile['hostname'] = self.ip_list[ii]
            profile['instance_id'] = self.instance_list[ii]
            disp = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
            dispatchers.append(disp)
        return dispatchers

    def delete_machine(self):
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        if len(self.instance_list) <= 100:
            request.set_InstanceIds(self.instance_list)
            request.set_Force(True)
            response = self.client.do_action_with_exception(request)
        else:
            iteration = len(self.instance_list) // 100
            for i in range(iteration):
                request.set_InstanceIds(self.instance_list[i*100:(i+1)*100])
                request.set_Force(True)
                response = self.client.do_action_with_exception(request)
            if len(self.instance_list) - iteration * 100 != 0:
                request.set_InstanceIds(self.instance_list[iteration*100:])
                request.set_Force(True)
                response = self.client.do_action_with_exception(request)
        self.instance_list = []
        self.ip_list = []
        os.remove('machine_record.json')
        dlog.debug("Successfully free the machine!")

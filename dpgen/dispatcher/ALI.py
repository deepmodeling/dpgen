from aliyunsdkecs.request.v20140526.DescribeInstancesRequest import DescribeInstancesRequest
from aliyunsdkcore.client import AcsClient
from aliyunsdkcore.acs_exception.exceptions import ClientException
from aliyunsdkcore.acs_exception.exceptions import ServerException
from aliyunsdkecs.request.v20140526.RunInstancesRequest import RunInstancesRequest
from aliyunsdkecs.request.v20140526.DeleteInstancesRequest import DeleteInstancesRequest
from aliyunsdkecs.request.v20140526.DescribeAutoProvisioningGroupInstancesRequest import DescribeAutoProvisioningGroupInstancesRequest
from aliyunsdkecs.request.v20140526.CreateAutoProvisioningGroupRequest import CreateAutoProvisioningGroupRequest
from aliyunsdkecs.request.v20140526.DeleteAutoProvisioningGroupRequest import DeleteAutoProvisioningGroupRequest
from aliyunsdkecs.request.v20140526.ModifyAutoProvisioningGroupRequest import ModifyAutoProvisioningGroupRequest
from aliyunsdkecs.request.v20140526.DeleteLaunchTemplateRequest import DeleteLaunchTemplateRequest
from aliyunsdkvpc.request.v20160428.DescribeVpcsRequest import DescribeVpcsRequest
from aliyunsdkecs.request.v20140526.DescribeLaunchTemplatesRequest import DescribeLaunchTemplatesRequest
from aliyunsdkecs.request.v20140526.CreateLaunchTemplateRequest import CreateLaunchTemplateRequest
from aliyunsdkecs.request.v20140526.DescribeImagesRequest import DescribeImagesRequest
from aliyunsdkecs.request.v20140526.DescribeSecurityGroupsRequest import DescribeSecurityGroupsRequest
import time, json, os, glob, string, random
from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, JobRecord
from dpgen.dispatcher.SSHContext import SSHSession
from dpgen.dispatcher.DispatcherList import DispatcherList, Entity
from os.path import join
from dpgen  import dlog
from hashlib import sha1

# cloud_resources = {"AccessKey_ID":"",
#                    "AccessKey_Secret":"",
#                    "regionID": "cn-shenzhen",
#                    "img_name": "kit",
#                    "machine_type_price": [
#                        {"machine_type": "ecs.gn6v-c8g1.2xlarge", "price_limit": 20.00, "numb": 1, "priority": 0},
#                        {"machine_type": "ecs.gn5-c4g1.xlarge",   "price_limit": 20.00, "numb": 1, "priority": 1}
#                    ],
#                    "instance_name": "CH4_test_username", 
#                    "pay_strategy": "spot"
#                    "apg_id": apg_id, 
#                    "template_id": template_id, 
#                    "vsw_id": vsw_id, 
#                    "region_id": region_id,
#                    "client": client}


class ALI(DispatcherList):
    def __init__(self, mdata_machine, mdata_resources, work_path, run_tasks, group_size, cloud_resources=None):
        super().__init__(mdata_machine, mdata_resources, work_path, run_tasks, group_size, cloud_resources)
        self.client = AcsClient(cloud_resources["AccessKey_ID"], cloud_resources["AccessKey_Secret"], cloud_resources["regionID"])

    def init(self):
        self.prepare()
        for ii in range(self.nchunks):
            self.create(ii)

    def create(self, ii):
        '''case1: use existed machine(finished) to make_dispatcher
           case2: create one machine, then make_dispatcher, change status from unallocated to unsubmitted'''
        if self.dispatcher_list[ii]["dispatcher_status"] == "unsubmitted": return # void build twice(in restart situation)
        if self.dispatcher_list[ii]["entity"]: self.make_dispatcher(ii) # case1
        else: # case2
            if len(self.ip_pool) > 0:
                self.dispatcher_list[ii]["entity"] = Entity(self.ip_pool.pop(0), self.server_pool.pop(0))
                self.make_dispatcher(ii)
            else:
                self.server_pool = self.get_server_pool()
                self.ip_pool = self.get_ip(self.server_pool)
                # time.sleep(120)

    # Derivate
    def delete(self, ii):
        '''delete one machine'''
        request = DeleteInstancesRequest()
        request.set_accept_format('json')
        request.set_InstanceIds([self.dispatcher_list[ii]["entity"].instance_id])
        request.set_Force(True)
        response = self.client.do_action_with_exception(request)
        running_num = 0
        for jj in range(self.nchunks):
            if self.dispatcher_list[jj]["dispatcher_status"] == "running" or self.dispatcher_list[jj]["dispatcher_status"] == "unsubmitted":
                running_num += 1
        self.change_apg_capasity(running_num)

    # Derivate
    def catch_dispatcher_exception(self, ii):
        '''everything is okay: return 0
           ssh not active    : return 1
           machine callback  : return 2'''
        if self.check_spot_callback(self.dispatcher_list[ii]["entity"].instance_id): 
            dlog.info("machine %s callback, ip: %s" % (self.dispatcher_list[ii]["entity"].instance_id, self.dispatcher_list[ii]["entity"].ip))
            return 2
        elif not self.dispatcher_list[ii]["dispatcher"].session._check_alive():
            try:
                self.dispatcher_list[ii]["dispatcher"].session.ensure_alive()
                return 0
            except RuntimeError:
                return 1
        else: return 0
        
    def get_server_pool(self):
        running_server = self.describe_apg_instances()
        allocated_server = list(item["entity"].instance_id for item in self.dispatcher_list if item["entity"])
        return list(set(running_server) - set(allocated_server))

    def clean(self):
        self.delete_apg()
        self.delete_template()
        os.remove("apg_id.json")

    def prepare(self):
        restart = False
        if os.path.exists('apg_id.json'):
            with open('apg_id.json') as fp:
                apg = json.load(fp)
                self.cloud_resources["apg_id"] = apg["apg_id"]
            running_server_pool = []
            for ii in range(self.nchunks):
                fn = 'jr.%.06d.json' % ii
                if os.path.exists(os.path.join(os.path.abspath(self.work_path), fn)):
                    task_chunks_str = ['+'.join(ii) for ii in self.task_chunks]
                    task_hashes = [sha1(ii.encode('utf-8')).hexdigest() for ii in task_chunks_str]
                    cur_hash = task_hashes[ii]
                    job_record = JobRecord(self.work_path, self.task_chunks[ii], fn)
                    if not job_record.check_finished(task_hashes[ii]): 
                        if not self.check_spot_callback(job_record.record[cur_hash]['context']['instance_id']):
                            self.dispatcher_list[ii]["entity"] = Entity(job_record.record[cur_hash]['context']['ip'], job_record.record[cur_hash]['context']['instance_id'], job_record)
                            dlog.info("prepare: make_dispatcher")
                            self.make_dispatcher(ii)
                        else:
                            # dlog.info("prepare: callback")
                            # self.dispatcher_list[ii]["dispatcher_status"] == "unallocated"
                            os.remove(os.path.join(os.path.abspath(self.work_path), fn))
                    else:
                        self.dispatcher_list[ii]["dispatcher_status"] = "finished"
            # self.server_pool = list(set(self.get_server_pool()) - set(running_server_pool))
            self.server_pool = self.get_server_pool()
            self.ip_pool = self.get_ip(self.server_pool)
            print(self.ip_pool)
            restart = True
            time.sleep(120)
        img_id = self.get_image_id(self.cloud_resources["img_name"])
        sg_id, vpc_id = self.get_sg_vpc_id()
        self.cloud_resources["template_id"] = self.create_template(img_id, sg_id, vpc_id)
        self.cloud_resources["vsw_id"] = self.get_vsw_id(vpc_id)
        if not restart:
            dlog.info("begin to create apg")
            self.cloud_resources["apg_id"] = self.create_apg()
            time.sleep(150)
            self.server_pool = self.get_server_pool()
            self.ip_pool = self.get_ip(self.server_pool)
        else: dlog.info("restart dpgen")
        
    def delete_apg(self):
        request = DeleteAutoProvisioningGroupRequest()
        request.set_accept_format('json')
        request.set_AutoProvisioningGroupId(self.cloud_resources["apg_id"])
        request.set_TerminateInstances(True)
        try:
            response = self.client.do_action_with_exception(request)
        except ServerException as e:
            dlog.info(e)
        except ClientException as e:
            dlog.info(e)

    def create_apg(self):
        request = CreateAutoProvisioningGroupRequest()
        request.set_accept_format('json')
        request.set_TotalTargetCapacity(str(self.nchunks_limit))
        request.set_LaunchTemplateId(self.cloud_resources["template_id"])
        request.set_AutoProvisioningGroupName(self.cloud_resources["instance_name"] + ''.join(random.choice(string.ascii_uppercase) for _ in range(20)))
        request.set_AutoProvisioningGroupType("maintain")
        request.set_SpotAllocationStrategy("lowest-price")
        request.set_SpotInstanceInterruptionBehavior("terminate")
        request.set_SpotInstancePoolsToUseCount(1)
        request.set_ExcessCapacityTerminationPolicy("termination")
        request.set_TerminateInstances(True)
        request.set_PayAsYouGoTargetCapacity("0")
        request.set_SpotTargetCapacity(str(self.nchunks_limit))
        config = self.generate_config()
        request.set_LaunchTemplateConfigs(config)
        try:
            response = self.client.do_action_with_exception(request)
            response = json.loads(response)
            with open('apg_id.json', 'w') as fp:
                json.dump({'apg_id': response["AutoProvisioningGroupId"]}, fp, indent=4)
            return response["AutoProvisioningGroupId"]
        except ServerException as e:
            dlog.info(e)
        except ClientException as e:
            dlog.info(e)

    def describe_apg_instances(self):
        request = DescribeAutoProvisioningGroupInstancesRequest()
        request.set_accept_format('json')
        request.set_AutoProvisioningGroupId(self.cloud_resources["apg_id"])
        request.set_PageSize(100)
        iteration = self.nchunks // 100
        instance_list = []
        for i in range(iteration + 1):
            request.set_PageNumber(i+1)
            try:
                response = self.client.do_action_with_exception(request)
                response = json.loads(response)
                for ins in response["Instances"]["Instance"]:
                    instance_list.append(ins["InstanceId"])
            except ServerException as e:
                dlog.info(e)
            except ClientException as e:
                dlog.info(e)
        return instance_list
        
    def generate_config(self):
        machine_config = self.cloud_resources["machine_type_price"]
        config = []
        for conf in machine_config:
            for vsw in self.cloud_resources["vsw_id"]:
                tmp = {
                    "InstanceType": conf["machine_type"],
                    "MaxPrice": str(conf["price_limit"] * conf["numb"]),
                    "VSwitchId": vsw,
                    "WeightedCapacity": "1",
                    "Priority": str(conf["priority"])
                }
                config.append(tmp)
        return config

    def create_template(self, image_id, sg_id, vpc_id):
        request = CreateLaunchTemplateRequest()
        request.set_accept_format('json')
        request.set_LaunchTemplateName(''.join(random.choice(string.ascii_uppercase) for _ in range(20)))
        request.set_ImageId(image_id)
        request.set_ImageOwnerAlias("self")
        request.set_PasswordInherit(True)
        request.set_InstanceType("ecs.c6.large")
        request.set_InstanceName(self.cloud_resources["instance_name"])
        request.set_SecurityGroupId(sg_id)
        request.set_VpcId(vpc_id)
        request.set_SystemDiskCategory("cloud_efficiency")
        request.set_SystemDiskSize(40)
        request.set_IoOptimized("optimized")
        request.set_InstanceChargeType("PostPaid")
        request.set_NetworkType("vpc")
        request.set_SpotStrategy("SpotWithPriceLimit")
        request.set_SpotPriceLimit(100)
        try:
            response = self.client.do_action_with_exception(request)
            response = json.loads(response)
            return response["LaunchTemplateId"]
        except ServerException as e:
            dlog.info(e)
        except ClientException as e:
            dlog.info(e)
            
    def delete_template(self):
        request = DeleteLaunchTemplateRequest()
        request.set_accept_format('json')
        request.set_LaunchTemplateId(self.cloud_resources["template_id"])
        response = self.client.do_action_with_exception(request)
        
    def get_image_id(self, img_name):
        request = DescribeImagesRequest()
        request.set_accept_format('json')
        request.set_ImageOwnerAlias("self")
        response = self.client.do_action_with_exception(request)
        response = json.loads(response)
        for img in response["Images"]["Image"]:
            if img["ImageName"] == img_name:
                return img["ImageId"]

    def get_sg_vpc_id(self):
        request = DescribeSecurityGroupsRequest()
        request.set_accept_format('json')
        response = self.client.do_action_with_exception(request)
        response = json.loads(response)
        for sg in response["SecurityGroups"]["SecurityGroup"]:
            if sg["SecurityGroupName"] == "sg":
                return sg["SecurityGroupId"], sg["VpcId"]

    def get_vsw_id(self, vpc_id):
        request = DescribeVpcsRequest()
        request.set_accept_format('json')
        request.set_VpcId(vpc_id)
        response = self.client.do_action_with_exception(request)
        response = json.loads(response)
        for vpc in response["Vpcs"]["Vpc"]:
            if vpc["VpcId"] == vpc_id:
                return vpc["VSwitchIds"]["VSwitchId"]

    def change_apg_capasity(self, capasity):
        request = ModifyAutoProvisioningGroupRequest()
        request.set_accept_format('json')
        request.set_AutoProvisioningGroupId(self.cloud_resources["apg_id"])
        request.set_TotalTargetCapacity(str(capasity))
        request.set_SpotTargetCapacity(str(capasity))
        request.set_PayAsYouGoTargetCapacity("0")
        response = self.client.do_action_with_exception(request)

    def check_spot_callback(self, instance_id):
        request = DescribeInstancesRequest()
        request.set_accept_format('json')
        request.set_InstanceIds([instance_id])
        status = False
        try:
            response = self.client.do_action_with_exception(request)
            response = json.loads(response)
            if len(response["Instances"]["Instance"]) == 1 and "Recycling" in response["Instances"]["Instance"][0]["OperationLocks"]["LockReason"]:
                status = True
            if instance_id not in self.describe_apg_instances():
                status = True
        except ServerException as e:
            dlog.info(e)
        except ClientException as e:
            dlog.info(e)
        return status

    def get_ip(self, instance_list):
        request = DescribeInstancesRequest()
        request.set_accept_format('json')
        ip_list = []
        if len(instance_list) == 0: return ip_list
        try:
            if len(instance_list) <= 10:
                for i in range(len(instance_list)):
                    request.set_InstanceIds([instance_list[i]])
                    response = self.client.do_action_with_exception(request)
                    response = json.loads(response)
                    ip_list.append(response["Instances"]["Instance"][0]["VpcAttributes"]["PrivateIpAddress"]['IpAddress'][0])
                    # ip_list.append(response["Instances"]["Instance"][0]["PublicIpAddress"]["IpAddress"][0])
            else:
                iteration = len(instance_list) // 10
                for i in range(iteration):
                    for j in range(10):
                        request.set_InstanceIds([instance_list[i*10+j]])
                        response = self.client.do_action_with_exception(request)
                        response = json.loads(response)
                        ip_list.append(response["Instances"]["Instance"][0]["VpcAttributes"]["PrivateIpAddress"]['IpAddress'][0])
                        # ip_list.append(response["Instances"]["Instance"][0]["PublicIpAddress"]["IpAddress"][0])
                if len(instance_list) - iteration * 10 != 0:
                    for j in range(len(instance_list) - iteration * 10):
                        request.set_InstanceIds([instance_list[iteration*10+j]])
                        response = self.client.do_action_with_exception(request)
                        response = json.loads(response)
                        ip_list.append(response["Instances"]["Instance"][0]["VpcAttributes"]["PrivateIpAddress"]['IpAddress'][0])
                        # ip_list.append(response["Instances"]["Instance"][0]["PublicIpAddress"]["IpAddress"][0])
            return ip_list
        except: return []


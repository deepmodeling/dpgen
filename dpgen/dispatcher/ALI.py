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
from dpgen.dispatcher.DispatcherList import DispatcherList
from os.path import join
from dpgen  import dlog
from hashlib import sha1

# def manual_delete(stage):
#     with open('machine-ali.json') as fp1:
#         mdata = json.load(fp1)
#         adata = mdata[stage][0]['machine']['ali_auth']
#         mdata_resources = mdata[stage][0]['resources']
#         mdata_machine = mdata[stage][0]['machine']
#         ali = ALI(adata, mdata_resources, mdata_machine, 0)
#         with open('apg_id.json', 'r') as fp2:
#             apg_id = json.load(fp2)
#             ali.apg_id = apg_id['apg_id']
#             ali.delete_apg()

# def manual_create(stage, machine_number):
#     with open('machine-ali.json') as fp:
#         mdata = json.load(fp)
#         adata = mdata[stage][0]['machine']['ali_auth']
#         mdata_resources = mdata[stage][0]['resources']
#         mdata_machine = mdata[stage][0]['machine']
#         ali = ALI(adata, mdata_resources, mdata_machine, machine_number)
#         ali.create_ess()
#         print(ali.ip_list)

# cloud_resources = {"AccessKey_ID":"",
#                    "AccessKey_Secret":"",
#                    "regionID": "cn-shenzhen",
#                    "img_name": "kit",
#                    "machine_type_price": [
#                        {"machine_type": "ecs.gn6v-c8g1.2xlarge", "price_limit": 20.00, "numb": 1, "priority": 0},
#                        {"machine_type": "ecs.gn5-c4g1.xlarge",   "price_limit": 20.00, "numb": 1, "priority": 1}
#                    ],
#                    "instance_name": "CH4", 
#                    "pay_strategy": "spot"
#                    "apg_id": apg_id, 
#                    "template_id": template_id, 
#                    "vsw_id": vsw_id, 
#                    "region_id": region_id,
#                    "client": client}


class ALI(DispatcherList):
    def __init__(self, mdata_machine, mdata_resources, run_tasks, group_size, work_path, cloud_resources=None):
        super().__init__()
        self.cloud_resources["client"] = AcsClient(cloud_resources["AccessKey_ID"], cloud_resources["AccessKey_Secret"], cloud_resources["regionID"])

    def init(self):
        self.prepare()
        for ii in range(self.nchunks):
            self.create(ii)

    def create(self, ii):
        '''case1: jr.json existed and job not finished, use jr.json to rebuild dispatcher
           case2: create one machine, then make_dispatcher, change status from unallocated to unsubmitted
           case3: use existed machine(finish) to make_dispatcher'''
        if not os.path.exists(os.path.join(os.path.abspath(self.work_path), "jr.%.06d.json" % ii)):
            # if entity has value -> case3
            # else create machine -> case2
            # job_record = JobRecord(self.work_path, self.task_chunks[ii], fname = "jr.%.06d.json" % ii)
            # self.dispatcher_list[ii]["entity"] = Entity(ip, instance_id, job_record)
            # self.make_dispatcher(ii)
            if self.dispatcher_list[ii]["entity"].ip:
                self.make_dispatcher(ii)
            else:
                if len(self.server_pool) > 0:
                    self.dispatcher_list[ii]["entity"] = Entity(self.ip_pool.pop(0), self.server_pool.pop(0))
                    self.make_dispatcher(ii)
                else:
                    self.server_pool = self.describe_apg_instances()
                    self.ip_pool = self.get_ip(self.server_pool)
        else:
            # case1
            task_chunks_str = ['+'.join(ii) for ii in self.task_chunks]
            task_hashes = [sha1(ii.encode('utf-8')).hexdigest() for ii in task_chunks_str]
            job_record = JobRecord(self.work_path, self.task_chunks[ii], fname = "jr.%.06d.json" % ii)
            if not job_record.check_finished(task_hashes[ii]): 
                with open(os.path.join(self.work_path, "jr.%.06d.json" % ii)) as fp:
                    jr = json.load(fp)
                    self.dispatcher_list[ii]["entity"] = Entity(jr[cur_hash]['context']['ip'], jr[cur_hash]['context']['instance_id'], job_record)
                    self.make_dispatcher(ii)
    
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
            if self.dispatcher_list[jj]["dispatcher_status"] == "running":
                running_num += 1
        self.change_apg_capasity(running_num)

    def clean(self):
        self.delete_apg()
        self.delete_template()

    # Derivate
    def catch_dispatcher_exception(self):
        pass

    def prepare(self):
        img_id = self.get_image_id(self.cloud_resources["img_name"])
        sg_id, vpc_id = self.get_sg_vpc_id()
        self.cloud_resources["template_id"] = self.create_template(img_id, sg_id, vpc_id)
        self.cloud_resources["vsw_id"] = self.get_vsw_id(vpc_id)
        self.cloud_resources["apg_id"] = self.create_apg()
        dlog.info("begin to create apg, please wait")
        time.sleep(120)
        self.server_pool = self.describe_apg_instances()
        self.ip_pool = self.get_ip(self.server_pool)

    def delete_apg(self):
        request = DeleteAutoProvisioningGroupRequest()
        request.set_accept_format('json')
        request.set_AutoProvisioningGroupId(self.apg_id)
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
        request.set_TotalTargetCapacity(str(self.nchunks))
        request.set_LaunchTemplateId(self.template_id)
        request.set_AutoProvisioningGroupName(self.adata["instance_name"] + ''.join(random.choice(string.ascii_uppercase) for _ in range(20)))
        request.set_AutoProvisioningGroupType("maintain")
        request.set_SpotAllocationStrategy("lowest-price")
        request.set_SpotInstanceInterruptionBehavior("terminate")
        request.set_SpotInstancePoolsToUseCount(1)
        request.set_ExcessCapacityTerminationPolicy("termination")
        request.set_TerminateInstances(True)
        request.set_PayAsYouGoTargetCapacity("0")
        request.set_SpotTargetCapacity(str(self.nchunks))
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
        
    def update_server_list(self):
        instance_list = self.describe_apg_instances()
        return list(set(instance_list) - set(self.instance_list))

    def describe_apg_instances(self):
        request = DescribeAutoProvisioningGroupInstancesRequest()
        request.set_accept_format('json')
        request.set_AutoProvisioningGroupId(self.apg_id)
        request.set_PageSize(100)
        iteration = self.nchunks // 100
        instance_list = []
        for i in range(iteration + 1):
            request.set_PageNumber(i+1)
            response = self.client.do_action_with_exception(request)
            response = json.loads(response)
            for ins in response["Instances"]["Instance"]:
                instance_list.append(ins["InstanceId"])
        return instance_list
        
    def generate_config(self):
        machine_config = self.adata["machine_type_price"]
        config = []
        for conf in machine_config:
            for vsw in self.vsw_id:
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
        request.set_InstanceName(self.adata["instance_name"])
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
        request.set_LaunchTemplateId(self.template_id)
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
        request.set_AutoProvisioningGroupId(self.apg_id)
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
        except:
            pass
        return status
































from dpgen.dispatcher.ALI import ALI
import unittest

class TestALI(unittest.TestCase):
    def setUp(self):
        mdata_machine = {
            "batch": "shell",
            "hostname": "",
            "password": "975481DING!",
            "port": 22,
            "username": "root",
            "work_path": "/root/dpgen_work",
        }
        mdata_resources = {
            "allow_failure": True,
            "ratio_failue": 0.05,
            "task_per_node": 16,
            "with_mpi": True,
            "source_list": ["/opt/intel/parallel_studio_xe_2018/psxevars.sh"],
            "envs" : {"PATH" : "/root/deepmd-pkg/vasp.5.4.4/bin:$PATH"}
        }
        run_tasks = ["000", "001", "002", "003"]
        group_size = 1
        work_path = ""
        cloud_resources = {
            "cloud_platform": "ali",
            "AccessKey_ID":"LTAI4Fo5VhHMJUcEUgazCRMf",
            "AccessKey_Secret":"5rPSLRplhPI2UJycvHS9wsTY0NnBmT",
            "regionID": "cn-beijing",
            "img_name": "vasp",
            "machine_type_price": [
                {"machine_type": "ecs.c5.4xlarge", "price_limit": 0.05, "numb": 16, "priority": 0},
                {"machine_type": "ecs.c6.4xlarge", "price_limit": 0.05, "numb": 16, "priority": 1}
            ],
            "instance_name": "CH4_test_dingzhaohan", 
            "pay_strategy": "spot"
        }
        ali = ALI(mdata_machine, mdata_resources, run_tasks, group_size, work_path, cloud_resources)

    def tearDown(self):
        pass

    def test_create():
        ali.create()
        pass

    def test_delete():
        pass

    def test_get_server_pool():
        pass
    
    def test_clean():
        pass

    def test_prepare():
        pass

    def test_create_apg():
        pass

    def test_delete_apg():
        pass

    def test_describe_apg_instances():
        pass

    def test_generate_config():
        pass

    def test_create_template():
        pass

    def test_delete_template():
        pass

    def test_get_image_id():
        pass

    def test_get_sg_vpc_id():
        pass

    def test_get_vsw_id():

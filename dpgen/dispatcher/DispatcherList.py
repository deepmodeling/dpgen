from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, JobRecord
from paramiko.ssh_exception import NoValidConnectionsError
import os, time
from dpgen import dlog
class Entity():
    def __init__(self, ip, instance_id, job_record=None, job_handler=None):
        self.ip = ip
        self.instance_id = instance_id
        self.job_record = job_record
        self.job_handler = job_handler

class DispatcherList():
    def __init__(self, mdata_machine, mdata_resources, work_path, run_tasks, group_size, cloud_resources=None):
        self.mdata_machine = mdata_machine
        self.mdata_resources = mdata_resources
        self.task_chunks = _split_tasks(run_tasks, group_size)
        self.nchunks = len(self.task_chunks)
        self.nchunks_limit = int(self.mdata_machine.get("machine_upper_bound", self.nchunks))
        if(self.nchunks_limit > self.nchunks):
            self.nchunks_limit = self.nchunks    
        self.work_path = work_path
        self.cloud_resources = cloud_resources
        self.server_pool = []
        self.ip_pool = []
        self.dispatcher_list = list({"dispatcher": None, 
                                     "dispatcher_status": "unallocated",
                                     "entity": None} for ii in range(self.nchunks))
    # Derivate
    def init(self):
        # do something necessary
        for ii in range(self.nchunks):
            self.create(ii)

    # Base
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
                 mark_failure = False,
                 outlog = 'log',
                 errlog = 'err'):
        ratio_failure = self.mdata_resources.get("ratio_failure", 0)
        while True:
            if self.check_all_dispatchers_finished(ratio_failure):
                self.clean()
                break
            self.exception_handling(ratio_failure)
            jj = self.nchunks - 1
            for ii in range(self.nchunks):
                dispatcher_status = self.check_dispatcher_status(ii)
                if dispatcher_status == "unsubmitted":
                    dlog.info(self.dispatcher_list[ii]["entity"].ip)
                    self.dispatcher_list[ii]["entity"].job_handler = self.dispatcher_list[ii]["dispatcher"].submit_jobs(resources,
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
                    self.dispatcher_list[ii]["entity"].job_record = self.dispatcher_list[ii]["entity"].job_handler["job_record"]
                    self.dispatcher_list[ii]["dispatcher_status"] = "running"
                elif dispatcher_status == "finished" and self.dispatcher_list[ii]["entity"]: 
                    # no jobs in queue, delete current machine
                    # else add current machine to server_pool
                    entity = self.dispatcher_list[ii]["entity"]
                    status_list = [item["dispatcher_status"] for item in self.dispatcher_list]
                    flag = "unallocated" in status_list
                    if not flag:
                        self.delete(ii)
                        self.dispatcher_list[ii]["entity"] = None
                    else:
                        self.dispatcher_list[ii]["entity"] = None
                        self.server_pool.append(entity.instance_id)
                        self.ip_pool.append(entity.ip)
                        while(jj>=ii):
                            if(self.dispatcher_list[jj]["dispatcher_status"] == "unallocated"):
                                self.create(jj)
                                if(self.dispatcher_list[jj]["dispatcher_status"] == "unsubmitted"):
                                    dlog.info(self.dispatcher_list[jj]["entity"].ip)
                                    self.dispatcher_list[jj]["entity"].job_handler = self.dispatcher_list[jj]["dispatcher"].submit_jobs(resources,
                                                                                                                                        command,
                                                                                                                                        work_path,
                                                                                                                                        self.task_chunks[jj],
                                                                                                                                        group_size,
                                                                                                                                        forward_common_files,
                                                                                                                                        forward_task_files,
                                                                                                                                        backward_task_files,
                                                                                                                                        forward_task_deference,
                                                                                                                                        outlog,
                                                                                                                                        errlog)
                                    self.dispatcher_list[jj]["entity"].job_record = self.dispatcher_list[jj]["entity"].job_handler["job_record"]
                                    self.dispatcher_list[jj]["dispatcher_status"] = "running"
                                break
                            jj -=1
                elif dispatcher_status == "running":
                    pass
                elif dispatcher_status == "unallocated":
                    # if len(server_pool) > 0: make_dispatcher
                    # else: pass
                    self.create(ii)
                    if self.dispatcher_list[ii]["dispatcher_status"] == "unsubmitted":
                        dlog.info(self.dispatcher_list[ii]["entity"].ip)
                        self.dispatcher_list[ii]["entity"].job_handler = self.dispatcher_list[ii]["dispatcher"].submit_jobs(resources,
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
                        self.dispatcher_list[ii]["entity"].job_record = self.dispatcher_list[ii]["entity"].job_handler["job_record"]
                        self.dispatcher_list[ii]["dispatcher_status"] = "running"
                elif dispatcher_status == "terminated":
                    pass
            self.update()
            time.sleep(10)    

    # Derivate
    def create(self, ii):
        '''case1: use existed machine(finished) to make_dispatcher
           case2: create one machine, then make_dispatcher, change status from unallocated to unsubmitted'''
        pass
    
    # Derivate
    def delete(self, ii):
        '''delete one machine
            if entity is none, means this machine is used by another dispatcher, shouldn't be deleted'''
        pass

    # Derivate, delete config like templates, etc.
    def clean(self):
        pass

    # Derivate
    def update():
        pass

    # Base
    def check_all_dispatchers_finished(self, ratio_failure):    
        status_list = [item["dispatcher_status"] for item in self.dispatcher_list]
        finished_num = status_list.count("finished")
        if finished_num / self.nchunks < (1 - ratio_failure): return False
        else: return True

    # Base
    def exception_handling(self, ratio_failure):
        status_list = [item["dispatcher_status"] for item in self.dispatcher_list]
        terminated_num = status_list.count("terminated")
        if terminated_num / self.nchunks > ratio_failure:
            # self.dispatcher_list = [lambda item["dispatcher_status"]: "finished" for item in self.dispatcher_list if item["dispatcher_status"] == "terminated"]
            for ii in range(self.nchunks):
                if self.dispatcher_list[ii]["dispatcher_status"] == "terminated":
                    self.dispatcher_list[ii]["dispatcher_status"] = "unallocated"
    # Base
    def make_dispatcher(self, ii):
        entity = self.dispatcher_list[ii]["entity"]
        profile = self.mdata_machine.copy()
        profile['hostname'] = entity.ip
        profile['instance_id'] = entity.instance_id
        count = 0
        flag = 0
        while count < 3:
            try:
                self.dispatcher_list[ii]["dispatcher"] = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
                self.dispatcher_list[ii]["dispatcher_status"] = "unsubmitted"
                flag = 1
                break
            except Exception:
                count += 1
                time.sleep(60)
        if not flag:
            # give up this machine, wait other machine in sever_pool.
            # this machine will be append into server_pool next time when update apg_instances.
            self.dispatcher_list[ii]["entity"] = None


    # Base
    def check_dispatcher_status(self, ii, allow_failure=False):
        '''catch running dispatcher exception
           if no exception occured, check finished'''
        if self.dispatcher_list[ii]["dispatcher_status"] == "running":
            status = self.catch_dispatcher_exception(ii)
            if status == 0:
                # param clean: delete remote work_dir or not.
                clean = self.mdata_resources.get("clean", False)
                try:
                    # avoid raising ssh exception in download proceess
                    finished = self.dispatcher_list[ii]["dispatcher"].all_finished(self.dispatcher_list[ii]["entity"].job_handler, allow_failure, clean)
                    if finished:
                        self.dispatcher_list[ii]["dispatcher_status"] = "finished"
                except Exception:
                    pass                    
            elif status == 1:
                # self.dispatcher_list[ii]["dispatcher_status"] = "terminated"
                pass
            elif status == 2:
                self.dispatcher_list[ii]["dispatcher"] = None
                self.dispatcher_list[ii]["dispatcher_status"] = "terminated"
                self.dispatcher_list[ii]["entity"] = None
                os.remove(os.path.join(self.work_path, "jr.%.06d.json" % ii))
        return self.dispatcher_list[ii]["dispatcher_status"]

    # Derivate
    def catch_dispatcher_exception(self, ii):
        '''everything is okay: return 0
           ssh not active    : return 1
           machine callback  : return 2'''
        pass







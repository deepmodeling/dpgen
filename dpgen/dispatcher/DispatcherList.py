from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, JobRecord

class Entity():
    def __init__(self, ip, instance_id, job_handler):
        self.ip = ip
        self.instance_id = instance_id
        self.job_handler = job_handler
        
class DispatcherList():
    def __init__(self, mdata_machine, mdata_resources, nchunks, cloud_resources=None):
        self.mdata_machine = mdata_machine
        self.mdata_resources = mdata_resources
        self.nchunks = nchunks
        self.cloud_resources = cloud_resources
        self.dispatcher_list = list({"dispatcher": None, 
                                     "dispatcher_status": "unallocated",
                                     "entity": None} for ii in range(nchunks))
    # Base
    def init(self):
        for dispatcher in self.dispatcher_list:
            self.create()

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
        task_chunks = _split_tasks(tasks, group_size)
        while True:
            if self.check_all_dispatchers_finished():
                break
            ratio_failure = self.mdata_resources["ratio_failure"]
            self.exception_handling(ratio_failure)
            for ii in range(self.nchunks):
                dispatcher_status = self.check_dispatcher_status(ii)
                if dispatcher_status == "unsubmitted":
                    self.dispatcher_list[ii]["job_handler"] = self.dispatcher_list[ii]["dispatcher"].submit_jobs(resources,
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
                    self.dispatcher_list[ii]["dispatcher_status"] == "running"
                elif dispatcher_status == "finished":
                    # to do 
                    # no jobs in queue, delete current machine
                    # else allocate current machine to unalloc dispatcher
                elif dispatcher_status == "running":
                    pass
                elif dispatcher_status == "unallocated":
                    # to do: if we can create machine, then make dispatcher
                    #        else pass
                elif dispatcher_status == "terminated":
                    pass

    # Derivate
    def check_restart(self):
        '''use some methods to decide restart or not'''
        pass

    # Derivate
    def create(self, num):
        '''create machines
           assignment the entities
           make dispatchers
           change status from unallocated to unsubmitted'''
        pass

    # Derivate
    def delete(self):
        '''delete one machine'''
        pass

    # Base
    def check_all_dispatchers_finished(self, ratio_failure=0):
        exception_num = 0
        for ii in range(self.nchunks):
            if self.dispatcher_list[ii]["dispatcher_status"] == "terminated":
                exception_num += 1
        if exception_num / self.nchunks > ratio_failure: return False
        else return True

    # Base
    def exception_handling(self, ratio_failure):
        terminated_num = 0
        for ii in range(self.nchunks):
            if self.dispatcher_list[ii]["dispatcher_status"] == "terminated":
                terminated_num += 1
        if terminated_num / self.nchunks > ratio_failure:
            self.resubmit()

    # Derivate
    def resubmit(self):
        '''create machines
           make dispatcher
           change status from terminated to unsubmitted'''
        pass

    # Base
    def make_dispatcher(self, entity):
        '''use entity to distinguish machine, for example if ip isn't None, means we can make_dispatcher
           change status from unallocated to unsubmitted'''
        if entity["ip"]:
            profile = self.mdata_machine.copy()
            profile['hostname'] = entity["ip"]
            profile['instance_id'] = entity["instance_id"]
            dispatcher = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
            self.dispatcher_list[ii]["dispatcher_status"] = "unsubmitted"

    # Base
    def check_dispatcher_status(self, ii):
        '''catch running dispatcher exception
           if no exception occured, check finished'''
        if self.dispatcher_list[ii]["dispatcher_status"] == "running":
            if not self.catch_dispatcher_exception():
                if self.dispatcher_list[ii]["dispatcher"].all_finished():
                    self.dispatcher_list[ii]["dispatcher_status"] = "finished"
            else:
                self.dispatcher_list[ii]["dispatcher_status"] = "terminated"
        return self.dispatcher_list[ii]["dispatcher_status"]

    # Derivate
    def catch_dispatcher_exception(self):
        pass







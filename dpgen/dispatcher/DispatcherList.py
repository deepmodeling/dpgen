from dpgen.dispatcher.Dispatcher import Dispatcher, _split_tasks, JobRecord

class Entity():
    def __init__(self, ip, instance_id, job_record=None, job_handler=None):
        self.ip = ip
        self.instance_id = instance_id
        self.job_record = job_record
        self.job_handler = job_handler

class DispatcherList():
    def __init__(self, mdata_machine, mdata_resources, run_tasks, group_size, work_path, cloud_resources=None):
        self.mdata_machine = mdata_machine
        self.mdata_resources = mdata_resources
        self.task_chunks = _split_tasks(run_tasks, group_size)
        self.nchunks = len(self.task_chunks)
        self.work_path = work_path
        self.cloud_resources = cloud_resources
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
        ratio_failure = self.mdata_resources.get("ratio_failue", 0)
        while True:
            if self.check_all_dispatchers_finished(ratio_failue):
                self.clean()
                break
            self.exception_handling(ratio_failure)
            for ii in range(self.nchunks):
                dispatcher_status = self.check_dispatcher_status(ii)
                if dispatcher_status == "unsubmitted":
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
                elif dispatcher_status == "finished":
                    # to do 
                    # no jobs in queue, delete current machine
                    # else allocate current machine to unalloc dispatcher
                    flag = True
                    for jj in range(ii, self.nchunks):
                        if self.dispatcher_list[jj]["dispatcher_status"] == "unallocated":
                            flag = False
                            entity = self.dispatcher_list[ii]["entity"]
                            self.dispatcher_list[jj]["entity"] = Entity(entity.ip, entity.instance_id)
                            self.create(ii)
                            break
                    if flag: self.delete(ii)
                elif dispatcher_status == "running":
                    pass
                elif dispatcher_status == "unallocated":
                    # to do: if we can create machine, then make dispatcher
                    #        else pass
                    self.create(ii)
                elif dispatcher_status == "terminated":
                    pass

    # Derivate
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
            pass
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
        pass

    # Derivate, delete config like templates, etc.
    def clean(self):
        pass

    # Base
    def check_all_dispatchers_finished(self, ratio_failure):
        finished_num = 0
        for ii in range(self.nchunks):
            if self.dispatcher_list[ii]["dispatcher_status"] == "finished":
                finished_num += 1
        if finished_num / self.nchunks < (1 - ratio_failure): return False
        else: return True

    # Base
    def exception_handling(self, ratio_failure):
        terminated_num = 0
        for ii in range(self.nchunks):
            if self.dispatcher_list[ii]["dispatcher_status"] == "terminated":
                terminated_num += 1
                if terminated_num / self.nchunks > ratio_failure:

                    # os.remove(os.path.join(self.work_path, "jr.%.06d.json" % ii))
                    self.create(ii)

    # Base
    def make_dispatcher(self, ii):
        '''use entity to distinguish machine, for example if ip isn't None, means we can make_dispatcher
           change status from unallocated to unsubmitted'''
        entity = self.dispatcher_list[ii]["entity"]
        if entity.ip:
            profile = self.mdata_machine.copy()
            profile['hostname'] = entity.ip
            profile['instance_id'] = entity.instance_id
            self.dispatcher_list[ii]["dispatcher"] = Dispatcher(profile, context_type='ssh', batch_type='shell', job_record='jr.%.06d.json' % ii)
            self.dispatcher_list[ii]["dispatcher_status"] = "unsubmitted"

    # Base
    def check_dispatcher_status(self, ii, allow_failue=False):
        '''catch running dispatcher exception
           if no exception occured, check finished'''
        if self.dispatcher_list[ii]["dispatcher_status"] == "running":
            if self.catch_dispatcher_exception(ii) == 0:
                # param clean: delete remote work_dir or not.
                clean = self.mdata_resources.get("clean", False)
                if self.dispatcher_list[ii]["dispatcher"].all_finished(self.dispatcher_list[ii]["entity"].job_handler, mark_failure, clean):
                    self.dispatcher_list[ii]["dispatcher_status"] = "finished"
            elif self.catch_dispatcher_exception(ii) == 1:
                # self.dispatcher_list[ii]["dispatcher_status"] = "terminated"
                pass
            elif self.catch_dispatcher_exception(ii) == 2:
                self.dispatcher_list[ii]["dispatcher_status"] = "terminated"
                os.remove(os.path.join(self.work_path, "jr.%.06d.json" % ii))
        return self.dispatcher_list[ii]["dispatcher_status"]

    # Derivate
    def catch_dispatcher_exception(self, ii):
        '''everything is okay: return 0
           ssh not active    : return 1
           machine callback  : return 2'''
        pass







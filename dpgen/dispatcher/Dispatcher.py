import os,sys,time,random,json,glob
from dpgen.dispatcher.LocalContext import LocalSession
from dpgen.dispatcher.LocalContext import LocalContext
from dpgen.dispatcher.LazyLocalContext import LazyLocalContext
from dpgen.dispatcher.SSHContext import SSHSession
from dpgen.dispatcher.SSHContext import SSHContext
from dpgen.dispatcher.Slurm import Slurm
from dpgen.dispatcher.LSF import LSF
from dpgen.dispatcher.PBS import PBS
from dpgen.dispatcher.Shell import Shell
from dpgen.dispatcher.AWS import AWS 
from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog
from hashlib import sha1

def _split_tasks(tasks,
                 group_size):
    ntasks = len(tasks)
    ngroups = ntasks // group_size
    if ngroups * group_size < ntasks:
        ngroups += 1
    chunks = [[]] * ngroups
    tot = 0
    for ii in range(ngroups) :
        chunks[ii] = (tasks[ii::ngroups])
        tot += len(chunks[ii])
    assert(tot == len(tasks))
    return chunks

    
class Dispatcher(object):
    def __init__ (self,
                  remote_profile,
                  context_type = 'local',
                  batch_type = 'slurm', 
                  job_record = 'jr.json'):
        self.remote_profile = remote_profile
        if context_type == 'local':
            self.session = LocalSession(remote_profile)
            self.context = LocalContext
            self.uuid_names = False
        elif context_type == 'lazy-local':
            self.session = None
            self.context = LazyLocalContext
            self.uuid_names = True
        elif context_type == 'ssh':
            self.session = SSHSession(remote_profile)
            self.context = SSHContext
            self.uuid_names = False
        else :
            raise RuntimeError('unknown context')
        if batch_type == 'slurm':
            self.batch = Slurm            
        elif batch_type == 'lsf':
            self.batch = LSF
        elif batch_type == 'pbs':
            self.batch = PBS
        elif batch_type == 'shell':
            self.batch = Shell
        elif batch_type == 'aws':
            self.batch = AWS
        else :
            raise RuntimeError('unknown batch ' + batch_type)
        self.jrname = job_record

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
                 errlog = 'err') :
        job_handler = self.submit_jobs(resources,
                                       command,
                                       work_path,
                                       tasks,
                                       group_size,
                                       forward_common_files,
                                       forward_task_files,
                                       backward_task_files,
                                       forward_task_deference,
                                       outlog,
                                       errlog)
        while not self.all_finished(job_handler) :
            time.sleep(10)
        # delete path map file when job finish
        # _pmap.delete()


    def submit_jobs(self,
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
                    errlog = 'err') :
        self.backward_task_files = backward_task_files
        # task_chunks = [
        #     [os.path.basename(j) for j in tasks[i:i + group_size]] \
        #     for i in range(0, len(tasks), group_size)
        # ]
        task_chunks = _split_tasks(tasks, group_size)    
        task_chunks_str = ['+'.join(ii) for ii in task_chunks]
        task_hashes = [sha1(ii.encode('utf-8')).hexdigest() for ii in task_chunks_str]
        job_record = JobRecord(work_path, task_chunks, fname = self.jrname)
        nchunks = len(task_chunks)
        
        job_list = []
        for ii in range(nchunks) :            
            cur_chunk = task_chunks[ii]
            cur_hash = task_hashes[ii]
            if not job_record.check_finished(cur_hash):                
                # chunk is not finished
                # check if chunk is submitted
                submitted = job_record.check_submitted(cur_hash)
                if not submitted:
                    job_uuid = None
                else :
                    job_uuid = job_record.get_uuid(cur_hash)
                    dlog.debug("load uuid %s for chunk %s" % (job_uuid, cur_hash))
                # communication context, bach system
                context = self.context(work_path, self.session, job_uuid)
                batch = self.batch(context, uuid_names = self.uuid_names)
                rjob = {'context':context, 'batch':batch}
                # upload files
                if not rjob['context'].check_file_exists('tag_upload'):
                    rjob['context'].upload('.',
                                           forward_common_files)
                    rjob['context'].upload(cur_chunk,
                                           forward_task_files, 
                                           dereference = forward_task_deference)
                    rjob['context'].write_file('tag_upload', '')
                    dlog.debug('uploaded files for %s' % task_chunks_str[ii])
                # submit new or recover old submission
                if not submitted:
                    rjob['batch'].submit(cur_chunk, command, res = resources, outlog=outlog, errlog=errlog)
                    job_uuid = rjob['context'].job_uuid
                    dlog.debug('assigned uuid %s for %s ' % (job_uuid, task_chunks_str[ii]))
                    dlog.info('new submission of %s for chunk %s' % (job_uuid, cur_hash))
                else:
                    rjob['batch'].submit(cur_chunk, command, res = resources, outlog=outlog, errlog=errlog, restart = True)
                    dlog.info('restart from old submission %s for chunk %s' % (job_uuid, cur_hash))
                # record job and its remote context
                job_list.append(rjob)
                job_record.record_remote_context(cur_hash,                                                 
                                                 context.local_root, 
                                                 context.remote_root, 
                                                 job_uuid)
            else :
                # finished job, append a None to list
                job_list.append(None)
        job_record.dump()
        assert(len(job_list) == nchunks)
        job_handler = {
            'task_chunks': task_chunks,
            'job_list': job_list,
            'job_record': job_record,
            'command': command,
            'resources': resources,
            'outlog': outlog,
            'errlog': errlog,
            'backward_task_files': backward_task_files
        }
        return job_handler


    def all_finished(self, 
                     job_handler):
        task_chunks = job_handler['task_chunks']
        task_chunks_str = ['+'.join(ii) for ii in task_chunks]
        task_hashes = [sha1(ii.encode('utf-8')).hexdigest() for ii in task_chunks_str]
        job_list = job_handler['job_list']
        job_record = job_handler['job_record']
        command = job_handler['command']
        resources = job_handler['resources']
        outlog = job_handler['outlog']
        errlog = job_handler['errlog']
        backward_task_files = job_handler['backward_task_files']
        dlog.debug('checking jobs')
        nchunks = len(task_chunks)
        for idx in range(nchunks) :
            cur_hash = task_hashes[idx]
            rjob = job_list[idx]
            if not job_record.check_finished(cur_hash) :
                # chunk not finished according to record
                status = rjob['batch'].check_status()
                job_uuid = rjob['context'].job_uuid
                dlog.debug('checked job %s' % job_uuid)
                if status == JobStatus.terminated :
                    job_record.increase_nfail(cur_hash)
                    if job_record.check_nfail(cur_hash) > 3:
                        raise RuntimeError('Job %s failed for more than 3 times' % job_uuid)
                    dlog.info('job %s terminated, submit again'% job_uuid)
                    dlog.debug('try %s times for %s'% (job_record.check_nfail(cur_hash), job_uuid))
                    rjob['batch'].submit(task_chunks[idx], command, res = resources, outlog=outlog, errlog=errlog,restart=True)
                elif status == JobStatus.finished :
                    dlog.info('job %s finished' % job_uuid)
                    rjob['context'].download(task_chunks[idx], backward_task_files)
                    rjob['context'].clean()
                    job_record.record_finish(cur_hash)
                    job_record.dump()
        job_record.dump()
        return job_record.check_all_finished()


class JobRecord(object):
    def __init__ (self, path, task_chunks, fname = 'job_record.json', ip=None):
        self.path = os.path.abspath(path)
        self.fname = os.path.join(self.path, fname)
        self.task_chunks = task_chunks
        if not os.path.exists(self.fname):
            self._new_record()
        else :
            self.load()

    def check_submitted(self, chunk_hash):
        self.valid_hash(chunk_hash)
        return self.record[chunk_hash]['context'] is not None

    def record_remote_context(self, 
                              chunk_hash, 
                              local_root, 
                              remote_root, 
                              job_uuid):
        self.valid_hash(chunk_hash)
        self.record[chunk_hash]['context'] = [local_root, remote_root, job_uuid]

    def get_uuid(self, chunk_hash):
        self.valid_hash(chunk_hash)
        return self.record[chunk_hash]['context'][2]

    def check_finished(self, chunk_hash):
        self.valid_hash(chunk_hash)
        return self.record[chunk_hash]['finished']

    def check_all_finished(self):
        flist = [self.record[ii]['finished'] for ii in self.record]
        return all(flist)

    def record_finish(self, chunk_hash):
        self.valid_hash(chunk_hash)
        self.record[chunk_hash]['finished'] = True

    def check_nfail(self,chunk_hash):
        self.valid_hash(chunk_hash)
        return self.record[chunk_hash]['fail_count']

    def increase_nfail(self,chunk_hash):
        self.valid_hash(chunk_hash)
        self.record[chunk_hash]['fail_count'] += 1

    def valid_hash(self, chunk_hash):
        if chunk_hash not in self.record.keys():
            raise RuntimeError('chunk hash %s not in record, a invalid record may be used, please check file %s' % (chunk_hash, self.fname))

    def dump(self):
        with open(self.fname, 'w') as fp:
            json.dump(self.record, fp, indent=4)

    def load(self):
        with open(self.fname) as fp:
            self.record = json.load(fp)

    def _new_record(self):
        task_chunks_str=['+'.join(ii) for ii in self.task_chunks]
        task_hash = [sha1(ii.encode('utf-8')).hexdigest() for ii in task_chunks_str]
        self.record = {}
        for ii,jj in zip(task_hash, self.task_chunks):
            self.record[ii] = {
                'context': None,
                'finished': False,
                'fail_count': 0,
                'task_chunk': jj,
            }


def make_dispatcher(mdata_machine, mdata_resource, run_tasks, group_size):
    if 'ali_auth' in mdata_machine:
        from dpgen.dispatcher.ALI import ALI
        nchunks = len(_split_tasks(run_tasks, group_size))
        dispatcher = ALI(mdata_machine['ali_auth'], mdata_resource, mdata_machine, nchunks)
        dispatcher.init()
        return dispatcher
    else:    
        try:
            hostname = mdata['hostname']
            context_type = 'ssh'
        except:
            context_type = 'local'
        try:
            batch_type = mdata['batch']
        except:
            dlog.info('cannot find key "batch" in machine file, try to use deprecated key "machine_type"')
            batch_type = mdata['machine_type']
        try:
            lazy_local = mdata['lazy_local']
        except:
            lazy_local = False
        if lazy_local and context_type == 'local':
            dlog.info('Dispatcher switches to the lazy local mode')
            context_type = 'lazy-local'
        disp = Dispatcher(mdata, context_type=context_type, batch_type=batch_type)
        return disp

import os,sys,time,random

from dpgen.dispatcher.LocalContext import LocalSession
from dpgen.dispatcher.LocalContext import LocalContext
from dpgen.dispatcher.SSHContext import SSHSession
from dpgen.dispatcher.SSHContext import SSHContext
from dpgen.dispatcher.Slurm import Slurm
from dpgen.dispatcher.Shell import Shell
from dpgen.dispatcher.JobStatus import JobStatus
from dpgen import dlog
from hashlib import sha1
from monty.serialization import dumpfn,loadfn


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
                  batch_type = 'slurm'):
        self.remote_profile = remote_profile
        if context_type == 'local':
            self.session = LocalSession(remote_profile)
            self.context = LocalContext
        elif context_type == 'ssh':
            self.session = SSHSession(remote_profile)
            self.context = SSHContext
        else :
            raise RuntimeError('unknown context')
        if batch_type == 'slurm':
            self.batch = Slurm            
        elif batch_type == 'shell':
            self.batch = Shell
        else :
            raise RuntimeError('unknown batch ' + batch_type)


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
        # task_chunks = [
        #     [os.path.basename(j) for j in tasks[i:i + group_size]] \
        #     for i in range(0, len(tasks), group_size)
        # ]
        task_chunks = _split_tasks(tasks, group_size)    
        _pmap=PMap(work_path)
        path_map=_pmap.load()
        _fr = FinRecord(work_path, len(task_chunks))        

        job_list = []
        task_chunks_=['+'.join(ii) for ii in task_chunks]
        job_fin = _fr.get_record()
        assert(len(job_fin) == len(task_chunks))
        for ii,chunk in enumerate(task_chunks) :
            if not job_fin[ii] :
                # map chunk info. to uniq id    
                chunk_sha1 = sha1(task_chunks_[ii].encode('utf-8')).hexdigest() 
                # if hash in map, recover job, else start a new job
                if chunk_sha1 in path_map:
                    job_uuid = path_map[chunk_sha1][1].split('/')[-1]
                    dlog.debug("load uuid %s for chunk %s" % (job_uuid, task_chunks_[ii]))
                else:
                    job_uuid = None
                # communication context, bach system
                context = self.context(self.session, work_path, job_uuid)
                batch = self.batch(context)
                rjob = {'context':context, 'batch':batch}
                # upload files
                if not rjob['context'].check_file_exists('tag_upload'):
                    rjob['context'].upload('.',
                                           forward_common_files)
                    rjob['context'].upload(chunk,
                                           forward_task_files, 
                                           dereference = forward_task_deference)
                    rjob['context'].write_file('tag_upload', '')
                    dlog.debug('uploaded files for %s' % task_chunks_[ii])
                # submit new or recover old submission
                if job_uuid is None:
                    rjob['batch'].submit(chunk, command, res = resources, outlog=outlog, errlog=errlog)
                    dlog.debug('assigned uudi %s for %s ' % (rjob['context'].job_uuid, task_chunks_[ii]))
                    dlog.info('new submission of %s' % rjob['context'].job_uuid)
                else:
                    rjob['batch'].submit(chunk, command, res = resources, outlog=outlog, errlog=errlog, restart = True)
                    dlog.info('restart from old submission %s ' % job_uuid)
                # record job and its hash
                job_list.append(rjob)
                path_map[chunk_sha1] = [context.local_root,context.remote_root]
            else :
                # finished job, append a None to list
                job_list.append(None)
        _pmap.dump(path_map)

        assert(len(job_list) == len(task_chunks))
        fcount = [0]*len(job_list)
        while not all(job_fin) :
            dlog.debug('checking jobs')
            for idx,rjob in enumerate(job_list) :
                if not job_fin[idx] :
                    status = rjob['batch'].check_status()
                    job_uuid = rjob['context'].job_uuid
                    if status == JobStatus.terminated :
                        fcount[idx] += 1
                        if fcount[idx] > 3:
                            raise RuntimeError('Job %s failed for more than 3 times' % job_uuid)
                        dlog.info('job %s terminated, submit again'% job_uuid)
                        dlog.debug('try %s times for %s'% (fcount[idx], job_uuid))
                        rjob['batch'].submit(task_chunks[idx], command, res = resources, outlog=outlog, errlog=errlog,restart=True)
                    elif status == JobStatus.finished :
                        dlog.info('job %s finished' % job_uuid)
                        rjob['context'].download(task_chunks[idx], backward_task_files)
                        rjob['context'].clean()
                        job_fin[idx] = True
                        _fr.write_record(job_fin)
            time.sleep(10)
        # delete path map file when job finish
        _pmap.delete()


class FinRecord(object):
    def __init__ (self, path, njobs, fname = 'fin.record'):
        self.path = os.path.abspath(path)
        self.fname = os.path.join(self.path, fname)
        self.njobs = njobs

    def get_record(self):
        if not os.path.exists(self.fname):
            return [False] * self.njobs
        else :
            with open(self.fname) as fp:
                return [bool(int(ii)) for ii in fp.read().split()]

    def write_record(self, job_fin):
        with open(self.fname, 'w') as fp:
            for ii in job_fin:
                if ii:
                    fp.write('1 ')
                else:
                    fp.write('0 ')


class PMap(object):
   '''
   Path map class to operate {read,write,delte} the pmap.json file
   '''

   def __init__(self,path,fname="pmap.json"):
       self.f_path_map=os.path.join(path,fname)

   def load(self):
      f_path_map=self.f_path_map
      if os.path.isfile(f_path_map):
         path_map=loadfn(f_path_map)
      else:
         path_map={}
      return path_map

   def dump(self,pmap,indent=4):
      f_path_map=self.f_path_map
      dumpfn(pmap,f_path_map,indent=indent)

   def delete(self):
      f_path_map=self.f_path_map
      try:
         os.remove(f_path_map)
      except:
         pass

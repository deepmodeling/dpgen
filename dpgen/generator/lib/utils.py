#!/usr/bin/env python3

import os, re, shutil, logging
import glob

iter_format = "%06d"
task_format = "%02d"
log_iter_head = "iter " + iter_format + " task " + task_format + ": "

def make_iter_name (iter_index) :
    return "iter." + (iter_format % iter_index)

def create_path (path) :
    path += '/'
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)        
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
            counter += 1
    os.makedirs (path)

def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

def copy_file_list (file_list, from_path, to_path) :
    for jj in file_list : 
        if os.path.isfile(os.path.join(from_path, jj)) :
            shutil.copy (os.path.join(from_path, jj), to_path)
        elif os.path.isdir(os.path.join(from_path, jj)) :
            shutil.copytree (os.path.join(from_path, jj), os.path.join(to_path, jj))

def cmd_append_log (cmd,
                    log_file) :
    ret = cmd
    ret = ret + " 1> " + log_file
    ret = ret + " 2> " + log_file
    return ret

def log_iter (task, ii, jj) :
    logging.info ((log_iter_head + "%s") % (ii, jj, task))

def repeat_to_length(string_to_expand, length):
    ret = ""
    for ii in range (length) : 
        ret += string_to_expand
    return ret

def log_task (message) :
    header = repeat_to_length (" ", len(log_iter_head % (0, 0)))
    logging.info (header + message)

def record_iter (record, ii, jj) :
    with open (record, "a") as frec :
        frec.write ("%d %d\n" % (ii, jj))

def symlink_user_forward_files(mdata, task_type, work_path, task_format = None):
    '''
    Symlink user-defined forward_common_files
    Current path should be work_path, such as 00.train
    
    Parameters
    ---------
    mdata : dict
        machine parameters
    task_type: str
        task_type, such as "train"
    work_path : str
        work_path, such as "iter.000001/00.train"
    Returns
    -------
    None
    '''
    user_forward_files = mdata.get(task_type + "_" + "user_forward_files", [])
    #Angus: In the future, we may unify the task format.
    if task_format is None:
        task_format = {"train" : "0*", "model_devi" : "task.*", "fp": "task.*"}
        #"init_relax" : "sys-*", "init_md" : "sys-*/scale*/00*"
    for file in user_forward_files:
        assert os.path.isfile(file)  ,\
            "user_forward_file %s of %s stage doesn't exist. " % (file, task_type)
        tasks = glob.glob(os.path.join(work_path, task_format[task_type]))
        for task in tasks:
            if os.path.isfile(os.path.join(task, os.path.basename(file))):
                os.remove(os.path.join(task, os.path.basename(file)))
            os.symlink(file, os.path.join(task, os.path.basename(file)))
    return
    
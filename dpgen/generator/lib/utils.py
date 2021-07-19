#!/usr/bin/env python3

import os, re, shutil, logging
import json

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

def _get_param_alias(jdata,
                     names) :
    for ii in names :
        if ii in jdata :
            return jdata[ii]
    raise ValueError("one of the keys %s should be in jdata %s" % (str(names), (json.dumps(jdata, indent=4))))

def expand_idx (in_list) :
    ret = []
    for ii in in_list :
        if type(ii) == int :
            ret.append(ii)
        elif type(ii) == str:
            step_str = ii.split(':')
            if len(step_str) > 1 :
                step = int(step_str[1])
            else :
                step = 1
            range_str = step_str[0].split('-')
            assert(len(range_str)) == 2
            ret += range(int(range_str[0]), int(range_str[1]), step)
    return ret

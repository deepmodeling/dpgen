#!/usr/bin/env python3

import os,sys,json,glob,argparse
import numpy as np
import subprocess as sp

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def collect_data(target_folder, param_file, output, verbose = True) :
    target_folder = os.path.abspath(target_folder)
    output = os.path.abspath(output)
    tool_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'template')    
    command_cvt_2_raw = os.path.join(tool_path, 'tools.vasp', 'convert2raw.py')
    command_cvt_2_raw += ' data.configs'
    command_shuffle_raw = os.path.join(tool_path, 'tools.raw', 'shuffle_raw.py')
    command_raw_2_set = os.path.join(tool_path, 'tools.raw', 'raw_to_set.sh')
    # goto input    
    cwd = os.getcwd()
    os.chdir(target_folder)
    jdata = json.load(open(param_file))
    sys = jdata['sys_configs']
    if verbose :
        max_str_len = max([len(str(ii)) for ii in sys])
        ptr_fmt = '%%%ds   %%6d' % (max_str_len+5)
    # collect systems from iter dirs
    coll_sys = [[] for ii in sys]    
    numb_sys = len(sys)
    iters = glob.glob('iter.[0-9]*[0-9]')
    iters.sort()    
    for ii in iters :
        iter_data = glob.glob(os.path.join(ii, '02.fp', 'data.[0-9]*[0-9]'))
        iter_data.sort()
        for jj in iter_data :
            sys_idx = int(os.path.basename(jj).split('.')[-1])
            coll_sys[sys_idx].append(jj)
    # create output dir
    os.makedirs(output, exist_ok = True)
    # loop over systems
    for idx,ii in enumerate(coll_sys) :
        if len(ii) == 0 :
            continue
        # link iter data dirs
        out_sys_path = os.path.join(output, 'system.%03d' % idx)
        os.makedirs(out_sys_path, exist_ok=True)
        cwd_ = os.getcwd()
        os.chdir(out_sys_path)
        for jj in ii :
            in_sys_path = os.path.join(target_folder, jj)
            in_iter = in_sys_path.split('/')[-3]
            in_base = in_sys_path.split('/')[-1]
            out_file = in_iter + '.' + in_base
            if os.path.exists(out_file) :
                os.remove(out_file)
            os.symlink(in_sys_path, out_file)
        # cat data.configs
        data_configs = glob.glob(os.path.join('iter.[0-9]*[0-9].data.[0-9]*[0-9]', 'orig', 'data.configs'))
        data_configs.sort()
        os.makedirs('orig', exist_ok = True)        
        with open(os.path.join('orig', 'data.configs'), 'w') as outfile:
            for fname in data_configs:
                with open(fname) as infile:
                    outfile.write(infile.read())        
        # convert to raw
        os.chdir('orig')
        sp.check_call(command_cvt_2_raw, shell = True)
        os.chdir('..')
        # shuffle raw
        sp.check_call(command_shuffle_raw + ' orig ' + ' . > /dev/null', shell = True)
        if os.path.exists('type.raw') :
            os.remove('type.raw')
        os.symlink(os.path.join('orig', 'type.raw'), 'type.raw')
        # raw to sets
        sp.check_call(command_raw_2_set + ' > /dev/null', shell = True)
        # print summary
        if verbose :
            ndata = file_len('box.raw')
            print(ptr_fmt % (str(sys[idx]), ndata))
        # ch dir
        os.chdir(cwd_)


def _main()   :
    parser = argparse.ArgumentParser(description='Collect data from DP-GEN iterations')
    parser.add_argument("JOB_DIR", type=str, 
                        help="the directory of the DP-GEN job")
    parser.add_argument("OUTPUT", type=str, 
                        help="the output directory of data")
    parser.add_argument('-p',"--parameter", type=str, default = 'param.json',
                        help="the json file provides DP-GEN paramters, should be located in JOB_DIR")
    parser.add_argument('-v',"--verbose", action = 'store_true',
                        help="print number of data in each system")
    args = parser.parse_args()
    
    collect_data(args.JOB_DIR, args.parameter, args.OUTPUT, args.verbose)        

if __name__ == '__main__':
    _main()

    

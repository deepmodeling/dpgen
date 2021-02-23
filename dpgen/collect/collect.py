#!/usr/bin/env python3

import os,sys,json,glob,argparse,dpdata
import numpy as np
from dpgen.generator.run import data_system_fmt

def collect_data(target_folder, param_file, output, 
                 verbose = True, 
                 shuffle = True, 
                 merge = True) :
    target_folder = os.path.abspath(target_folder)
    output = os.path.abspath(output)
    # goto input    
    cwd = os.getcwd()
    os.chdir(target_folder)
    jdata = json.load(open(param_file))
    sys_configs_prefix = jdata.get('sys_configs_prefix', '')
    sys_configs = jdata.get('sys_configs', [])
    if verbose :
        max_str_len = max([len(str(ii)) for ii in sys_configs])
        max_form_len = 16
        ptr_fmt = '%%%ds %%%ds natoms %%6d nframes %%6d' % (max_str_len+5, max_form_len)
    # init systems
    init_data = []
    init_data_prefix = jdata.get('init_data_prefix', '')
    init_data_sys = jdata.get('init_data_sys', [])
    for ii in init_data_sys:
        init_data.append(dpdata.LabeledSystem(os.path.join(init_data_prefix, ii), fmt='deepmd/npy'))        
    # collect systems from iter dirs    
    coll_data = {}
    numb_sys = len(sys_configs)
    model_devi_jobs = jdata.get('model_devi_jobs', {})
    numb_jobs = len(model_devi_jobs)
    iters = ['iter.%06d' % ii for ii in range(numb_jobs)]
    # loop over iters to collect data
    for ii in range(len(iters)) :
        iter_data = glob.glob(os.path.join(iters[ii], '02.fp', 'data.[0-9]*[0-9]'))
        iter_data.sort()
        for jj in iter_data :
            sys = dpdata.LabeledSystem(jj, fmt = 'deepmd/npy')
            if merge:
                sys_str = sys.formula
            else:
                sys_str = (os.path.basename(jj).split('.')[-1])            
            if sys_str in coll_data.keys():
                coll_data[sys_str].append(sys)
            else:
                coll_data[sys_str] = sys
    # print information
    if verbose:
        for ii in range(len(init_data)):
            print(ptr_fmt % (str(init_data_sys[ii]), 
                             init_data[ii].formula,
                             init_data[ii].get_natoms(),
                             init_data[ii].get_nframes() ))            
        keys = list(coll_data.keys())
        keys.sort()
        for ii in keys:
            if merge:
                sys_str = ii
            else :
                sys_str = str(sys_configs[int(ii)])
            print(ptr_fmt % (sys_str, 
                             coll_data[ii].formula,
                             coll_data[ii].get_natoms(),
                             coll_data[ii].get_nframes() ))
    # shuffle system data
    if shuffle:
        for kk in coll_data.keys():
            coll_data[kk].shuffle()
    # create output dir
    os.chdir(cwd)
    os.makedirs(output, exist_ok = True)
    # dump init data
    for idx,ii in enumerate(init_data):
        out_dir = 'init.' + (data_system_fmt % idx)
        ii.to('deepmd/npy', os.path.join(output, out_dir))
    # dump iter data
    for kk in coll_data.keys():
        out_dir = 'sys.%s' % kk
        nframes = coll_data[kk].get_nframes()
        coll_data[kk].to('deepmd/npy', os.path.join(output, out_dir), set_size = nframes)
        # coll_data[kk].to('deepmd/npy', os.path.join(output, out_dir))

def gen_collect(args):
    collect_data(args.JOB_DIR, args.parameter, args.OUTPUT, 
                 verbose = args.verbose,
                 shuffle = args.shuffle,
                 merge = args.merge)

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
    parser.add_argument('-m',"--merge", action = 'store_true',
                             help="merge the systems with the same chemical formula")
    parser.add_argument('-s',"--shuffle", action = 'store_true',
                             help="shuffle the data systems")
    args = parser.parse_args()
    gen_collect(args)

if __name__ == '__main__':
    _main()

    

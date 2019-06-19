#!/usr/bin/env python3

import os, glob, argparse, json, re
import lib.lammps as lammps
import lib.vasp as vasp

global_task_name = '01.eos'

def comput_lmp_eos(conf_dir, task_name) :
    conf_path = re.sub('confs', global_task_name, conf_dir)
    conf_path = os.path.abspath(conf_path)
    conf_path = os.path.join(conf_path, task_name)
    vol_paths = glob.glob(os.path.join(conf_path, 'vol-*'))
    vol_paths.sort()
    for ii in vol_paths :
        log_lammps = os.path.join(ii, 'log.lammps')
        natoms, epa, vpa = lammps.get_nev(log_lammps)
        print(vpa, epa)

def comput_vasp_eos(jdata, conf_dir) :
    kspacing = jdata['vasp_params']['kspacing']
    conf_path = re.sub('confs', global_task_name, conf_dir)
    conf_path = os.path.abspath(conf_path)
    conf_path = os.path.join(conf_path, 'vasp-k%.2f' % kspacing)
    vol_paths = glob.glob(os.path.join(conf_path, 'vol-*'))
    vol_paths.sort()
    for ii in vol_paths :
        outcar = os.path.join(ii, 'OUTCAR')
        natoms, epa, vpa = vasp.get_nev(outcar)
        print(vpa, epa)

def _main():
    parser = argparse.ArgumentParser(
        description="cmpt 01.eos")
    parser.add_argument('TASK', type=str,
                        choices = ['vasp', 'deepmd', 'meam'], 
                        help='the task of generation, vasp or lammps')
    parser.add_argument('PARAM', type=str,
                        help='json parameter file')
    parser.add_argument('CONF', type=str,
                        help='the path to conf')
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

    if args.TASK == 'vasp':
        comput_vasp_eos(jdata, args.CONF)               
    elif args.TASK == 'deepmd' :
        comput_lmp_eos(args.CONF, args.TASK)
    elif args.TASK == 'meam' :
        comput_lmp_eos(args.CONF, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)

if __name__ == '__main__' :
    _main()
    

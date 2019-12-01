#!/usr/bin/env python3

import os, glob, argparse, json, re
import dpgen.auto_test.lib.util as util
import dpgen.auto_test.lib.lammps as lammps
import dpgen.auto_test.lib.vasp as vasp

global_task_name = '01.eos'

def comput_lmp_eos(jdata,conf_dir, task_name) :
    conf_path = re.sub('confs', global_task_name, conf_dir)
    conf_path = os.path.abspath(conf_path)
    conf_path = os.path.join(conf_path, task_name)
    vol_paths = glob.glob(os.path.join(conf_path, 'vol-*'))
    vol_paths.sort()
    result = os.path.join(conf_path,'result')
    print('Vpa(A^3)\tEpA(eV)')
    with open(result,'w') as fp:
        fp.write('conf_dir:%s\n VpA(A^3)  EpA(eV)\n'% (conf_dir))
        for ii in vol_paths :
            log_lammps = os.path.join(ii, 'log.lammps')
            natoms, epa, vpa = lammps.get_nev(log_lammps)
            print(vpa, epa)
            fp.write('%7.3f  %8.4f\n' % (vpa,epa))
    fp.close()
    if 'upload_username' in jdata.keys() and task_name =='deepmd':
        upload_username=jdata['upload_username']
        util.insert_data('eos','deepmd',upload_username,result)

def comput_vasp_eos(jdata, conf_dir) :
    conf_path = re.sub('confs', global_task_name, conf_dir)
    conf_path = os.path.abspath(conf_path)

    if 'relax_incar' in jdata.keys():
        vasp_str='vasp-relax_incar'
    else:
        kspacing = jdata['vasp_params']['kspacing']
        vasp_str='vasp-k%.2f' % kspacing
    task_path = os.path.join(conf_path, vasp_str)
    vol_paths = glob.glob(os.path.join(task_path, 'vol-*'))
    vol_paths.sort()
    result = os.path.join(task_path,'result')
    print('Vpa(A^3)\tEpA(eV)')
    with open(result,'w') as fp:
        fp.write('conf_dir:%s\n VpA(A^3)  EpA(eV)\n'% (conf_dir))
        for ii in vol_paths :
            outcar = os.path.join(ii, 'OUTCAR')
            natoms, epa, vpa = vasp.get_nev(outcar)
            print(vpa, epa)
            fp.write('%7.3f  %8.4f\n' % (vpa,epa))
    fp.close()
    if 'upload_username' in jdata.keys():
        upload_username=jdata['upload_username']
        util.insert_data('eos','vasp',upload_username,result)

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
        comput_lmp_eos(jdata,args.CONF, args.TASK)
    elif args.TASK == 'meam' :
        comput_lmp_eos(jdata,args.CONF, args.TASK)
    else :
        raise RuntimeError("unknow task ", args.TASK)

if __name__ == '__main__' :
    _main()

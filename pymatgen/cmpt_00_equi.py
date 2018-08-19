#!/usr/bin/env python3

import os, glob, argparse, json, warnings, re
import numpy as np
import lib.lammps as lammps
import lib.vasp as vasp

global_equi_name = '00.equi'

def comput_lmp_nev(conf_dir) :
    conf_path = re.sub('confs', global_equi_name, conf_dir)
    conf_path = os.path.abspath(conf_path)
    conf_path = os.path.join(conf_path, 'lmp')
    log_lammps = os.path.join(conf_path, 'log.lammps')
    if os.path.isfile(log_lammps):
        natoms, epa, vpa = lammps.get_nev(log_lammps)
        return natoms, epa, vpa
    else :
        return None, None, None

def comput_vasp_nev(jdata, conf_dir) :
    kspacing = jdata['vasp_params']['kspacing']
    conf_path = re.sub('confs', global_equi_name, conf_dir)
    conf_path = os.path.abspath(conf_path)
    vasp_path = os.path.join(conf_path, 'vasp-k%.2f' % kspacing)
    outcar = os.path.join(vasp_path, 'OUTCAR')
    # tag_fin = os.path.join(vasp_path, 'tag_finished')
    if not os.path.isfile(outcar) :
        warnings.warn("cannot find OUTCAR in "+vasp_path+" skip")
    elif not vasp.check_finished(outcar):
        warnings.warn("incomplete job "+vasp_path+" use the last frame")
    if os.path.isfile(outcar):
        natoms, epa, vpa = vasp.get_nev(outcar)
        return natoms, epa, vpa
    else :
        return None, None, None

def _main():
    parser = argparse.ArgumentParser(
        description="cmpt 00.equi")
    parser.add_argument('PARAM', type=str,
                        help='the json param')
    parser.add_argument('CONF', type=str,
                        help='the dir of conf')
    args = parser.parse_args()
    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)

    ln, le, lv = comput_lmp_nev(args.CONF)
    vn, ve, vv = comput_vasp_nev(jdata, args.CONF)
    if le == None or ve == None or lv == None or vv == None:
        print("%s" % args.CONF)
    else :
        print("%s\t %8.4f %8.4f %8.4f  %7.3f %7.3f %7.3f" % (args.CONF, ve, le, (le-ve), vv, lv, (le-ve)))

if __name__ == '__main__' :
    _main()
    

#!/usr/bin/env python3 

import glob, os, shutil

def copy(kspacing, static = False) :
    test_dir = os.path.abspath(os.path.dirname(__file__))
    if static :
        vasp_name = 'vasp-static-k%.2f' % kspacing
    else :
        vasp_name = 'vasp-k%.2f' % kspacing    
    confs = glob.glob(os.path.join(test_dir, "confs"))
    params = glob.glob(os.path.join(test_dir, "param.json"))
    cmpts = glob.glob(os.path.join(test_dir, "cmpt_*.py"))
    gens = glob.glob(os.path.join(test_dir, "gen_*.py"))
    vasp = glob.glob(os.path.join(test_dir, "0*/*/*/" + vasp_name))
    cwd = os.getcwd()
    os.chdir(test_dir)
    vasp_dirs = glob.glob( "0*/*/*/" + vasp_name)
    os.chdir(cwd)
    for res,ii in zip(vasp,vasp_dirs) :
        dir_ii = os.path.dirname(ii)
        os.makedirs(dir_ii, exist_ok = True)
        os.chdir(dir_ii)
        if os.path.islink(vasp_name) :
            os.remove(vasp_name)
        os.symlink(os.path.relpath(res), vasp_name)
        os.chdir(cwd)

    os.makedirs('results', exist_ok = True)
    for ii in cmpts :
        fname_ii = os.path.basename(ii)
        if os.path.islink(fname_ii) :
            os.remove(fname_ii)
        os.symlink(os.path.relpath(ii), fname_ii)
    for ii in gens :
        fname_ii = os.path.basename(ii)
        if os.path.islink(fname_ii) :
            os.remove(fname_ii)
        os.symlink(os.path.relpath(ii), fname_ii)
    for ii in confs :
        fname_ii = os.path.basename(ii)
        if os.path.islink(fname_ii) :
            os.remove(fname_ii)
        os.symlink(os.path.relpath(ii), fname_ii)
    for ii in params :
        fname_ii = os.path.basename(ii)
        if os.path.isfile(fname_ii) :
            os.remove(fname_ii)
        shutil.copy2(os.path.relpath(ii), fname_ii)

copy(0.08)
copy(0.08, static = True)
copy(0.16)

#!/usr/bin/env python3 

import os
import re
import sys
import argparse
import glob
import json
import random
import logging
import warnings
import shutil
import time
import dpdata
import numpy as np
from dpgen import dlog
import subprocess as sp
import dpgen.data.tools.hcp as hcp
import dpgen.data.tools.fcc as fcc
import dpgen.data.tools.bcc as bcc
import dpgen.data.tools.diamond as diamond
import dpgen.data.tools.sc as sc
from distutils.version import LooseVersion
from dpgen.generator.lib.vasp import incar_upper
from pymatgen.core import Structure
from pymatgen.io.vasp import Incar
from dpgen.remote.decide_machine import  decide_fp_machine
from dpgen import ROOT_PATH
from dpgen.dispatcher.Dispatcher import Dispatcher, make_dispatcher, make_submission



def create_path (path,back=False) :
    if  path[-1] != "/":
        path += '/'
    if os.path.isdir(path) : 
        if back:
           dirname = os.path.dirname(path)        
           counter = 0
           while True :
               bk_dirname = dirname + ".bk%03d" % counter
               if not os.path.isdir(bk_dirname) : 
                   shutil.move (dirname, bk_dirname) 
                   break
               counter += 1
           os.makedirs (path)
           return path
        else:
           return path

    os.makedirs (path)
    return path

def replace (file_name, pattern, subst) :
    file_handel = open (file_name, 'r')
    file_string = file_handel.read ()
    file_handel.close ()
    file_string = ( re.sub (pattern, subst, file_string) )
    file_handel = open (file_name, 'w')
    file_handel.write (file_string)
    file_handel.close ()

"""
0, make unit cell
1, copy
2, place element
3, relax
4, perturb
"""
global_dirname_02 = '00.place_ele'
global_dirname_03 = '01.scale_pert'
global_dirname_04 = '02.md'

def out_dir_name(jdata) :
    cell_type = jdata['cell_type']
    elements = jdata['elements']
    super_cell = jdata['super_cell']    
    from_poscar = False
    if 'from_poscar' in jdata :
        from_poscar = jdata['from_poscar']
        from_poscar_path = jdata['from_poscar_path']

    if from_poscar:
        poscar_name = os.path.basename(from_poscar_path)
        cell_str = "%02d" % (super_cell[0])
        for ii in range(1,len(super_cell)) :
            cell_str = cell_str + ("x%02d" % super_cell[ii])
        return poscar_name + '.' + cell_str
    else :
        ele_str = ""
        for ii in elements:
            ele_str = ele_str + ii.lower()
        cell_str = "%02d" % (super_cell[0])
        for ii in range(1,len(super_cell)) :
            cell_str = cell_str + ("x%02d" % super_cell[ii])
        return ele_str + '.' + cell_type + '.' + cell_str

def class_cell_type(jdata) :
    ct = jdata['cell_type']
    if ct == "hcp" :
        cell_type = hcp
    elif ct == "fcc" :
        cell_type = fcc
    elif ct == "diamond" :
        cell_type = diamond
    elif ct == "sc" :
        cell_type = sc
    elif ct == "bcc" :
        cell_type = bcc
    else :
        raise RuntimeError("unknown cell type %s" % ct)
    return cell_type

def poscar_ele(poscar_in, poscar_out, eles, natoms) :
    ele_line = ""
    natom_line = ""
    for ii in eles :
        ele_line += str(ii) + " "
    for ii in natoms :
        natom_line += str(ii) + " "
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
        lines[5] = ele_line + "\n"
        lines[6] = natom_line + "\n"
    with open(poscar_out, 'w') as fout :
        fout.write("".join(lines))

def poscar_natoms(lines) :
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    return numb_atoms

def poscar_shuffle(poscar_in, poscar_out) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    numb_atoms = poscar_natoms(lines)
    idx = np.arange(8, 8+numb_atoms)
    np.random.shuffle(idx)
    out_lines = lines[0:8]
    for ii in range(numb_atoms) :
        out_lines.append(lines[idx[ii]])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(out_lines))

def poscar_scale_direct (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = poscar_natoms(lines)
    pscale = float(lines[1])
    pscale = pscale * scale
    lines[1] = str(pscale) + "\n"
    return lines

def poscar_scale_cartesian (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = poscar_natoms(lines)
    # scale box
    for ii in range(2,5) :
        boxl = lines[ii].split()
        boxv = [float(ii) for ii in boxl]
        boxv = np.array(boxv) * scale
        lines[ii] = "%.16e %.16e %.16e\n" % (boxv[0], boxv[1], boxv[2])
    # scale coord
    for ii in range(8, 8+numb_atoms) :
        cl = lines[ii].split()
        cv = [float(ii) for ii in cl]
        cv = np.array(cv) * scale
        lines[ii] = "%.16e %.16e %.16e\n" % (cv[0], cv[1], cv[2])
    return lines    

def poscar_scale (poscar_in, poscar_out, scale) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    if 'D' == lines[7][0] or 'd' == lines[7][0]: 
        lines = poscar_scale_direct(lines, scale)
    elif 'C' == lines[7][0] or 'c' == lines[7][0] :
        lines = poscar_scale_cartesian(lines, scale)
    else :
        raise RuntimeError("Unknow poscar style at line 7: %s" % lines[7])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(lines))

def make_unit_cell (jdata) :
    latt = jdata['latt']
    out_dir = jdata['out_dir']
    path_uc = os.path.join(out_dir, global_dirname_02)
    cell_type = class_cell_type(jdata)

    cwd = os.getcwd()    
    # for ii in scale :
    # path_work = create_path(os.path.join(path_uc, '%.3f' % ii))
    path_work = create_path(path_uc)    
    os.chdir(path_work)
    with open('POSCAR.unit', 'w') as fp:
        fp.write (cell_type.poscar_unit(latt))
    os.chdir(cwd)        

def make_super_cell (jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    path_uc = os.path.join(out_dir, global_dirname_02)
    path_sc = os.path.join(out_dir, global_dirname_02)
    assert(os.path.isdir(path_uc)), "path %s should exists" % path_uc
    assert(os.path.isdir(path_sc)), "path %s should exists" % path_sc

    # for ii in scale :
    from_path = path_uc
    from_file = os.path.join(from_path, 'POSCAR.unit')
    to_path = path_sc
    to_file = os.path.join(to_path, 'POSCAR')

    #minor bug for element symbol behind the coordinates
    from_struct=Structure.from_file(from_file)
    from_struct.make_supercell(super_cell)
    from_struct.to('poscar',to_file)

def make_super_cell_poscar(jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    path_sc = os.path.join(out_dir, global_dirname_02)    
    create_path(path_sc)
    from_poscar_path = jdata['from_poscar_path']
    assert(os.path.isfile(from_poscar_path)), "file %s should exists" % from_poscar_path
    
    from_file = os.path.join(path_sc, 'POSCAR.copied')
    shutil.copy2(from_poscar_path, from_file)
    to_path = path_sc
    to_file = os.path.join(to_path, 'POSCAR')
  
    #minor bug for element symbol behind the coordinates
    from_struct=Structure.from_file(from_file)
    from_struct.make_supercell(super_cell)
    from_struct.to('poscar',to_file)  

    # make system dir (copy)
    lines = open(to_file, 'r').read().split('\n')
    natoms_str = lines[6]
    natoms_list = [int(ii) for ii in natoms_str.split()]
    dlog.info(natoms_list)
    comb_name = "sys-"
    for idx,ii in enumerate(natoms_list) :
        comb_name += "%04d" % ii
        if idx != len(natoms_list)-1 :
            comb_name += "-"
    path_work = os.path.join(path_sc, comb_name)
    create_path(path_work)
    cwd = os.getcwd()
    to_file = os.path.abspath(to_file)
    os.chdir(path_work)
    try:
        os.symlink(os.path.relpath(to_file), 'POSCAR')
    except FileExistsError:
        pass
    os.chdir(cwd)

def make_combines (dim, natoms) :
    if dim == 1 :
        return [[natoms]]
    else :
        res = []
        for ii in range(natoms+1) :
            rest = natoms - ii
            tmp_combines = make_combines(dim-1, rest)
            for jj in tmp_combines :
                jj.append(ii)
            if len(res) == 0 :
                res = tmp_combines
            else :
                res += tmp_combines
        return res

def place_element (jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    cell_type = class_cell_type(jdata)
    natoms = np.cumprod(super_cell)[-1] * cell_type.numb_atoms()
    elements = jdata['elements']
    path_sc = os.path.join(out_dir, global_dirname_02)
    path_pe = os.path.join(out_dir, global_dirname_02)    
    combines = np.array(make_combines(len(elements), natoms), dtype = int)
    
    assert(os.path.isdir(path_pe))
    cwd = os.getcwd()
    for ii in combines :
        if any(ii == 0) :
            continue
        comb_name = "sys-"
        for idx,jj in enumerate(ii) :            
            comb_name += "%04d" % jj
            if idx != len(ii)-1 :
                comb_name += "-"
        path_pos_in = path_sc
        path_work = os.path.join(path_pe, comb_name)
        create_path(path_work)
        pos_in = os.path.join(path_pos_in, 'POSCAR')
        pos_out = os.path.join(path_work, 'POSCAR')
        poscar_ele(pos_in, pos_out, elements, ii)
        poscar_shuffle(pos_out, pos_out)

def make_vasp_relax (jdata, mdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    cwd = os.getcwd()
    work_dir = os.path.join(out_dir, global_dirname_02)
    assert (os.path.isdir(work_dir))
    work_dir = os.path.abspath(work_dir)

    if os.path.isfile(os.path.join(work_dir, 'INCAR' )) :
        os.remove(os.path.join(work_dir, 'INCAR' ))
    if os.path.isfile(os.path.join(work_dir, 'POTCAR')) :
        os.remove(os.path.join(work_dir, 'POTCAR'))
    shutil.copy2( jdata['relax_incar'], 
                 os.path.join(work_dir, 'INCAR'))
    is_cvasp = False
    if 'cvasp' in mdata['fp_resources'].keys():
        is_cvasp = mdata['fp_resources']['cvasp']
    if is_cvasp:
        cvasp_file=os.path.join(ROOT_PATH,'generator/lib/cvasp.py')
        shutil.copyfile(cvasp_file, os.path.join(work_dir, 'cvasp.py'))
    out_potcar = os.path.join(work_dir, 'POTCAR')
    with open(out_potcar, 'w') as outfile:
        for fname in potcars:
            with open(fname) as infile:
                outfile.write(infile.read())
    
    os.chdir(work_dir)
    
    sys_list = glob.glob('sys-*')
    for ss in sys_list:
        os.chdir(ss)
        ln_src = os.path.relpath(os.path.join(work_dir,'INCAR'))
        try:
           os.symlink(ln_src, 'INCAR')
        except FileExistsError:
           pass
        ln_src = os.path.relpath(os.path.join(work_dir,'POTCAR'))
        try:
           os.symlink(ln_src, 'POTCAR')
        except FileExistsError:
           pass
        os.chdir(work_dir)
    os.chdir(cwd)

def make_scale(jdata):
    out_dir = jdata['out_dir']
    scale = jdata['scale']    
    skip_relax = jdata['skip_relax']    

    cwd = os.getcwd()
    init_path = os.path.join(out_dir, global_dirname_02)
    init_path = os.path.abspath(init_path)
    work_path = os.path.join(out_dir, global_dirname_03)
    os.chdir(init_path)
    init_sys = glob.glob("sys-*")
    init_sys.sort()
    os.chdir(cwd)

    create_path(work_path)
    for ii in init_sys :
        for jj in scale :
            if skip_relax :
                pos_src = os.path.join(os.path.join(init_path, ii), 'POSCAR')
                assert(os.path.isfile(pos_src))
            else :
                try:
                    pos_src = os.path.join(os.path.join(init_path, ii), 'CONTCAR')
                    assert(os.path.isfile(pos_src))
                except:
                    raise RuntimeError("not file %s, vasp relaxation should be run before scale poscar")
            scale_path = os.path.join(work_path, ii)
            scale_path = os.path.join(scale_path, "scale-%.3f" % jj)
            create_path(scale_path)
            os.chdir(scale_path) 
            poscar_scale(pos_src, 'POSCAR', jj)
            os.chdir(cwd)

def pert_scaled(jdata) :
    out_dir = jdata['out_dir']
    scale = jdata['scale']    
    pert_box = jdata['pert_box']
    pert_atom = jdata['pert_atom']
    pert_numb = jdata['pert_numb']
    from_poscar = False 
    if 'from_poscar' in jdata :
        from_poscar = jdata['from_poscar']
    
    cwd = os.getcwd()
    path_sp = os.path.join(out_dir, global_dirname_03)
    assert(os.path.isdir(path_sp))
    os.chdir(path_sp)
    sys_pe = glob.glob('sys-*')
    sys_pe.sort()
    os.chdir(cwd)    

    pert_cmd = os.path.dirname(__file__)
    pert_cmd = os.path.join(pert_cmd, 'tools')
    pert_cmd = os.path.join(pert_cmd, 'create_random_disturb.py')
    pert_cmd = 'python3 ' + pert_cmd + ' -etmax %f -ofmt vasp POSCAR %d %f > /dev/null' %(pert_box, pert_numb, pert_atom)
    for ii in sys_pe :
        for jj in scale :
            path_work = path_sp
            path_work = os.path.join(path_work, ii)
            path_work = os.path.join(path_work, 'scale-%.3f' % jj)
            assert(os.path.isdir(path_work))
            os.chdir(path_work)
            sp.check_call(pert_cmd, shell = True)
            for kk in range(pert_numb) :
                pos_in = 'POSCAR%d.vasp' % (kk+1)
                dir_out = '%06d' % (kk+1)
                create_path(dir_out)
                pos_out = os.path.join(dir_out, 'POSCAR')
                if not from_poscar:
                    poscar_shuffle(pos_in, pos_out)
                else :
                    shutil.copy2(pos_in, pos_out)
                os.remove(pos_in)
            kk = -1
            pos_in = 'POSCAR'
            dir_out = '%06d' % (kk+1)
            create_path(dir_out)
            pos_out = os.path.join(dir_out, 'POSCAR')
            if not from_poscar:
                poscar_shuffle(pos_in, pos_out)
            else :
                shutil.copy2(pos_in, pos_out)
            os.chdir(cwd)

def make_vasp_md(jdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    scale = jdata['scale']   
    pert_numb = jdata['pert_numb'] 
    md_nstep = jdata['md_nstep']

    cwd = os.getcwd()
    path_ps = os.path.join(out_dir, global_dirname_03)
    path_ps = os.path.abspath(path_ps)
    assert(os.path.isdir(path_ps))
    os.chdir(path_ps)
    sys_ps = glob.glob('sys-*')
    sys_ps.sort()
    os.chdir(cwd) 
    path_md = os.path.join(out_dir, global_dirname_04)
    path_md = os.path.abspath(path_md)
    create_path(path_md)
    shutil.copy2(jdata['md_incar'], 
                 os.path.join(path_md, 'INCAR'))
    out_potcar = os.path.join(path_md, 'POTCAR')
    with open(out_potcar, 'w') as outfile:
        for fname in potcars:
            with open(fname) as infile:
                outfile.write(infile.read())
    os.chdir(path_md)
    os.chdir(cwd)    

    for ii in sys_ps :
        for jj in scale :
            for kk in range(pert_numb) :
                path_work = path_md
                path_work = os.path.join(path_work, ii)
                path_work = os.path.join(path_work, "scale-%.3f" % jj)
                path_work = os.path.join(path_work, "%06d" % kk)
                create_path(path_work)
                os.chdir(path_work)                
                path_pos = path_ps
                path_pos = os.path.join(path_pos, ii)
                path_pos = os.path.join(path_pos, "scale-%.3f" % jj)
                path_pos = os.path.join(path_pos, "%06d" % kk)
                init_pos = os.path.join(path_pos, 'POSCAR')
                shutil.copy2 (init_pos, 'POSCAR')
                file_incar = os.path.join(path_md, 'INCAR')
                file_potcar = os.path.join(path_md, 'POTCAR')
                try:
                    os.symlink(os.path.relpath(file_incar), 'INCAR')
                except FileExistsError:
                    pass
                try:
                    os.symlink(os.path.relpath(file_potcar), 'POTCAR')
                except FileExistsError:
                    pass
                 
                os.chdir(cwd)                

def coll_vasp_md(jdata) :
    out_dir = jdata['out_dir']
    md_nstep = jdata['md_nstep']
    scale = jdata['scale']    
    pert_numb = jdata['pert_numb']
    coll_ndata = jdata['coll_ndata']

    cwd = os.getcwd()
    path_md = os.path.join(out_dir, global_dirname_04)
    path_md = os.path.abspath(path_md)
    assert(os.path.isdir(path_md)), "md path should exists"
    os.chdir(path_md)
    sys_md = glob.glob('sys-*')
    sys_md.sort()

    for ii in sys_md :
        os.chdir(ii)
        # convert outcars
        valid_outcars = []
        for jj in scale :
            for kk in range(pert_numb) :
                path_work = os.path.join("scale-%.3f" % jj, "%06d" % kk)
                outcar = os.path.join(path_work, 'OUTCAR')
                #dlog.info("OUTCAR",outcar)
                if os.path.isfile(outcar) :
                    #dlog.info("*"*40)
                    with open(outcar, 'r') as fin:
                        nforce = fin.read().count('TOTAL-FORCE')
                    #dlog.info("nforce is", nforce)
                    #dlog.info("md_nstep", md_nstep)
                    if nforce == md_nstep :
                        valid_outcars.append(outcar)
                    else:
                        dlog.info("WARNING : in directory %s nforce in OUTCAR is not equal to settings in INCAR"%(os.getcwd()))
        arg_cvt = " "
        if len(valid_outcars) == 0:
            raise RuntimeError("MD dir: %s: find no valid outcar in sys %s, "
                               "check if your vasp md simulation is correctly done" 
                               % (path_md, ii)) 

        flag=True
        if ("type_map" in jdata) and isinstance(jdata["type_map"], list):
            type_map = jdata["type_map"]
        else:
            type_map = None 
        for oo in valid_outcars :
            if flag:
                _sys = dpdata.LabeledSystem(oo, type_map= type_map)
                if len(_sys)>0:
                   all_sys=_sys
                   flag=False
                else:
                   pass
            else:
                _sys = dpdata.LabeledSystem(oo, type_map= type_map)
                if len(_sys)>0:
                   all_sys.append(_sys)
        # create deepmd data
        if all_sys.get_nframes() >= coll_ndata :
            all_sys = all_sys.sub_system(np.arange(coll_ndata))
        all_sys.to_deepmd_raw('deepmd')
        all_sys.to_deepmd_npy('deepmd', set_size = all_sys.get_nframes())
        os.chdir(path_md)
    os.chdir(cwd)

def _vasp_check_fin (ii) :
    if os.path.isfile(os.path.join(ii, 'OUTCAR')) :
        with open(os.path.join(ii, 'OUTCAR'), 'r') as fp :
            content = fp.read()
            count = content.count('Elapse')
            if count != 1 :
                return False
    else :
        return False
    return True

def run_vasp_relax(jdata, mdata):
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']
    #machine_type = mdata['fp_machine']['machine_type']
    work_dir = os.path.join(jdata['out_dir'], global_dirname_02)
    
    forward_files = ["POSCAR", "INCAR", "POTCAR"]
    backward_files = ["OUTCAR","CONTCAR"]
    forward_common_files = []
    if 'cvasp' in mdata['fp_resources']:
        if mdata['fp_resources']['cvasp']:
            forward_common_files=['cvasp.py']
    relax_tasks = glob.glob(os.path.join(work_dir, "sys-*"))
    relax_tasks.sort()
    #dlog.info("work_dir",work_dir)
    #dlog.info("relax_tasks",relax_tasks)
    if len(relax_tasks) == 0:
        return

    relax_run_tasks = relax_tasks
    #for ii in relax_tasks : 
    #    if not _vasp_check_fin(ii):
    #        relax_run_tasks.append(ii)
    run_tasks = [os.path.basename(ii) for ii in relax_run_tasks]

    api_version = mdata.get('api_version', '0.9')
    if LooseVersion(api_version) < LooseVersion('1.0'):
        warnings.warn(f"the dpdispatcher will be updated to new version."
            f"And the interface may be changed. Please check the documents for more details")
        dispatcher = make_dispatcher(mdata['fp_machine'], mdata['fp_resources'], work_dir, run_tasks, fp_group_size)
        dispatcher.run_jobs(fp_resources,
                       [fp_command],
                       work_dir,
                       run_tasks,
                       fp_group_size,
                       forward_common_files,
                       forward_files,
                       backward_files)

    elif LooseVersion(api_version) >= LooseVersion('1.0'):
        submission = make_submission(
            mdata['fp_machine'],
            mdata['fp_resources'],
            commands=[fp_command],
            work_path=work_dir,
            run_tasks=run_tasks,
            group_size=fp_group_size,
            forward_common_files=forward_common_files,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog = 'fp.log',
            errlog = 'fp.log')
        submission.run_submission()


def run_vasp_md(jdata, mdata):
    fp_command = mdata['fp_command']
    fp_group_size = mdata['fp_group_size']
    fp_resources = mdata['fp_resources']
    #machine_type = mdata['fp_machine']['machine_type']
    work_dir = os.path.join(jdata['out_dir'], global_dirname_04)
    scale = jdata['scale']   
    pert_numb = jdata['pert_numb'] 
    md_nstep = jdata['md_nstep']

    forward_files = ["POSCAR", "INCAR", "POTCAR"]
    backward_files = ["OUTCAR"]
    forward_common_files = []
    if 'cvasp' in mdata['fp_resources']:
        if mdata['fp_resources']['cvasp']:
            forward_common_files=['cvasp.py']

    path_md = work_dir
    path_md = os.path.abspath(path_md)
    cwd = os.getcwd()
    assert(os.path.isdir(path_md)), "md path should exists"
    md_tasks = glob.glob(os.path.join(work_dir, 'sys-*/scale*/00*'))
    md_tasks.sort()

    if len(md_tasks) == 0:
        return

    md_run_tasks = md_tasks
    #for ii in md_tasks : 
    #    if not _vasp_check_fin(ii):
    #        md_run_tasks.append(ii)

    run_tasks = [ii.replace(work_dir+"/", "") for ii in md_run_tasks]
    #dlog.info("md_work_dir", work_dir)
    #dlog.info("run_tasks",run_tasks)
    api_version = mdata.get('api_version', '0.9')
    if LooseVersion(api_version) < LooseVersion('1.0'):
        warnings.warn(f"the dpdispatcher will be updated to new version."
            f"And the interface may be changed. Please check the documents for more details")
        dispatcher = make_dispatcher(mdata['fp_machine'], mdata['fp_resources'], work_dir, run_tasks, fp_group_size)
        dispatcher.run_jobs(fp_resources,
                       [fp_command],
                       work_dir,
                       run_tasks,
                       fp_group_size,
                       forward_common_files,
                       forward_files,
                       backward_files)

    elif LooseVersion(api_version) >= LooseVersion('1.0'):
        submission = make_submission(
            mdata['fp_machine'],
            mdata['fp_resources'],
            commands=[fp_command],
            work_path=work_dir,
            run_tasks=run_tasks,
            group_size=fp_group_size,
            forward_common_files=forward_common_files,
            forward_files=forward_files,
            backward_files=backward_files,
            outlog = 'fp.log',
            errlog = 'fp.log')
        submission.run_submission()

def gen_init_bulk(args) :
    try:
       import ruamel
       from monty.serialization import loadfn,dumpfn
       warnings.simplefilter('ignore', ruamel.yaml.error.MantissaNoDotYAML1_1Warning)
       jdata=loadfn(args.PARAM)
       if args.MACHINE is not None:
          mdata=loadfn(args.MACHINE)
    except:
       with open (args.PARAM, 'r') as fp :
           jdata = json.load (fp)
       if args.MACHINE is not None:
          with open (args.MACHINE, "r") as fp:
              mdata = json.load(fp)

    if args.MACHINE is not None:
       # Selecting a proper machine
       mdata = decide_fp_machine(mdata)
       #disp = make_dispatcher(mdata["fp_machine"])

    # Decide work path
    out_dir = out_dir_name(jdata)
    jdata['out_dir'] = out_dir
    dlog.info ("# working dir %s" % out_dir)
    # Decide whether to use a given poscar
    from_poscar = False 
    if 'from_poscar' in jdata :
        from_poscar = jdata['from_poscar']
    # Verify md_nstep
    md_nstep_jdata = jdata["md_nstep"]
    try:
        md_incar = jdata['md_incar']
        if os.path.isfile(md_incar):
            standard_incar = incar_upper(Incar.from_file(md_incar))
            nsw_flag = False
            if "NSW" in standard_incar:
                    nsw_flag = True
                    nsw_steps = standard_incar['NSW']
            #dlog.info("nsw_steps is", nsw_steps)
            #dlog.info("md_nstep_jdata is", md_nstep_jdata)
            if nsw_flag:
                if (nsw_steps != md_nstep_jdata):
                    dlog.info("WARNING: your set-up for MD steps in PARAM and md_incar are not consistent!")
                    dlog.info("MD steps in PARAM is %d"%(md_nstep_jdata))
                    dlog.info("MD steps in md_incar is %d"%(nsw_steps))
                    dlog.info("DP-GEN will use settings in md_incar!")
                    jdata['md_nstep'] = nsw_steps
    except KeyError:
        pass
    ## correct element name 
    temp_elements = []
    for ele in jdata['elements']:
        temp_elements.append(ele[0].upper() + ele[1:])
    jdata['elements'] = temp_elements
    dlog.info("Elements are %s"% ' '.join(jdata['elements']))

    ## Iteration 
    stage_list = [int(i) for i in jdata['stages']]
    for stage in stage_list:
        if stage == 1 :
            dlog.info("Current stage is 1, relax")
            create_path(out_dir)
            shutil.copy2(args.PARAM, os.path.join(out_dir, 'param.json'))
            if from_poscar :
                make_super_cell_poscar(jdata)
            else :
                make_unit_cell(jdata)
                make_super_cell(jdata)
                place_element(jdata)
            if args.MACHINE is not None:
               make_vasp_relax(jdata, mdata)
               run_vasp_relax(jdata, mdata)
            else:
               make_vasp_relax(jdata, {"fp_resources":{}})
        elif stage == 2 :
            dlog.info("Current stage is 2, perturb and scale")
            make_scale(jdata)
            pert_scaled(jdata)
        elif stage == 3 :
            dlog.info("Current stage is 3, run a short md")
            make_vasp_md(jdata)
            if args.MACHINE is not None:
               run_vasp_md(jdata, mdata)
        elif stage == 4 :
            dlog.info("Current stage is 4, collect data")
            coll_vasp_md(jdata)
        else :
            raise RuntimeError("unknown stage %d" % stage)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generating initial data for bulk systems.")
    parser.add_argument('PARAM', type=str,
                        help="parameter file, json/yaml format")
    parser.add_argument('MACHINE', type=str,default=None,nargs="?",
                        help="machine file, json/yaml format")
    args = parser.parse_args()
    gen_init_bulk(args)

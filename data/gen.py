#!/usr/bin/env python3 

import os,json,shutil
import numpy as np
import subprocess as sp
import tools.hcp as hcp
import tools.fcc as fcc

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
    return path

"""
0, make unit cell
1, copy
2, place element
3, relax
4, perturb
"""
global_dirname_00 = '00.unit_cell'
global_dirname_01 = '01.super_cell'
global_dirname_02 = '02.place_ele'

def class_cell_type(jdata) :
    ct = jdata['cell_type']
    if ct == "hcp" :
        cell_type = hcp
    elif ct == "fcc" :
        cell_type = fcc
    else :
        raise RuntimeError("unknow cell type %s" % ct)
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

def poscar_shuffle(poscar_in, poscar_out) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    idx = np.arange(8, 8+numb_atoms)
    np.random.shuffle(idx)
    out_lines = lines[0:8]
    for ii in range(numb_atoms) :
        out_lines.append(lines[idx[ii]])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(out_lines))

def make_unit_cell (jdata) :
    latt = jdata['latt']
    scale = jdata['scale']
    out_dir = jdata['out_dir']
    path_uc = os.path.join(out_dir, global_dirname_00)
    cell_type = class_cell_type(jdata)

    cwd = os.getcwd()    
    create_path(path_uc)
    for ii in scale :
        path_work = create_path(os.path.join(path_uc, '%.3f' % ii))
        os.chdir(path_work)
        with open('POSCAR', 'w') as fp:
            fp.write (cell_type.poscar_unit(latt * ii))
        os.chdir(cwd)        

def make_super_cell (jdata) :
    scale = jdata['scale']
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    path_uc = os.path.join(out_dir, global_dirname_00)
    path_sc = os.path.join(out_dir, global_dirname_01)

    create_path(path_sc)
    for ii in scale :
        from_path = os.path.join(path_uc, '%.3f' % ii)
        from_file = os.path.join(from_path, 'POSCAR')
        to_path = create_path(os.path.join(path_sc, '%.3f' % ii))
        to_file = os.path.join(to_path, 'POSCAR')
        cmd = "./tools/copy.py -n %d %d %d " % (super_cell[0], super_cell[1], super_cell[2]) + \
              from_file + " " + \
              to_file
        sp.check_call(cmd, shell = True)


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
    scale = jdata['scale']
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    cell_type = class_cell_type(jdata)
    natoms = np.cumprod(super_cell)[-1] * cell_type.numb_atoms()
    elements = jdata['elements']
    path_sc = os.path.join(out_dir, global_dirname_01)
    path_pe = os.path.join(out_dir, global_dirname_02)    
    combines = np.array(make_combines(len(elements), natoms), dtype = int)
    
    create_path(path_pe)
    cwd = os.getcwd()
    for ii in combines :
        if any(ii == 0) :
            continue
        comb_name = ""
        for idx,jj in enumerate(ii) :            
            comb_name += "%04d" % jj
            if idx != len(ii)-1 :
                comb_name += "-"
        for jj in scale :
            path_pos_in = os.path.join(path_sc, '%.3f' % jj)
            path_work = os.path.join(path_pe, comb_name)
            path_work = os.path.join(path_work, '%.3f' % jj)
            create_path(path_work)
            pos_in = os.path.join(path_pos_in, 'POSCAR')
            pos_out = os.path.join(path_work, 'POSCAR')
            poscar_ele(pos_in, pos_out, elements, ii)
            poscar_shuffle(pos_out, pos_out)

with open ('param.json', 'r') as fp :
    jdata = json.load (fp)

# make_unit_cell(jdata)
# make_super_cell(jdata)
place_element(jdata)
# poscar_shuffle('POSCAR', 'POSCAR.out')

# sp.check_call("./tools/copy.py -n 2 2 2 POSCAR POSCAR.out", shell = True)
    

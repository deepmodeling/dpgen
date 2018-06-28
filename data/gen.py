#!/usr/bin/env python3 

import os,json,shutil,re,glob,argparse
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
global_dirname_00 = '00.unit_cell'
global_dirname_01 = '01.super_cell'
global_dirname_02 = '02.place_ele'
global_dirname_03 = '03.scale_pert'
global_dirname_04 = '04.md'

def out_dir_name(jdata) :
    cell_type = jdata['cell_type']
    elements = jdata['elements']
    super_cell = jdata['super_cell']    

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
    if 'Direct' in lines[7] : 
        lines = poscar_scale_direct(lines, scale)
    elif 'Cartesian' in lines[7] :
        lines = poscar_scale_cartesian(lines, scale)
    else :
        raise RuntimeError("Unknow poscar style at line 7: %s" % lines[7])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(lines))

def make_unit_cell (jdata) :
    latt = jdata['latt']
    out_dir = jdata['out_dir']
    path_uc = os.path.join(out_dir, global_dirname_00)
    cell_type = class_cell_type(jdata)

    cwd = os.getcwd()    
    # for ii in scale :
    # path_work = create_path(os.path.join(path_uc, '%.3f' % ii))
    path_work = create_path(path_uc)    
    os.chdir(path_work)
    with open('POSCAR', 'w') as fp:
        fp.write (cell_type.poscar_unit(latt))
    os.chdir(cwd)        

def make_super_cell (jdata) :
    out_dir = jdata['out_dir']
    super_cell = jdata['super_cell']
    path_uc = os.path.join(out_dir, global_dirname_00)
    path_sc = os.path.join(out_dir, global_dirname_01)

    # for ii in scale :
    from_path = path_uc
    from_file = os.path.join(from_path, 'POSCAR')
    to_path = create_path(path_sc)
    to_file = os.path.join(to_path, 'POSCAR')
    cmd = "./tools/poscar_copy.py -n %d %d %d " % (super_cell[0], super_cell[1], super_cell[2]) + \
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

def make_vasp_relax (jdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    encut = jdata['encut']
    kspacing = jdata['kspacing']
    kgamma = jdata['kgamma']
    cwd = os.getcwd()
    vasp_dir = os.path.join(cwd, 'vasp.in')

    work_dir = os.path.join(out_dir, global_dirname_02)
    assert (os.path.isdir(work_dir))
    work_dir = os.path.abspath(work_dir)
    if os.path.isfile(os.path.join(work_dir, 'INCAR' )) :
        os.remove(os.path.join(work_dir, 'INCAR' ))
    if os.path.isfile(os.path.join(work_dir, 'POTCAR')) :
        os.remove(os.path.join(work_dir, 'POTCAR'))
    shutil.copy2(os.path.join(vasp_dir, 'INCAR.rlx' ), 
                 os.path.join(work_dir, 'INCAR'))
    out_potcar = os.path.join(work_dir, 'POTCAR')
    with open(out_potcar, 'w') as outfile:
        for fname in potcars:
            with open(fname) as infile:
                outfile.write(infile.read())
    
    os.chdir(work_dir)
    replace('INCAR', 'ENCUT=.*', 'ENCUT=%f' % encut)
    replace('INCAR', 'ISIF=.*', 'ISIF=3')
    replace('INCAR', 'KSPACING=.*', 'KSPACING=%f' % kspacing)
    if kgamma :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=T')
    else :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=F')
    
    sys_list = glob.glob('sys-*')
    for ss in sys_list:
        os.chdir(ss)
        ln_src = os.path.relpath(os.path.join(work_dir,'INCAR'))
        os.symlink(ln_src, 'INCAR')
        ln_src = os.path.relpath(os.path.join(work_dir,'POTCAR'))
        os.symlink(ln_src, 'POTCAR')
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
            pos_src = os.path.join(os.path.join(init_path, ii), 'CONTCAR')
            if not os.path.isfile(pos_src):
                if skip_relax :
                    pos_src = os.path.join(os.path.join(init_path, ii), 'POSCAR')
                    assert(os.path.isfile(pos_src))
                else :
                    raise RuntimeError("not file %s, vasp relaxation should be run before scale poscar")
            scale_path = os.path.join(work_path, ii)
            scale_path = os.path.join(scale_path, "sys-%.3f" % jj)
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
    
    cwd = os.getcwd()
    path_sp = os.path.join(out_dir, global_dirname_03)
    assert(os.path.isdir(path_sp))
    os.chdir(path_sp)
    sys_pe = glob.glob('sys-*')
    sys_pe.sort()
    os.chdir(cwd)    

    pert_cmd = cwd
    pert_cmd = os.path.join(pert_cmd, 'tools')
    pert_cmd = os.path.join(pert_cmd, 'create_random_disturb.py')
    pert_cmd += ' -etmax %f -ofmt vasp POSCAR %d %f > /dev/null' %(pert_box, pert_numb, pert_atom)
    for ii in sys_pe :
        for jj in scale :
            path_work = path_sp
            path_work = os.path.join(path_work, ii)
            path_work = os.path.join(path_work, 'sys-%.3f' % jj)
            assert(os.path.isdir(path_work))
            os.chdir(path_work)
            sp.check_call(pert_cmd, shell = True)
            for kk in range(pert_numb) :
                pos_in = 'POSCAR%d.vasp' % (kk+1)
                dir_out = '%06d' % (kk+1)
                create_path(dir_out)
                pos_out = os.path.join(dir_out, 'POSCAR')
                poscar_shuffle(pos_in, pos_out)
                os.remove(pos_in)
            kk = -1
            pos_in = 'POSCAR'
            dir_out = '%06d' % (kk+1)
            create_path(dir_out)
            pos_out = os.path.join(dir_out, 'POSCAR')
            poscar_shuffle(pos_in, pos_out)
            os.chdir(cwd)

def make_vasp_md(jdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    scale = jdata['scale']    
    encut = jdata['encut']
    kspacing = jdata['kspacing']
    kgamma = jdata['kgamma']
    pert_numb = jdata['pert_numb']
    md_temp = jdata['md_temp']
    md_nstep = jdata['md_nstep']

    cwd = os.getcwd()
    vasp_dir = os.path.join(cwd, 'vasp.in')
    vasp_dir = os.path.join(cwd, vasp_dir)
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
    shutil.copy2(os.path.join(vasp_dir, 'INCAR.md'), 
                 os.path.join(path_md, 'INCAR'))
    out_potcar = os.path.join(path_md, 'POTCAR')
    with open(out_potcar, 'w') as outfile:
        for fname in potcars:
            with open(fname) as infile:
                outfile.write(infile.read())
    os.chdir(path_md)
    replace('INCAR', 'ENCUT=.*', 'ENCUT=%f' % encut)
    replace('INCAR', 'ISIF=.*', 'ISIF=2')
    replace('INCAR', 'KSPACING=.*', 'KSPACING=%f' % kspacing)
    if kgamma :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=T')
    else :
        replace('INCAR', 'KGAMMA=.*', 'KGAMMA=F')    
    replace('INCAR', 'NSW=.*', 'NSW=%d' % md_nstep)
    replace('INCAR', 'TEBEG=.*', 'TEBEG=%d' % md_temp)
    replace('INCAR', 'TEEND=.*', 'TEEND=%d' % md_temp)
    os.chdir(cwd)    

    for ii in sys_ps :
        for jj in scale :
            for kk in range(pert_numb) :
                path_work = path_md
                path_work = os.path.join(path_work, ii)
                path_work = os.path.join(path_work, "sys-%.3f" % jj)
                path_work = os.path.join(path_work, "%06d" % kk)
                create_path(path_work)
                os.chdir(path_work)                
                path_pos = path_ps
                path_pos = os.path.join(path_pos, ii)
                path_pos = os.path.join(path_pos, "sys-%.3f" % jj)
                path_pos = os.path.join(path_pos, "%06d" % kk)
                init_pos = os.path.join(path_pos, 'POSCAR')
                shutil.copy2 (init_pos, 'POSCAR')
                file_incar = os.path.join(path_md, 'INCAR')
                file_potcar = os.path.join(path_md, 'POTCAR')
                os.symlink(os.path.relpath(file_incar), 'INCAR')
                os.symlink(os.path.relpath(file_potcar), 'POTCAR')
                os.chdir(cwd)                
                
def _main() :
    parser = argparse.ArgumentParser(
        description="gen init confs")
    parser.add_argument('PARAM', type=str, 
                        help="parameter file, json format")
    parser.add_argument('STAGE', type=int,
                        help="the stage of init, can be 1, 2 or 3. "
                        "1: setup vasp jobs for relaxation. "
                        "2: perturb system. "
                        "3: setup vasp jobs for MD of perturbed system. ")
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)
    out_dir = out_dir_name(jdata)
    jdata['out_dir'] = out_dir
    print ("# working dir %s" % out_dir)

    stage = args.STAGE

    if stage == 1 :
        create_path(out_dir)
        make_unit_cell(jdata)
        make_super_cell(jdata)
        place_element(jdata)
        make_vasp_relax(jdata)
    elif stage == 2 :
        make_scale(jdata)
        pert_scaled(jdata)
    elif stage == 3 :
        make_vasp_md(jdata)
    else :
        raise RuntimeError("unknow stage %d" % stage)

# # 00
# make_unit_cell(jdata)
# # 01
# make_super_cell(jdata)
# # 02
# place_element(jdata)
# make_vasp_relax(jdata)
# # 03
# make_scale(jdata)
# pert_scaled(jdata)

# poscar_shuffle('POSCAR', 'POSCAR.out')
# poscar_scale ('POSCAR', 'POSCAR.out', 2)

# sp.check_call("./tools/copy.py -n 2 2 2 POSCAR POSCAR.out", shell = True)
    
if __name__ == "__main__":
    _main()

#!/usr/bin/env python3 

import os,json,shutil,re,glob,argparse
import numpy as np
import subprocess as sp
import tools.hcp as hcp
import tools.fcc as fcc
import tools.diamond as diamond
import tools.sc as sc
import tools.bcc as bcc
from pymatgen.core.surface import SlabGenerator,generate_all_slabs, Structure

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
1       make unit cell
        copy to make super cell
        place element
        make vasp relax
1a      vasp relax
2       scale system
        perturb system
3       make vasp md
3a      vasp md
4       collect md data
"""
global_dirname_02 = '00.place_ele'
global_dirname_03 = '01.scale_pert'
global_dirname_04 = '02.md'

def out_dir_name(jdata) :
    cell_type = jdata['cell_type']
    elements = jdata['elements']
    super_cell = jdata['super_cell']    

    ele_str = "surf."
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

def _poscar_natoms(lines) :
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    return numb_atoms

def poscar_natoms(poscar_in) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
        return _poscar_natoms(lines)

def poscar_shuffle(poscar_in, poscar_out) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    numb_atoms = _poscar_natoms(lines)
    idx = np.arange(8, 8+numb_atoms)
    np.random.shuffle(idx)
    out_lines = lines[0:8]
    for ii in range(numb_atoms) :
        out_lines.append(lines[idx[ii]])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(out_lines))

def poscar_scale_direct (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = _poscar_natoms(lines)
    pscale = float(lines[1])
    pscale = pscale * scale
    lines[1] = str(pscale) + "\n"
    return lines

def poscar_scale_cartesian (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = _poscar_natoms(lines)
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

def poscar_elong (poscar_in, poscar_out, elong) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    if lines[7][0] != 'C' :
        raise RuntimeError("only works for Cartesian POSCAR")
    sboxz = lines[4].split()
    boxz = np.array([float(sboxz[0]), float(sboxz[1]), float(sboxz[2])])
    boxzl = np.linalg.norm(boxz)
    elong_ratio = elong / boxzl
    boxz = boxz * (1. + elong_ratio)
    lines[4] = '%.16e %.16e %.16e\n' % (boxz[0],boxz[1],boxz[2])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(lines))        

def make_super_cell (jdata) :
    latt = jdata['latt']
    out_dir = jdata['out_dir']
    lmp_cmd = jdata['lmp_cmd']
    path_sc = os.path.join(out_dir, global_dirname_02)
    cell_type = jdata['cell_type']
    z_min = jdata['z_min']
    super_cell = jdata['super_cell']

    cwd = os.getcwd()    
    path_surf = os.path.join(cwd, 'surf')
    path_surf = os.path.join(path_surf, cell_type)
    lmp_input = glob.glob(os.path.join(path_surf, "in.*[0-9]"))
    pcpy_cmd = os.path.join(cwd, 'tools')
    pcpy_cmd = os.path.join(pcpy_cmd, 'poscar_copy.py')
    cvt_cmd = os.path.join(cwd, 'tools')
    cvt_cmd = os.path.join(cvt_cmd, 'ovito_file_convert.py')
    
    path_work = create_path(path_sc)    
    path_work = os.path.abspath(path_work)
    os.chdir(path_work)
    for ii in lmp_input:
        surf_idx = ii.split('.')[-1]
        path_cur_surf = create_path('surf-'+surf_idx)
        os.chdir(path_cur_surf)
        if os.path.isfile('in.lmp') :
            os.remove('in.lmp')
        shutil.copy2(ii, 'in.lmp')
        replace('in.lmp', 'LATT', str(latt))
        tmp_cmd = lmp_cmd + " -i in.lmp" 
        proc = sp.Popen(tmp_cmd, shell = True,
                        stdout=sp.PIPE, 
                        stderr=sp.PIPE) 
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        for jj in stdout.split('\n'): 
            if "atom z-dim" in jj:
                z_dim = float(jj.split(' ')[2])
        nz = 1
        while z_dim * nz < z_min:
            nz += 1        
        # print(nz, z_dim, nz * z_dim)
        sp.check_call(cvt_cmd + ' -m vasp dump.atom POSCAR.unit', 
                      shell = True)
        sp.check_call(pcpy_cmd + ' -n %d %d %d POSCAR.unit POSCAR ' % (super_cell[0], super_cell[1], nz), 
                      shell = True)
        os.chdir(path_work)
    os.chdir(cwd)        

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

def make_super_cell_pymatgen (jdata) :
    make_unit_cell(jdata)

    out_dir = jdata['out_dir']
    path_uc = os.path.join(out_dir, global_dirname_02)
    from_path = path_uc
    from_file = os.path.join(from_path, 'POSCAR.unit')
    ss = Structure.from_file(from_file)

    all_millers = jdata['millers']
    path_sc = os.path.join(out_dir, global_dirname_02)
    z_min = jdata['z_min']
    super_cell = jdata['super_cell']

    cwd = os.getcwd()    
    path_work = (path_sc)    
    path_work = os.path.abspath(path_work)
    pcpy_cmd = os.path.join(cwd, 'tools')
    pcpy_cmd = os.path.join(pcpy_cmd, 'poscar_copy.py')
    os.chdir(path_work)
    for miller in all_millers:
        miller_str=""
        for ii in miller :
            miller_str += str(ii)        
        path_cur_surf = create_path('surf-'+miller_str)
        os.chdir(path_cur_surf)
        slabgen = SlabGenerator(ss, miller, z_min, 1e-3)
        all_slabs = slabgen.get_slabs() 
        print("Miller %s: The slab has %s termination, use the first one" %(str(miller), len(all_slabs)))
        all_slabs[0].to('POSCAR', 'POSCAR')
        if super_cell[0] > 1 or super_cell[1] > 1 :
            sp.check_call(pcpy_cmd + ' -n %d %d %d POSCAR POSCAR ' % (super_cell[0], super_cell[1], 1), 
                          shell = True)
        os.chdir(path_work)
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
    elements = jdata['elements']
    path_sc = os.path.join(out_dir, global_dirname_02)
    path_pe = os.path.join(out_dir, global_dirname_02)    
    path_sc = os.path.abspath(path_sc)
    path_pe = os.path.abspath(path_pe)
    
    assert(os.path.isdir(path_sc))
    assert(os.path.isdir(path_pe))
    cwd = os.getcwd()
    os.chdir(path_sc)
    surf_list = glob.glob('surf-*')
    surf_list.sort()
    os.chdir(cwd)

    for ss in surf_list:
        path_surf = os.path.join(path_sc, ss)        
        pos_in = os.path.join(path_surf, 'POSCAR')
        natoms = poscar_natoms(pos_in)
        combines = np.array(make_combines(len(elements), natoms), dtype = int)
        for ii in combines :
            if any(ii == 0) :
                continue
            comb_name = "sys-"
            for idx,jj in enumerate(ii) :            
                comb_name += "%04d" % jj
                if idx != len(ii)-1 :
                    comb_name += "-"
            path_work = os.path.join(path_surf, comb_name)
            create_path(path_work)
            pos_out = os.path.join(path_work, 'POSCAR')
            poscar_ele(pos_in, pos_out, elements, ii)
            poscar_shuffle(pos_out, pos_out)

def make_vasp_relax (jdata) :
    out_dir = jdata['out_dir']
    potcars = jdata['potcars']
    encut = jdata['encut']
    kspacing = jdata['kspacing_relax']
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
    
    sys_list = glob.glob(os.path.join('surf-*', 'sys-*'))
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
    cvt_cmd = os.path.join(cwd, 'tools')
    cvt_cmd = os.path.join(cvt_cmd, 'ovito_file_convert.py -m vasp ')
    init_path = os.path.join(out_dir, global_dirname_02)
    init_path = os.path.abspath(init_path)
    work_path = os.path.join(out_dir, global_dirname_03)
    os.chdir(init_path)
    init_sys = glob.glob(os.path.join('surf-*', 'sys-*'))
    init_sys.sort()
    os.chdir(cwd)

    create_path(work_path)
    for ii in init_sys :
        for jj in scale :
            pos_cont = os.path.join(os.path.join(init_path, ii), 'CONTCAR')
            if not os.path.isfile(pos_cont):
                if skip_relax :
                    pos_cont = os.path.join(os.path.join(init_path, ii), 'POSCAR')
                    assert(os.path.isfile(pos_cont))
                else :
                    raise RuntimeError("not file %s, vasp relaxation should be run before scale poscar")
            pos_src = os.path.join(os.path.join(init_path, ii), 'POSCAR.rlxed')
            sp.check_call(cvt_cmd + ' ' + pos_cont + ' ' + pos_src, shell = True)
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
    vacuum_resol = jdata['vacuum_resol']
    vacuum_max = jdata['vacuum_max']
    
    cwd = os.getcwd()
    path_sp = os.path.join(out_dir, global_dirname_03)
    assert(os.path.isdir(path_sp))
    path_sp = os.path.abspath(path_sp)
    os.chdir(path_sp)
    sys_pe = glob.glob(os.path.join('surf-*', 'sys-*'))
    sys_pe.sort()
    os.chdir(cwd)    

    pert_cmd = cwd
    pert_cmd = os.path.join(pert_cmd, 'tools')
    pert_cmd = os.path.join(pert_cmd, 'create_random_disturb.py')
    pert_cmd += ' -etmax %f -ofmt vasp POSCAR %d %f > /dev/null' %(pert_box, pert_numb, pert_atom)    
    for ii in sys_pe :
        for jj in scale :
            path_scale = path_sp
            path_scale = os.path.join(path_scale, ii)
            path_scale = os.path.join(path_scale, 'scale-%.3f' % jj)
            assert(os.path.isdir(path_scale))
            os.chdir(path_scale)
            poscar_in = os.path.join(path_scale, 'POSCAR')
            assert(os.path.isfile(poscar_in))
            for ll in np.arange(vacuum_resol, vacuum_max, vacuum_resol) :
                path_elong = path_scale
                path_elong = os.path.join(path_elong, 'elong-%3.3f' % ll) 
                create_path(path_elong)
                os.chdir(path_elong)
                poscar_elong(poscar_in, 'POSCAR', ll)                
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
                
def _main() :
    parser = argparse.ArgumentParser(
        description="gen init confs")
    parser.add_argument('PARAM', type=str, 
                        help="parameter file, json format")
    parser.add_argument('STAGE', type=int,
                        help="the stage of init, can be 1 or 2 "
                        "1: Setup vasp jobs for relaxation. "
                        "2: Collect vasp relaxed confs (if relax is not skiped). Perturb system. "
    )
    args = parser.parse_args()

    with open (args.PARAM, 'r') as fp :
        jdata = json.load (fp)
    out_dir = out_dir_name(jdata)
    jdata['out_dir'] = out_dir
    pymatgen_surf = False 
    if 'pymatgen_surf' in jdata:
        pymatgen_surf = jdata['pymatgen_surf']
    print ("# working dir %s" % out_dir)

    stage = args.STAGE

    if stage == 1 :
        create_path(out_dir)
        if pymatgen_surf :
            make_super_cell_pymatgen(jdata)
        else :
            make_super_cell(jdata)
        place_element(jdata)
        make_vasp_relax(jdata)
    # elif stage == 0 :
    #     # create_path(out_dir)
    #     # make_super_cell(jdata)
    #     # place_element(jdata)
    #     # make_vasp_relax(jdata)
    #     # make_scale(jdata)
    #     # pert_scaled(jdata)
    #     # poscar_elong('POSCAR', 'POSCAR.out', 3)
    #     pert_scaled(jdata)
    elif stage == 2 :
        make_scale(jdata)
        pert_scaled(jdata)
    else :
        raise RuntimeError("unknow stage %d" % stage)
    
if __name__ == "__main__":
    _main()

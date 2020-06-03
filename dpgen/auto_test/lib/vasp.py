#!/usr/bin/python3
import os
import warnings
import numpy as np
import dpgen.auto_test.lib.lammps as lammps
import dpgen.auto_test.lib.util as util
from dpgen.generator.lib.vasp import incar_upper
from pymatgen.io.vasp import Incar,Kpoints,Potcar

class OutcarItemError(Exception):
    pass

# def get_poscar(conf_dir) :
#     conf_path = os.path.abspath(conf_dir)
#     poscar_out = os.path.join(conf_path, 'POSCAR')
#     contcar = os.path.join(conf_path, 'POSCAR')
#     lmp_dump = glob.glob(os.path.join(conf_path, 'dump.*'))
#     lmp_dump.sort()
#     if os.path.isfile(poscar_out) :
#         pass
#     elif os.path.isfile(contcar) :
#         os.symlink(contcar, poscar_out)
#     elif len(lmp_dump) == 1 :
#         dump_file = lmp_dump[0]
#         lammps.poscar_from_last_dump(dump_file, task_poscar, deepmd_type_map)

def regulate_poscar(poscar_in, poscar_out) :
    with open(poscar_in, 'r') as fp:
        lines = fp.read().split('\n')
    names = lines[5].split()
    counts = [int(ii) for ii in lines[6].split()]
    assert(len(names) == len(counts))
    uniq_name = []
    for ii in names :
        if not (ii in uniq_name) :
            uniq_name.append(ii)
    uniq_count = np.zeros(len(uniq_name), dtype = int)
    for nn,cc in zip(names,counts) :
        uniq_count[uniq_name.index(nn)] += cc
    natoms = np.sum(uniq_count)
    posis = lines[8:8+natoms]
    all_lines = []
    for ele in uniq_name:
        ele_lines = []
        for ii in posis :
            ele_name = ii.split()[-1]
            if ele_name == ele :
                ele_lines.append(ii)
        all_lines += ele_lines
    all_lines.append('')
    ret = lines[0:5]
    ret.append(" ".join(uniq_name))
    ret.append(" ".join([str(ii) for ii in uniq_count]))
    ret.append("Direct")
    ret += all_lines
    with open(poscar_out, 'w') as fp:
        fp.write("\n".join(ret))

def sort_poscar(poscar_in, poscar_out, new_names) :
    with open(poscar_in, 'r') as fp:
        lines = fp.read().split('\n')
    names = lines[5].split()
    counts = [int(ii) for ii in lines[6].split()]
    new_counts = np.zeros(len(counts), dtype = int)
    for nn,cc in zip(names,counts) :
        new_counts[new_names.index(nn)] += cc
    natoms = np.sum(new_counts)
    posis = lines[8:8+natoms]
    all_lines = []
    for ele in new_names:
        ele_lines = []
        for ii in posis :
            ele_name = ii.split()[-1]
            if ele_name == ele :
                ele_lines.append(ii)
        all_lines += ele_lines
    all_lines.append('')
    ret = lines[0:5]
    ret.append(" ".join(new_names))
    ret.append(" ".join([str(ii) for ii in new_counts]))
    ret.append("Direct")
    ret += all_lines
    with open(poscar_out, 'w') as fp:
        fp.write("\n".join(ret))

def perturb_xz (poscar_in, poscar_out, pert = 0.01) :
    with open(poscar_in, 'r') as fp:
        lines = fp.read().split('\n')
    zz = lines[4]
    az = [float(ii) for ii in zz.split()]
    az[0] += pert
    zz = [str(ii) for ii in az]
    zz = " ".join(zz)
    lines[4] = zz
    with open(poscar_out, 'w') as fp:
        fp.write("\n".join(lines))

def reciprocal_box(box) :
    rbox = np.linalg.inv(box)
    rbox = rbox.T
    # rbox = rbox / np.linalg.det(box)
    # print(np.matmul(box, rbox.T))
    # print(rbox)
    return rbox

def make_kspacing_kpoints(poscar, kspacing, kgamma) :
    if type(kspacing) is not list:
        kspacing = [kspacing, kspacing, kspacing]
    with open(poscar, 'r') as fp:
        lines = fp.read().split('\n')
    scale = float(lines[1])
    box = []
    for ii in range(2,5) :
        box.append([float(jj) for jj in lines[ii].split()[0:3]])
    box = np.array(box)
    box *= scale
    rbox = reciprocal_box(box)
    kpoints = [max(1,(np.ceil(2 * np.pi * np.linalg.norm(ii) / ks).astype(int))) for ii,ks in zip(rbox,kspacing)]
    ret = make_vasp_kpoints(kpoints, kgamma)
    return ret

def get_energies (fname) :
    if not check_finished(fname):
        warnings.warn("incomplete outcar: "+fname)
    with open(fname, 'r') as fp:
        lines = fp.read().split('\n')
    try :
        ener = _get_energies(lines)
        return ener
    except OutcarItemError :
        return None

def get_boxes (fname) :
    if not check_finished(fname):
        warnings.warn("incomplete outcar: "+fname)
    with open(fname, 'r') as fp:
        lines = fp.read().split('\n')
    try :
        ener = _get_boxes(lines)
        return ener
    except OutcarItemError :
        return None

def get_nev(fname) :
    if not check_finished(fname):
        warnings.warn("incomplete outcar: "+fname)
    with open(fname, 'r') as fp:
        lines = fp.read().split('\n')
    try:
        natoms = _get_natoms(lines)
        vol = _get_volumes(lines)[-1]
        ener = _get_energies(lines)[-1]
        return natoms, ener/natoms, vol/natoms
    except OutcarItemError:
        raise OutcarItemError("cannot find the result, please check the OUTCAR")
    # print(fname, natoms, vol, ener)

def get_stress(fname) :
    if not check_finished(fname):
        warnings.warn("incomplete outcar: "+fname)
    with open(fname, 'r') as fp:
        lines = fp.read().split('\n')
    try:
        stress = _get_stress(lines)[-1]
        return stress
    except OutcarItemError:
        return None

def check_finished(fname) :
    with open(fname, 'r') as fp:
        return 'Elapsed time (sec):' in fp.read()

def _get_natoms(lines) :
    ipt = None
    for ii in lines:
        if 'ions per type' in ii :
            ipt = [int(jj) for jj in ii.split()[4:]]
            return sum(ipt)
    raise OutcarItemError("cannot find item 'ions per type'")

def _get_energies(lines) :
    items = []
    for ii in lines:
        if 'free  energy   TOTEN' in ii:
            items.append(float (ii.split()[4]))
    if len(items) == 0:
        raise OutcarItemError("cannot find item 'free  energy   TOTEN'")
    return items

def _split_box_line(line) :
    return [float(line[0:16]), float(line[16:29]), float(line[29:42])]

def _get_boxes(lines) :
    items = []
    for idx,ii in enumerate(lines):
        tmp_box = []
        if 'direct lattice vectors' in ii :
            tmp_box.append(_split_box_line(lines[idx+1]))
            tmp_box.append(_split_box_line(lines[idx+2]))
            tmp_box.append(_split_box_line(lines[idx+3]))
            items.append(tmp_box)
    return np.array(items)

def _get_volumes(lines) :
    items = []
    for ii in lines:
        if 'volume of cell' in ii:
            items.append(float (ii.split()[4]))
    if len(items) == 0:
        raise OutcarItemError("cannot find item 'volume of cell'")
    return items

def _get_stress(lines) :
    items = []
    for ii in lines:
        if 'in kB' in ii:
            sv = [float(jj) for jj in ii.split()[2:8]]
            tmp = sv[4]
            sv[4] = sv[5]
            sv[5] = tmp
            items.append(util.voigt_to_stress(sv))
    if len(items) == 0:
        raise OutcarItemError("cannot find item 'in kB'")
    return items

def _compute_isif (relax_ions,
                  relax_shape,
                  relax_volume) :
    if   (relax_ions) and (not relax_shape) and (not relax_volume) :
        isif = 2
    elif (relax_ions) and (relax_shape) and (relax_volume) :
        isif = 3
    elif (relax_ions) and (relax_shape) and (not relax_volume) :
        isif = 4
    elif (not relax_ions) and (relax_shape) and (not relax_volume) :
        isif = 5
    elif (not relax_ions) and (relax_shape) and (relax_volume) :
        isif = 6
    elif (not relax_ions) and (not relax_shape) and (relax_volume) :
        isif = 7
    else :
        raise ValueError("unknow relax style")
    return isif

def make_vasp_static_incar (ecut, ediff,
                            npar, kpar,
                            kspacing = 0.5, kgamma = True,
                            ismear = 1, sigma = 0.2) :
    isif = 2
    ret = ''
    ret += 'PREC=A\n'
    ret += 'ENCUT=%d\n' % ecut
    ret += '# ISYM=0\n'
    ret += 'ALGO=normal\n'
    ret += 'EDIFF=%e\n' % ediff
    ret += 'EDIFFG=-0.01\n'
    ret += 'LREAL=A\n'
    ret += 'NPAR=%d\n' % npar
    ret += 'KPAR=%d\n' % kpar
    ret += "\n"
    ret += 'ISMEAR=%d\n' % ismear
    ret += 'SIGMA=%f\n' % sigma
    ret += "\n"
    ret += 'ISTART=0\n'
    ret += 'ICHARG=2\n'
    ret += 'NELMIN=6\n'
    ret += 'ISIF=%d\n' % isif
    ret += 'IBRION=-1\n'
    ret += "\n"
    ret += 'NSW=0\n'
    ret += "\n"
    ret += 'LWAVE=F\n'
    ret += 'LCHARG=F\n'
    ret += 'PSTRESS=0\n'
    ret += "\n"
    if kspacing is not None :
        ret += 'KSPACING=%f\n' % kspacing
    if kgamma is not None :
        if kgamma:
            ret += 'KGAMMA=T\n'
        else :
            ret += 'KGAMMA=F\n'
    return ret

def make_vasp_relax_incar (ecut, ediff,
                           relax_ion, relax_shape, relax_volume,
                           npar, kpar,
                           kspacing = 0.5, kgamma = True,
                           ismear = 1, sigma = 0.22) :
    isif = _compute_isif(relax_ion, relax_shape, relax_volume)
    ret = ''
    ret += 'PREC=A\n'
    ret += 'ENCUT=%d\n' % ecut
    ret += '# ISYM=0\n'
    ret += 'ALGO=normal\n'
    ret += 'EDIFF=%e\n' % ediff
    ret += 'EDIFFG=-0.01\n'
    ret += 'LREAL=A\n'
    ret += 'NPAR=%d\n' % npar
    ret += 'KPAR=%d\n' % kpar
    ret += "\n"
    ret += 'ISMEAR=%d\n' % ismear
    ret += 'SIGMA=%f\n' % sigma
    ret += "\n"
    ret += 'ISTART=0\n'
    ret += 'ICHARG=2\n'
    ret += 'NELM=100\n'
    ret += 'NELMIN=6\n'
    ret += 'ISIF=%d\n' % isif
    ret += 'IBRION=2\n'
    ret += "\n"
    ret += 'NSW=50\n'
    ret += "\n"
    ret += 'LWAVE=F\n'
    ret += 'LCHARG=F\n'
    ret += 'PSTRESS=0\n'
    ret += "\n"
    if kspacing is not None :
        ret += 'KSPACING=%f\n' % kspacing
    if kgamma is not None :
        if kgamma:
            ret += 'KGAMMA=T\n'
        else :
            ret += 'KGAMMA=F\n'
    return ret

def make_vasp_phonon_incar (ecut, ediff,
                            npar, kpar,
                            kspacing = 0.5, kgamma = True,
                            ismear = 1, sigma = 0.2) :
    isif = 2
    ret = ''
    ret += 'PREC=A\n'
    ret += 'ENCUT=%d\n' % ecut
    ret += '# ISYM=0\n'
    ret += 'ALGO=normal\n'
    ret += 'EDIFF=%e\n' % ediff
    ret += 'EDIFFG=-0.01\n'
    ret += 'LREAL=A\n'
    #ret += 'NPAR=%d\n' % npar
    ret += 'KPAR=%d\n' % kpar
    ret += "\n"
    ret += 'ISMEAR=%d\n' % ismear
    ret += 'SIGMA=%f\n' % sigma
    ret += "\n"
    ret += 'ISTART=0\n'
    ret += 'ICHARG=2\n'
    ret += 'NELMIN=4\n'
    ret += 'ISIF=%d\n' % isif
    ret += 'IBRION=8\n'
    ret += "\n"
    ret += 'NSW=1\n'
    ret += "\n"
    ret += 'LWAVE=F\n'
    ret += 'LCHARG=F\n'
    ret += 'PSTRESS=0\n'
    ret += "\n"
    if kspacing is not None :
        ret += 'KSPACING=%f\n' % kspacing
    if kgamma is not None :
        if kgamma:
            ret += 'KGAMMA=T\n'
        else :
            ret += 'KGAMMA=F\n'
    return ret

def get_poscar_types (fname) :
    with open(fname, 'r') as fp :
        lines = fp.read().split('\n')
    return lines[5].split()

def get_poscar_natoms (fname) :
    with open(fname, 'r') as fp :
        lines = fp.read().split('\n')
    return [int(ii) for ii in lines[6].split()]

def _poscar_natoms(lines) :
    numb_atoms = 0
    for ii in lines[6].split() :
        numb_atoms += int(ii)
    return numb_atoms

def _poscar_scale_direct (str_in, scale) :
    lines = str_in.copy()
    numb_atoms = _poscar_natoms(lines)
    pscale = float(lines[1])
    pscale = pscale * scale
    lines[1] = str(pscale) + "\n"
    return lines

def _poscar_scale_cartesian (str_in, scale) :
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

def poscar_natoms(poscar_in) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    return _poscar_natoms(lines)

def poscar_scale (poscar_in, poscar_out, scale) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    if 'D' == lines[7][0] or 'd' == lines[7][0] :
        lines = _poscar_scale_direct(lines, scale)
    elif 'C' == lines[7][0] or 'c' == lines[7][0] :
        lines = _poscar_scale_cartesian(lines, scale)
    else :
        raise RuntimeError("Unknow poscar coord style at line 7: %s" % lines[7])
    with open(poscar_out, 'w') as fout:
        fout.write("".join(lines))

def poscar_vol (poscar_in) :
    with open(poscar_in, 'r') as fin :
        lines = list(fin)
    box = []
    for ii in range(2,5) :
        words = lines[ii].split()
        vec = [float(jj) for jj in words]
        box.append(vec)
    scale = float(lines[1].split()[0])
    box = np.array(box)
    box *= scale
    return np.linalg.det(box)

def _make_vasp_kp_gamma(kpoints):
    ret = ""
    ret += "Automatic mesh\n"
    ret += "0\n"
    ret += "Gamma\n"
    ret += "%d %d %d\n" % (kpoints[0], kpoints[1], kpoints[2])
    ret += "0  0  0\n"
    return ret

def _make_vasp_kp_mp(kpoints):
    ret = ""
    ret += "K-Points\n"
    ret += " 0\n"
    ret += "Monkhorst Pack\n"
    ret += "%d %d %d\n" % (kpoints[0], kpoints[1], kpoints[2])
    ret += " 0  0  0\n"
    return ret

def make_vasp_kpoints (kpoints, kgamma = False) :
    if kgamma :
        ret = _make_vasp_kp_gamma(kpoints)
    else :
        ret = _make_vasp_kp_mp(kpoints)
    return ret


def make_vasp_kpoints_from_incar(work_dir,jdata):
    cwd=os.getcwd()
    fp_aniso_kspacing = jdata.get('fp_aniso_kspacing')
    os.chdir(work_dir)
    # get kspacing and kgamma from incar
    assert(os.path.exists('INCAR'))
    with open('INCAR') as fp:
        incar = fp.read()
    standard_incar = incar_upper(Incar.from_string(incar))
    if fp_aniso_kspacing is None:
        try:
            kspacing = standard_incar['KSPACING']
        except KeyError:
            raise RuntimeError ("KSPACING must be given in INCAR")
    else :
        kspacing = fp_aniso_kspacing
    try:
        gamma = standard_incar['KGAMMA']
        if isinstance(gamma,bool):
            pass
        else:
            if gamma[0].upper()=="T":
                gamma=True
            else:
                gamma=False
    except KeyError:
        raise RuntimeError ("KGAMMA must be given in INCAR")
    # check poscar
    assert(os.path.exists('POSCAR'))
    # make kpoints
    ret=make_kspacing_kpoints('POSCAR', kspacing, gamma)
    kp=Kpoints.from_string(ret)
    kp.write_file("KPOINTS")
    os.chdir(cwd)

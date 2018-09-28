#!/usr/bin/python3 

import os
import numpy as np

def system_from_poscar(poscar) :
    lines = open(poscar, 'r').read().split('\n')
    system = {}
    system['atom_names'] = [str(ii) for ii in lines[5].split()]
    system['atom_numbs'] = [int(ii) for ii in lines[6].split()]
    scale = float(lines[1])
    cell = []
    for ii in range(2,5) :
        boxv = [float(jj) for jj in lines[ii].split()]
        boxv = np.array(boxv) * scale
        cell.append(boxv)
    system['cell'] = np.array(cell)
    natoms = sum(system['atom_numbs'])
    coord = []
    for ii in range(8, 8+natoms) :
        tmpv = [float(jj) for jj in lines[ii].split()]
        tmpv = np.array(tmpv) * scale
        coord.append(tmpv)
    system['coordinates'] = np.array(coord)
    return system

def make_vasp_kpoints (kpoints) :
    ret = ""
    ret += "Automatic mesh\n"
    ret += "0\n"
    ret += "Gamma\n"
    ret += "%d %d %d\n" % (kpoints[0], kpoints[1], kpoints[2])
    ret += "0  0  0\n"
    return ret

def _make_vasp_incar (ecut, ediff, npar, kpar, kspacing = 0.5, kgamma = True) :
    ret = ''
    ret += 'PREC=A\n'
    ret += 'ENCUT=%d\n' % ecut
    ret += 'ISYM=0\n'
    ret += 'ALGO=fast\n'
    ret += 'EDIFF=%e\n' % ediff
    ret += 'LREAL=A\n'
    ret += 'NPAR=%d\n' % npar
    ret += 'KPAR=%d\n' % kpar
    ret += "\n"
    ret += 'NELMIN=4\n'
    ret += 'ISIF=2\n'
    ret += '# ISMEAR=1\n'
    ret += '# SIGMA=0.25\n'
    ret += 'IBRION=-1\n'
    ret += "\n"
    ret += 'NSW=0\n'
    ret += "\n"
    ret += 'LWAVE=F\n'
    ret += 'LCHARG=F\n'
    ret += 'PSTRESS=0\n'
    ret += "\n"
    ret += 'KSPACING=%f\n' % kspacing
    if kgamma:
        ret += 'KGAMMA=.TRUE.\n'
    else :
        ret += 'KGAMMA=.FALSE.\n'
    return ret

def make_vasp_incar(fp_params) :
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    incar = _make_vasp_incar(ecut, ediff, npar, kpar, kspacing = kspacing, kgamma = False)
    return incar    
    
def make_vasp_kpoints_gamma (kpoints) :
    ret = ''
    ret += 'Automatic mesh\n'
    ret += '0\n'
    ret += 'Gamma\n'
    ret += '%d %d %d\n' % (kpoints[0], kpoints[1], kpoints[2])
    ret += '0  0  0\n'
    return ret

def make_vasp_kpoints (kpoints) :
    return make_vasp_kpoints_gamma(kpoints)

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

def _make_vasp_incar (ecut, ediff, npar, kpar, 
                      kspacing = 0.5, kgamma = True, 
                      smearing = None, sigma = None, 
                      metagga = None) :
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
    if smearing is not None :
        ret += 'ISMEAR=%d\n' % smearing
    if sigma is not None :
        ret += 'SIGMA=%f\n' % sigma
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
    if metagga is not None :
        ret += '\n'
        ret += 'LASPH=T\n'
        ret += 'METAGGA=%s\n' % metagga
    return ret


def _make_vasp_incar_dict (ecut, ediff, npar, kpar, 
                           kspacing = 0.5, kgamma = True, 
                           smearing = None, sigma = None, 
                           metagga = None) :
    incar_dict = {}
    incar_dict['PREC'] = 'A'
    incar_dict['ENCUT'] = ecut
    incar_dict['ISYM'] = 0
    incar_dict['ALGO'] = 'fast'
    incar_dict['EDIFF'] = ediff
    incar_dict['LREAL'] = 'A'
    incar_dict['NPAR'] = npar
    incar_dict['KPAR'] = kpar
    incar_dict['NELMIN'] = 4
    incar_dict['ISIF'] = 2
    if smearing is not None :
        incar_dict['ISMEAR'] = smearing        
    if sigma is not None :
        incar_dict['SIGMA'] = sigma
    incar_dict['IBRION'] = -1
    incar_dict['NSW'] = 0
    incar_dict['LWAVE'] = 'F'
    incar_dict['LCHARG'] = 'F'
    incar_dict['PSTRESS'] = 0
    incar_dict['KSPACING'] = kspacing
    if kgamma:
        incar_dict['KGAMMA'] = 'T'
    else :
        incar_dict['KGAMMA'] = 'F'
    if metagga is not None :
        incar_dict['LASPH'] = 'T'
        incar_dict['METAGGA'] = metagga
    return incar_dict


def _update_incar_dict(incar_dict_, user_dict) :
    if user_dict is None:
        return incar_dict_
    incar_dict = incar_dict_
    for ii in user_dict :
        ci = ii.upper()
        incar_dict[ci] = user_dict[ii]
    return incar_dict

def _write_incar_dict(incar_dict) :
    lines = []
    for key in incar_dict:
        lines.append('%s=%s' % (key, incar_dict[key]))
    return '\n'.join(lines)


def _make_smearing(fp_params) :
    smearing = None
    sigma = None
    if 'smearing' in fp_params :
        smearing = fp_params['smearing']
    if 'sigma' in fp_params :
        sigma = fp_params['sigma']
    if smearing == None :
        return None, sigma
    smearing_method = (smearing.split(':')[0]).lower()
    if smearing_method == 'mp' :
        order = 1
        if len(smearing.split(':')) == 2 :
            order = int(smearing.split(':')[1])
        return order, sigma
    elif smearing_method == 'gauss' :
        return 0, sigma
    elif smearing_method == 'fd' :
        return -1, sigma
    else :
        raise RuntimeError("unsuppported smearing method %s " % smearing_method)

def _make_metagga(fp_params) :
    metagga = None
    if 'metagga' in fp_params :
        metagga = fp_params['metagga']
    if metagga == 'NONE':
        metagga = None
    elif metagga not in [None,'SCAN', 'TPSS', 'RTPSS', 'M06L', 'MBJ'] :
        raise RuntimeError ("unknown metagga method " + metagga) 
    return metagga

def make_vasp_incar(fp_params) :
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    smearing, sigma = _make_smearing(fp_params)
    metagga = _make_metagga(fp_params)
    incar = _make_vasp_incar(ecut, ediff, npar, kpar, 
                             kspacing = kspacing, kgamma = False, 
                             smearing = smearing, sigma = sigma, 
                             metagga = metagga
    )
    return incar    
    
def make_vasp_incar_user_dict(fp_params) :
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    npar = fp_params['npar']
    kpar = fp_params['kpar']
    kspacing = fp_params['kspacing']
    if 'user_vasp_params' in fp_params :
        user_dict = fp_params['user_vasp_params']
    else :
        user_dict = None
    smearing, sigma = _make_smearing(fp_params)
    metagga = _make_metagga(fp_params)
    incar_dict = _make_vasp_incar_dict(ecut, ediff, npar, kpar, 
                                       kspacing = kspacing, kgamma = False, 
                                       smearing = smearing, sigma = sigma, 
                                       metagga = metagga
    )
    incar_dict = _update_incar_dict(incar_dict, user_dict)
    incar = _write_incar_dict(incar_dict)
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


if __name__ == '__main__' :
    import json
    jdata = json.load(open('param.json'))
    incar_1 = make_vasp_incar(jdata['fp_params'])
    incar_2 = make_vasp_incar_user_dict(jdata['fp_params'])
    with open('tmp1.out', 'w') as fp:
        fp.write(incar_1)
    with open('tmp2.out', 'w') as fp:
        fp.write(incar_2)


#!/usr/bin/python3 

import os
import numpy as np
from pymatgen.io.vasp import Incar

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

def write_incar_dict(incar_dict) :
    lines = []
    for key in incar_dict:
        if (type(incar_dict[key]) == bool):
            if incar_dict[key]:
                rs = 'T'
            else :
                rs = 'F'
        else :
            rs = str(incar_dict[key])
        lines.append('%s=%s' % (key, rs))
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
    incar = write_incar_dict(incar_dict)
    return incar
    
def incar_upper(dincar):
    standard_incar={}
    for key,val in dincar.items():
        standard_incar[key.upper()]=val
    return Incar(standard_incar)

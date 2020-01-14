#!/usr/bin/python3 

import os
import numpy as np

def _reciprocal_box(box) :
    rbox = np.linaly.inv(box)
    rbox = rbox.T

def _make_pwmat_kp_mp(kpoints) :
    ret = ""
    ret += "%d %d %d 0 0 0 " % (kpoints[0], kpoints[1], kpoints[2])
    return ret

def _make_kspacing_kpoints(config, kspacing) :
    with open(config, 'r') as fp:
        lines = fp.read().split('\n')
    box = []
    for idx, ii in enumerate(lines):
        if 'lattice' in ii or 'Lattice' in ii or 'LATTICE' in ii:
            for kk in range(idx+1,idx+1+3):
                vector=[float(jj) for jj in lines[kk].split()[0:3]]
                box.append(vector)
            box = np.array(box)
            rbox = _reciprocal_box(box)
    kpoints = [(np.ceil(2 * np.pi * np.linalg.norm(ii) / kspacing).astype(int)) for ii in rbox]
    ret = _make_pwmat_kp_mp(kpoints)
    return ret


def _make_pwmat_input_dict (node1, node2, atom_config, ecut, e_error, 
                            rho_error, icmix = None, smearing = None,
                            sigma = None,kspacing = 0.5, xcfunctional = None,
                            convergence = None, flag_symm = None) :
    input_dict = {}
    input_dict['node1'] = node1
    input_dict['node2'] = node2
    input_dict['IN.ATOM'] = atom_config
    input_dict['ECUT'] = ecut
    input_dict['E_ERROR'] = e_eroor
    input_dict['RHO_ERROR'] = rho_error
    if icmix is not None:
        if sigma is not None:
            if smearing is not None:
                SCF_ITER0_1 = "6 4 3 0.0000 " + str(sigma) + " " + str(smearing)
                SCF_ITER0_2 = "94 4 3 " + str(icmix) + " " + str(sigma) + " " + str(smearing)
            else:
                SCF_ITER0_1 = "6 4 3 0.0000 " + str(simga) + " 2" 
                SCF_ITER0_2 = "94 4 3 " + str(icmix) + " " + str(simga) + " 2" 

        else:
            if smearing is not None:
                SCF_ITER0_1 = "6 4 3 0.0000 0.025 " + str(smearing)
                SCF_ITER0_2 = "94 4 3 " + str(icmix) + " 0.025 " + str(smearing)
            else:
                SCF_ITER0_1 = "6 4 3 0.0000 0.025 2" 
                SCF_ITER0_2 = "94 4 3 " + str(icmix) + " 0.025 2" 
    else:
        if sigma is not None:
            if smearing is not None:
                SCF_ITER0_1 = "6 4 3 0.0000 " + str(sigma) + " " + str(smearing)
                SCF_ITER0_2 = "94 4 3 1.0000 " + str(sigma) + " " + str(smearing)
            else:
                SCF_ITER0_1 = "6 4 3 0.0000 " + str(sigma) + " 2"
                SCF_ITER0_2 = "94 4 3 1.0000 " + str(sigma) + " 2"
        else:
            if smearing is not None:
                SCF_ITER0_1 = "6 4 3 0.0000 0.025 " + str(smearing)
                SCF_ITER0_2 = "94 4 3 1.0000 0.025 " + str(smearing)
            else:
                SCF_ITER0_1 = "6 4 3 0.0000 0.025 2"
                SCF_ITER0_2 = "94 4 3 1.0000 0.025 2"
    input_dict['SCF_ITER0_1'] = SCF_ITER0_1
    input_dict['SCF_ITER0_2'] = SCF_ITER0_2
    if xcfunctional is not None:
        input_dict['XCFUNCTIONAL'] = xcfunctional
    if convergence is not None:
        input_dict['CONVERGENCE'] = convergence
    if flag_symm is not None :
        MP_N123 = _make_kspacing_kpoints(atom.config, kspacing)
        MP_N123 += str(flag_symm)
    else:
        MP_N123 = _make_kspacing_kpoints(atom.config, kspacing)
    input_dict['MP_N123'] = MP_N123
    input_dict['OUT.WG'] = 'F'
    input_dict['OUT.RHO'] = 'F'
    input_dict['OUT.MLMD'] = 'T'
    return input_dict

def _update_input_dict(input_dict_, user_dict) :
    if user_dict is None:
        return input_dict_
    input_dict = input_dict_
    for ii in user_dict :
        input_dict[ci] = user_dict[ii]
    return input_dict

def write_input_dict(input_dict) :
    lines = []
    for key in input_dict:
        if (type(input_dict[key]) == bool):
            if input_dict[key]:
                rs = 'T'
            else :
                rs = 'F'
        else :
            rs = str(input_dict[key])
        lines.append('%s=%s' % (key, rs))
    return '\n'.join(lines)


def _make_smearing(fp_params) :
    icmix = None
    smearing = None
    sigma = None
    if 'icmix' in fp_params :
        icmix = fp_params['icmix']
    if 'smearing' in fp_params :
        smearing = fp_params['smearing']
    if 'sigma' in fp_params :
        sigma = fp_params['sigma']
    if icmix == None:
        if smearing == None:
            if sigma == None:
                return None, None, None
            else:
                return None, None, sigma
        else:
            if sigma == None:
                return None, smearing, None
            else:
                return None, smearing, sigma
    else:
        if smearing == None:
            if sigma == None:
                return icmix, None, None
            else:
                return icmix, None, sigma
        else:
            if sigma == None:
                return icmix, smearing, None
            else:
                return icmix, smearing, sigma

def _make_xcfunctional(fp_params) :
    xcfunctional = None
    if 'xcfunctional' in fp_params :
        xcfunctional = fp_params['xcfunctional']
    if xcfunctional == 'NONE':
        xcfunctional = None
    elif xcfunctional not in [None,'LDA', 'PBE', 'HSE', 'PBESOL', 'PW91', 'TPSS', 'SCAN'] :
        raise RuntimeError ("unknown xcfunctional " + xcfunctional) 
    return xcfunctional

def _make_convergence(fp_params) :
    convergence = None
    if 'convergence' in fp_params :
        convergence = fp_params['convergence']
    if convergence == 'NONE' :
        convergence = None
    elif convergence not in [None, 'EASY', 'DIFFICULT'] :
        raise RuntimeError ("unknow convergence type " + convergence)
    return convergence

def _make_flag_symm(fp_params) :
    flag_symm = None
    if 'flag_symm' in fp_params :
        flag_symm = fp_params['flag_symm']
    if flag_symm == 'NONE' :
        flag_symm = None
    elif str(flag_symm) not in [None, '0', '1', '2', '3'] :
        raise RuntimeError ("unknow flag_symm type " + str(flag_symm))
    return flag_symm

def make_pwmat_input_user_dict(fp_params) :
    node1 = fp_params['node1']
    node2 = fp_params['node2']
    atom_config = fp_params['in.atom']
    ecut = fp_params['ecut']
    e_error = fp_params['e_error']
    rho_error = fp_params['rho_error']
    kspacing = fp_params['kspacing']
    if 'user_pwmat_params' in fp_params :
        user_dict = fp_params['user_pwmat_params']
    else :
        user_dict = None
    icmix, smearing, sigma = _make_smearing(fp_params)
    xcfunctional = _make_xcfunctional(fp_params)
    convergence = _make_convergence(fp_params)
    flag_symm = _make_flag_symm(fp_params)
    input_dict = _make_pwmat_input_dict(node1, node2, atom_config, ecut, e_error,
                                        rho_error, icmix = icmix, smearing = smearing,
                                        sigma = sigma, kspacing = kspacing, 
                                        xcfunctional = xcfunctional, convergence = convergence,
                                        flag_symm = flag_symm
    )
    input_dict = _update_input_dict(input_dict, user_dict)
    input = write_input_dict(input_dict)
    return input
    
def input_upper(dinput):
    standard_input={}
    for key,val in dinput.items():
        standard_input[key.upper()]=val
    return Input(standard_input)

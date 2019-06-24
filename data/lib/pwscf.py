#!/usr/bin/python3 

import numpy as np
# from lib.vasp import system_from_poscar

def _convert_dict(idict) :
    lines = []
    for key in idict.keys() :
        if type(idict[key]) == bool:
            if idict[key] :
                ws = '.TRUE.'
            else:
                ws = '.FALSE.'
        elif type(idict[key]) == str:
            ws = '\'' + idict[key] + '\''
        else :
            ws = str(idict[key])
        lines.append('%s=%s,' % (key, ws))
    return lines


def make_pwscf_01_runctrl_dict(sys_data, idict) :
    tot_natoms = sum(sys_data['atom_numbs'])
    ntypes = len(sys_data['atom_names'])
    lines = []
    lines.append('&control')
    lines += _convert_dict(idict['control'])
    lines.append('pseudo_dir=\'./\',')
    lines.append('/')
    lines.append('&system')
    lines += _convert_dict(idict['system'])
    lines.append('ibrav=0,')
    lines.append('nat=%d,' % tot_natoms)
    lines.append('ntyp=%d,' % ntypes)
    lines.append('/')
    if 'electrons' in idict :
        lines.append('&electrons')
        lines += _convert_dict(idict['electrons'])
    lines.append('/')
    lines.append('')
    return '\n'.join(lines)
    

def _make_pwscf_01_runctrl(sys_data, ecut, ediff, smearing, degauss) :
    tot_natoms = sum(sys_data['atom_numbs'])
    ntypes = len(sys_data['atom_names'])
    ret = ""
    ret += '&control\n'
    ret += "calculation='scf',\n"
    ret += "restart_mode='from_scratch',\n"
    ret += "outdir='./OUT',\n"
    ret += "tprnfor=.TRUE.,\n"
    ret += "tstress=.TRUE.,\n"
    ret += "disk_io='none',\n"
    ret += "pseudo_dir='./',\n"
    ret += "/\n"
    ret += "&system\n"
    ret += "vdw_corr='TS',\n"
    ret += "ecutwfc=%s,\n" % str(ecut)
    ret += "ts_vdw_econv_thr=%s,\n" % str(ediff)
    ret += "nosym=.TRUE.,\n"
    ret += "ibrav=0,\n"
    ret += "nat=%d,\n" % tot_natoms
    ret += "ntyp=%d,\n" % ntypes
    if degauss is not None :
        ret += 'degauss=%f,\n' % degauss
    if smearing is not None :
        ret += 'smearing=\'%s\',\n' % (smearing.lower())
    ret += "/\n"
    ret += "&electrons\n"
    ret += "conv_thr=%s,\n" % str(ediff)
    ret += "/\n"
    return ret

def _make_pwscf_02_species(sys_data, pps) :
    atom_names = (sys_data['atom_names'])
    if 'atom_masses' in sys_data:
        atom_masses = (sys_data['atom_masses'])
    else :
        atom_masses = [1 for ii in atom_names]
    ret = ""
    ret += "ATOMIC_SPECIES\n"
    ntypes = len(atom_names)
    assert(ntypes == len(atom_names))
    assert(ntypes == len(atom_masses))
    assert(ntypes == len(pps))
    for ii in range(ntypes) :
        ret += "%s %d %s\n" % (atom_names[ii], atom_masses[ii], pps[ii])
    return ret
        
def _make_pwscf_03_config(sys_data) :
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3,3])
    coordinates = sys_data['coords'][0]
    atom_names = (sys_data['atom_names'])
    atom_numbs = (sys_data['atom_numbs'])
    ntypes = len(atom_names)
    ret = ""
    ret += "CELL_PARAMETERS { angstrom }\n"
    for ii in range(3):
        for jj in range(3):
            ret += "%f " % cell[ii][jj]
        ret += "\n"
    ret += "\n"
    ret += "ATOMIC_POSITIONS { angstrom }\n"
    cc = 0
    for ii in range(ntypes):
        for jj in range(atom_numbs[ii]):            
            ret += "%s %f %f %f\n" % (atom_names[ii],
                                      coordinates[cc][0],
                                      coordinates[cc][1],
                                      coordinates[cc][2])
            cc += 1
    return ret

def _kshift(nkpt) :
    if (nkpt//2) * 2 == nkpt :
        return 1
    else :
        return 0
            
def _make_pwscf_04_kpoints(sys_data, kspacing):
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3,3])
    rcell = np.linalg.inv(cell)
    rcell = rcell.T
    kpoints = [(np.ceil(2 * np.pi * np.linalg.norm(ii) / kspacing).astype(int))
               for ii in rcell]
    ret = ""
    ret += "K_POINTS { automatic }\n"
    for ii in range(3) :
        ret += "%d " % kpoints[ii]
    for ii in range(3) :
        ret += "%d " % _kshift(kpoints[ii])
    ret += "\n"
    return ret

def _make_smearing(fp_params) :
    smearing = None
    degauss = None 
    if 'smearing' in fp_params :
        smearing = (fp_params['smearing']).lower()        
    if 'sigma' in fp_params :
        degauss = fp_params['sigma']
    if (smearing is not None) and (smearing.split(':')[0] == 'mp') :
        smearing = 'mp'
    if not (smearing in [None, 'gauss', 'mp', 'fd']) :
        raise RuntimeError("unknow smearing method " + smearing)
    return smearing, degauss

def make_pwscf_input(sys_data, fp_pp_files, fp_params, user_input = True) :
    if not user_input :
        ecut = fp_params['ecut']
        ediff = fp_params['ediff']
        smearing, degauss = _make_smearing(fp_params)
    kspacing = fp_params['kspacing']
    ret = ""
    if not user_input:
        ret += _make_pwscf_01_runctrl(sys_data, ecut, ediff, smearing, degauss)
    else :
        ret += make_pwscf_01_runctrl_dict(sys_data, fp_params)
    ret += "\n"
    ret += _make_pwscf_02_species(sys_data, fp_pp_files)
    ret += "\n"
    ret += _make_pwscf_03_config(sys_data)
    ret += "\n"
    ret += _make_pwscf_04_kpoints(sys_data, kspacing)
    ret += "\n"
    return ret

    

# sys_data = system_from_poscar('POSCAR')
# ret = ""
# ret += _make_pwscf_01_runctrl(sys_data, 20, 1e-8)
# ret += "\n"
# ret += _make_pwscf_02_species(sys_data, ['../pp/C_HSCV_PBE-1.0.UPF', '../H_HSCV_PBE-1.0.UPF', '../N_HSCV_PBE-1.0.UPF'])
# ret += "\n"
# ret += _make_pwscf_03_config(sys_data)
# ret += "\n"
# ret += _make_pwscf_04_kpoints(sys_data, 0.6)
# ret += "\n"

# open('input', 'w').write(ret)

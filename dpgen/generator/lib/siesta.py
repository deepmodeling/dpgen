import numpy as np
from dpdata.periodic_table import Element

def _make_siesta_01_common(sys_data, fp_params):
    tot_natoms = sum(sys_data['atom_numbs'])
    ntypes = len(sys_data['atom_names'])

    ret = ""
    ret += 'SystemName        system\n'
    ret += 'SystemLabel       system\n'
    ret += 'NumberOfAtoms     %d\n' % tot_natoms
    ret += 'NumberOfSpecies   %d\n' % ntypes
    ret += '\n'
    ret += 'WriteForces       T\n'
    ret += 'WriteCoorStep     T\n'
    ret += 'WriteCoorXmol     T\n'
    ret += 'WriteMDXmol       T\n'
    ret += 'WriteMDHistory    T\n\n'

    if 'ecut' in fp_params.keys():
        ecut = fp_params['ecut']
        ret += 'MeshCutoff            %s' % str(ecut)
        ret += ' Ry\n'
    if 'ediff' in fp_params.keys():
        ediff = fp_params['ediff']
        ret += 'DM.Tolerance          %e\n' % ediff
    if 'mixWeight' in fp_params.keys():
        mixingWeight = fp_params['mixingWeight']
        ret += 'DM.MixingWeight       %f\n' % mixingWeight
    if 'NumberPulay' in fp_params.keys():
        NumberPulay = fp_params['NumberPulay']
        ret += 'DM.NumberPulay         %d\n' % NumberPulay
    ret += 'DM.UseSaveDM          true\n'
    ret += 'XC.functional          GGA\n'
    ret += 'XC.authors             PBE\n'
    ret += 'MD.UseSaveXV           T\n\n'
    ret += 'DM.UseSaveDM           F\n'
    ret += 'WriteDM                F\n'
    ret += 'WriteDM.NetCDF         F\n'
    ret += 'WriteDMHS.NetCDF       F\n'
    return ret

def _make_siesta_02_species(sys_data, pps):
    atom_nums = sys_data['atom_numbs']
    atom_names = sys_data['atom_names']
    ntypes = len(atom_nums)
    assert (ntypes == len(atom_names))
    assert (ntypes == len(pps))
    ret = ''
    ret += '%block Chemical_Species_label\n'
    for i in range(0, len(atom_names)):
        ret += str(i + 1) + '\t' + str(Element(atom_names[i]).Z) + '\t' + atom_names[i] + '\n'
    ret += '%endblock Chemical_Species_label\n'
    return ret

# ## kpoints !!!
def _make_siesta_03_kpoint(sys_data, fp_param):
    if 'kspacing' in fp_param.keys():
        kspacing = fp_param['kspacing']
        cell = sys_data['cells'][0]
        cell = np.reshape(cell, [3, 3])
        rcell = np.linalg.inv(cell)
        rcell = rcell.T
        kpoints = [(np.ceil(2 * np.pi * np.linalg.norm(ii) / kspacing).astype(int))
                   for ii in rcell]

        ret = ""
        ret += '%block kgrid_Monkhorst_Pack\n'
        ret += '%d' % kpoints[0]
        ret += '\t0\t0\t0.0\n'

        ret += '0\t'
        ret += '%d' % kpoints[1]
        ret += '\t0\t0.0\n'

        ret += '0\t0\t'
        ret += '%d' % kpoints[2]
        ret += '\t0.0\n'

        ret += '%endblock kgrid_Monkhorst_Pack\n'
        return ret
    else:
        return ''

### coordinate
def _make_siesta_04_ucVectorCoord(sys_data):
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3, 3])
    coordinates = sys_data['coords'][0]
    atom_names = (sys_data['atom_names'])
    atom_numbs = (sys_data['atom_numbs'])
    ntypes = len(atom_names)
    ret = ""
    ret += "LatticeConstant 1.00 Ang\n"
    ret += "%block LatticeVectors\n"
    for ii in range(3):
        for jj in range(3):
            ret += "%f " % cell[ii][jj]
        ret += "\n"
    ret += "%endblock LatticeVectors\n"

    ret += "\n"
    ret += "AtomicCoordinatesFormat Ang\n"
    ret += "%block AtomicCoordinatesAndAtomicSpecies\n"
    cc = 0
    for ii in range(ntypes):
        for jj in range(atom_numbs[ii]):
            ret += "%f %f %f %d %s\n" % (coordinates[cc][0],
                                         coordinates[cc][1],
                                         coordinates[cc][2],
                                         ii + 1,
                                         atom_names[ii])
            cc += 1
    ret += "%endblock AtomicCoordinatesAndAtomicSpecies"
    return ret

def make_siesta_input(sys_data, fp_pp_files, fp_params):
    ret = ""
    ret += _make_siesta_01_common(sys_data, fp_params)
    ret += "\n"
    ret += _make_siesta_02_species(sys_data, fp_pp_files)
    ret += "\n"
    ret += _make_siesta_03_kpoint(sys_data, fp_params)
    ret += "\n"
    ret += _make_siesta_04_ucVectorCoord(sys_data)
    ret += "\n"
    return ret
    

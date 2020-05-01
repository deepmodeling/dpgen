import numpy as np
from dpdata.periodic_table import Element



def _make_siesta_01_common(sys_data, ecut, ediff, mixingWeight, NumberPulay):
    tot_natoms = sum(sys_data['atom_numbs'])
    ntypes = len(sys_data['atom_names'])
    ret = ""
    ret += 'SystemName        system\n'
    ret += 'SystemLabel       system\n'
    ret += 'NumberOfAtoms     %d' % tot_natoms
    ret += '\nNumberOfSpecies   %d\n' % ntypes
    ret += '\n'
    ret += 'WriteForces       T\n'
    ret += 'WriteCoorStep     T\n'
    ret += 'WriteCoorXmol     T\n'
    ret += 'WriteMDXmol       T\n'
    ret += 'WriteMDHistory    T\n\n'

    ret += 'MeshCutoff            %s' % str(ecut)
    ret += ' Ry\n'
    ret += 'DM.MixingWeight       %f\n' % mixingWeight
    ret += 'DM.Tolerance          %e\n' % ediff
    ret += 'DM.UseSaveDM          true\n'
    ret += 'DM.NumberPulay         %d\n' % NumberPulay
    ret += 'MD.UseSaveXV           T\n\n'

    ret += 'XC.functional          GGA\n'
    ret += 'XC.authors             PBE\n'
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


# ## kpoints !!! can not understand
def _make_siesta_03_kpoint(sys_data, kspacing):
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3, 3])
    ## np.linalg.inv()：矩阵求逆
    rcell = np.linalg.inv(cell)
    ## .T 矩阵转置
    rcell = rcell.T
    # np.ceil()是向上取整，与四舍五入无关 -5.6 --> -5
    # np.linalg.norm：进行范数运算，范数是对向量（或者矩阵）的度量，是一个标量（scalar）
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


### coordinate
def _make_siesta_04_ucVectorCoord(sys_data):
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3, 3])
    coordinates = sys_data['coords'][0]
    atom_names = (sys_data['atom_names'])
    atom_numbs = (sys_data['atom_numbs'])
    ntypes = len(atom_names)
    ret = ""
    ret += "LatticeConstant try_input_output.00 Ang\n"
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
    ecut = fp_params['ecut']
    ediff = fp_params['ediff']
    mixingWeight = fp_params['mixingWeight']
    NumberPulay = fp_params['NumberPulay']
    kspacing = fp_params['kspacing']
    ret = ""
    ret += _make_siesta_01_common(sys_data, ecut, ediff, mixingWeight, NumberPulay)
    ret += "\n"
    ret += _make_siesta_02_species(sys_data, fp_pp_files)
    ret += "\n"
    ret += _make_siesta_03_kpoint(sys_data, kspacing)
    ret += "\n"
    ret += _make_siesta_04_ucVectorCoord(sys_data)
    ret += "\n"
    return ret


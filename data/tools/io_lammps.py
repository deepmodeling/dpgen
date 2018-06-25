#!/usr/bin/env python3

"""

ASE Atoms convert to LAMMPS configuration
Some functions are adapted from ASE lammpsrun.py

"""

import numpy as np
from numpy.linalg import norm

import ase.io


def dir2car(v, A):
    """Direct to cartesian coordinates"""
    return np.dot(v, A)


def car2dir(v, Ainv):
    """Cartesian to direct coordinates"""
    return np.dot(v, Ainv)


def stress9_to_stress6(s9):
    # S6: xx yy zz yz xz xy
    s6 = np.zeros(6)
    s6[0] = s9[0][0]
    s6[1] = s9[1][1]
    s6[2] = s9[2][2]
    s6[3] = s9[1][2]
    s6[4] = s9[0][2]
    s6[5] = s9[0][1]
    return s6


def stress6_to_stress9(s6):
    # S6: xx yy zz yz xz xy
    s9 = np.zeros((3, 3))
    s9[0, :] = np.array([s6[0], s6[5], s6[4]])
    s9[1, :] = np.array([s6[5], s6[1], s6[3]])
    s9[2, :] = np.array([s6[4], s6[3], s6[2]])
    return s9


def is_upper_triangular(mat):
    """
    test if 3x3 matrix is upper triangular
    LAMMPS has a rule for cell matrix definition
    """
    def near0(x):
        """Test if a float is within .00001 of 0"""
        return abs(x) < 0.00001
    return near0(mat[1, 0]) and near0(mat[2, 0]) and near0(mat[2, 1])


def convert_cell(ase_cell):
    """
    Convert a parallel piped (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """

    # if ase_cell is lower triangular, cell is upper tri-angular
    cell = np.matrix.transpose(ase_cell) 

    if not is_upper_triangular(cell):
        # rotate bases into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = cell[:, 0]
        B = cell[:, 1]
        C = cell[:, 2]
        tri_mat[0, 0] = norm(A)
        Ahat = A / norm(A)
        AxBhat = np.cross(A, B) / norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(ase_cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans = trans / volume
        coord_transform = tri_mat * trans
        return tri_mat.T # return the lower-tri-angular
    else:
        return ase_cell


def convert_positions(pos0, cell0, cell_new, direct=False):
    if direct:
        pos = np.dot(pos0, cell_new)  # positions in new cellmatrix
    else:
        cell0_inv = np.linalg.inv(cell0)
        R = np.dot(cell_new, cell0_inv)
        pos = np.dot(pos0, R)
        '''
        print(R)
        print(R.T)
        print(np.linalg.inv(R))
        print(np.linalg.det(R))
        '''
    return pos


def convert_forces(forces0, cell0, cell_new):
    forces = convert_positions(forces0, cell0, cell_new, direct=False)
    return forces


def convert_stress(s6_0, cell0, cell_new):
    s9_0 = stress6_to_stress9(s6_0)

    cell0_inv = np.linalg.inv(cell0)
    R = np.dot(cell_new, cell0_inv)
    R_T = R.T

    s9 = np.dot(np.dot(R, s9_0), R_T)
    s6 = stress9_to_stress6(s9)
    return s6


def get_atoms_ntypes(atoms):
    atomic_numbers = atoms.numbers
    Uatomic_numbers = np.unique(atomic_numbers)  # unique atomic numbers
    ntypes = Uatomic_numbers.size
    return ntypes


def set_atoms_typeids(atoms):
    csymbols = atoms.get_chemical_symbols()
    U_csymbols = np.unique(csymbols)  # unique atomic numbers
    ntypes = U_csymbols.size
    typeids = {}
    for i in range(ntypes):
        typeids[U_csymbols[i]] = i + 1  # typeid start from 1
    return typeids


def set_atoms_typeids_with_atomic_numbers(atoms):
    # set typeids as atomic numbers, which can be robust when several
    # configuration with different atomic types were used
    # however, lammps do not allow it.
    csymbols = atoms.get_chemical_symbols()
    U_csymbols = np.unique(csymbols)  # unique atomic numbers
    ntypes = U_csymbols.size
    typeids = {}
    for i in range(ntypes):
        cs = U_csymbols[i]
        typeids[cs] = ase.data.atomic_numbers[cs]
    return typeids


def get_typeid(typeids, csymbol):
    return typeids[csymbol]


def ase2lammpsdata(atoms, typeids=None, fout='out.lmp'):
    # atoms: ase.Atoms
    # typeids: eg. {'Zr': 1, 'Nb': 2, 'Hf': 3}, should start with 1 and continuous
    # fout: output file name
    fw = open(fout, 'w')
    fw.write('# LAMMPS data written by PotGen-ASE\n')
    fw.write('\n')

    # write number of atoms
    natoms = atoms.get_number_of_atoms()
    fw.write('%d atoms\n' % natoms)
    fw.write('\n')

    # write number of types
    ntypes = get_atoms_ntypes(atoms)
    fw.write("%d atom types\n" % ntypes)
    fw.write('\n')

    # write cell information
    # transfer the cell into lammps' style
    cell0 = atoms.get_cell()
    cell = convert_cell(cell0)
    xhi = cell[0, 0]
    yhi = cell[1, 1]
    zhi = cell[2, 2]

    # low triangular matrix
    xy = cell[1, 0]
    xz = cell[2, 0]
    yz = cell[2, 1]

    fw.write("%f\t%f\t xlo xhi\n" % (0, xhi))
    fw.write("%f\t%f\t ylo yhi\n" % (0, yhi))
    fw.write("%f\t%f\t zlo zhi\n" % (0, zhi))
    fw.write("%f\t%f\t%f\t xy xz yz\n" % (xy, xz, yz))
    fw.write('\n')

    # write mases
    masses = np.unique(atoms.get_masses())
    fw.write('Masses\n')
    fw.write('\n')
    for i in range(ntypes):
        fw.write('%d\t%f\n' % (i + 1, masses[i]))
        fw.flush()
    fw.write('\n')

    # convert positions
    atoms.set_cell(cell, scale_atoms=True)
    pos = atoms.get_positions()
    '''
    pos0 = atoms.get_positions()
    pos = convert_positions(pos0, cell0, cell)   # positions in new cellmatrix
    '''
    # === Write postions ===
    fw.write('Atoms\n')
    fw.write('\n')
    symbols = atoms.get_chemical_symbols()
    if typeids is None:
        typeids = set_atoms_typeids(atoms)
    for i in range(natoms):
        cs = symbols[i]
        typeid = get_typeid(typeids, cs)  # typeid start from 1~N
        # typeid = ase.data.atomic_numbers[cs]  # typeid as their atomic
        # numbers
        fw.write('%d\t%d\t%f\t%f\t%f\n' %
                 (i + 1, typeid, pos[i][0], pos[i][1], pos[i][2]))
        fw.flush()
    fw.close()
    return

# test
if __name__ == "__main__":
    import sys
    fin = sys.argv[1]
    ATOMS = ase.io.read(fin)
    ase2lammpsdata(ATOMS)
    ase2lammpsdata(ATOMS, typeids={'Al': 1}, fout=fin + '.lmp')

    sep = '=' * 40

    pos0 = ATOMS.get_positions()
    cell0 = ATOMS.get_cell()
    print(cell0)
    print(pos0)
    print(sep)

    '''
    cell_new = np.eye(3)*4.05
    '''
    delta = 4.05 * 0.02

    '''
    cell_new = 4.05*np.array([[1, delta, 0],
                         [delta, 1, 0],
                         [0, 0, 1 / (1 - delta**2)]])
    '''
    cell_new = 4.05 * np.array([[1 + delta, 0, 0],
                                [0, 1 + delta, 0],
                                [0, 0, 1 / (1 + delta**2)]])

    pos = convert_positions(pos0, cell0, cell_new)
    print(cell0)
    print(cell_new)
    #print(np.linalg.det(cell0), np.linalg.det(cell_new))
    print(pos)
    print(sep)

    # anothother transoformation for LAMMPS low-tri-angle matrix
    cell_new = convert_cell(cell_new)
    pos = convert_positions(pos0, cell0, cell_new)
    print(cell0)
    print(cell_new)
    #print(np.linalg.det(cell0), np.linalg.det(cell_new))
    print(pos)
    print(sep)

    # test for stress tensor transformation
    stress0 = np.array([0.000593	, 0.000593,	0.000593	,
                        0.,	0.00,	0.00])
    stress_new = convert_stress(stress0, cell0, cell_new)
    print(stress0)
    print(stress_new)

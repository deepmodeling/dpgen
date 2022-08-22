#!/usr/bin/env python3

import sys
import os
import shutil
import glob
import argparse

import numpy as np
import ase.io
import dpgen.data.tools.io_lammps as io_lammps

from dpgen.generator.lib.abacus_scf import get_abacus_STRU, make_abacus_scf_stru

def create_disturbs_atomsk(fin, nfile, dmax=1.0, ofmt="lmp"):
    # removing the exists files
    flist = glob.glob('*.' + ofmt)
    for f in flist:
        os.remove(f)
    # Based on our tests, we find it always creates a disturb by
    # constant value of dmax for atomsk
    for i in range(1, nfile + 1):
        fout = fin + str(i) + '.' + ofmt
        cmd = "atomsk " + fin + " -disturb " + str(dmax) + " -wrap -ow " + fout
        os.system(cmd)
    return


def random_range(a, b, ndata=1):
    data = np.random.random(ndata) * (b - a) + a
    return data


def gen_random_disturb(dmax, a, b, dstyle='uniform'):
    d0 = np.random.rand(3) * (b - a) + a
    dnorm = np.linalg.norm(d0)
    if dstyle == 'normal':
        dmax = np.random.standard_normal(0, 0.5) * dmax
    elif dstyle == 'constant':
        pass
    else:
        # use if we just wanna a disturb in a range of [0, dmax),
        dmax = np.random.random() * dmax
    dr = dmax / dnorm * d0
    return dr


def create_disturbs_ase(fin, nfile, dmax=1.0, ofmt="lmp", dstyle='uniform', write_d=False):
    # removing the exists files
    flist = glob.glob('*.' + ofmt)
    for f in flist:
        os.remove(f)

    # read-in by ase
    atoms = ase.io.read(fin)
    natoms = atoms.get_number_of_atoms()
    pos0 = atoms.get_positions()

    # creat nfile ofmt files.
    for fid in range(1, nfile + 1):
        dpos = np.zeros((natoms, 3))
        atoms_d = atoms.copy()
        if write_d:
            fw = open('disp-' + str(fid) + '.dat', 'w')
        for i in range(natoms):
            # Use copy(), otherwise it will modify the input atoms every time.
            dr = gen_random_disturb(dmax, -0.5, 0.5, dstyle)
            '''
            if i == 1:
                print(dr)
                print(np.linalg.norm(dr))
            '''
            dpos[i, :] = dr
            if write_d:
                dnorm = np.linalg.norm(dr)
                fw.write('%d\t%f\t%f\t%f\t%f\n' %
                         (i + 1, dr[0], dr[1], dr[2], dnorm))
                fw.flush()
        pos = pos0 + dpos
        atoms_d.set_positions(pos)
        fout = fin + str(fid) + '.' + ofmt
        print("Creating %s ..." % fout)
        if ofmt in ['lmp', 'lammps_data']:
            # for lammps, use my personal output functions
            io_lammps.ase2lammpsdata(atoms_d, fout)
        else:
            ase.io.write(fout, atoms_d, ofmt, vasp5 = True)
        if write_d:
            fw.close()
    return


def gen_random_emat(etmax, diag=0):
    if np.abs(etmax) >= 1e-6:
        e = random_range(-etmax, etmax, 6)
    else:
        e = np.zeros(6)
    if diag != 0:
        # isotropic behavior
        e[3], e[4], e[5] = 0, 0, 0
    emat = np.array(
        [[e[0], 0.5 * e[5], 0.5 * e[4]],
         [0.5 * e[5], e[1], 0.5 * e[3]],
         [0.5 * e[4], 0.5 * e[3], e[2]]]
    )
    emat = emat + np.eye(3)
    return emat


def create_disturbs_ase_dev(fin, nfile, dmax=1.0, etmax=0.1, ofmt="lmp", dstyle='uniform', write_d=False, diag=0):
    # removing the exists files
    flist = glob.glob('*.' + ofmt)
    for f in flist:
        os.remove(f)

    # read-in by ase
    atoms = ase.io.read(fin)
    natoms = atoms.get_number_of_atoms()
    cell0 = atoms.get_cell()

    # creat nfile ofmt files.
    for fid in range(1, nfile + 1):
        # Use copy(), otherwise it will modify the input atoms every time.
        atoms_d = atoms.copy()

        # random flux for atomic positions
        if write_d:
            fw = open('disp-' + str(fid) + '.dat', 'w')
        dpos = np.zeros((natoms, 3))
        for i in range(natoms):
            dr = gen_random_disturb(dmax, -0.5, 0.5, dstyle)
            dpos[i, :] = dr
            if write_d:
                dnorm = np.linalg.norm(dr)
                fw.write('%d\t%f\t%f\t%f\t%f\n' %
                         (i + 1, dr[0], dr[1], dr[2], dnorm))
                fw.flush()

        # random flux for volumes
        cell = np.dot(cell0, gen_random_emat(etmax, diag))
        atoms_d.set_cell(cell, scale_atoms=True)
        if write_d:
            fout_c = 'cell-' + str(fid) + '.dat'
            np.savetxt(fout_c, cell, '%f')

        # determine new cell & atomic positions randomiziations
        pos = atoms_d.get_positions() + dpos
        atoms_d.set_positions(pos)

        # pre-converting the Atoms to be in low tri-angular cell matrix
        cell_new = io_lammps.convert_cell(cell)
        #pos_new = io_lammps.convert_positions(pos, cell, cell_new)
        atoms_d.set_cell(cell_new, scale_atoms=True)
        # atoms_d.set_positions(pos_new)

        # Writing it
        fout = fin + str(fid) + '.' + ofmt
        print("Creating %s ..." % fout)
        if ofmt in ['lmp', 'lammps_data']:
            # for lammps, use my personal output functions
            io_lammps.ase2lammpsdata(atoms_d, fout=fout)
        else:
            ase.io.write(fout, atoms_d, ofmt, vasp5 = True)
        if write_d:
            fw.close()
    return

def create_disturbs_abacus_dev(fin, nfile, dmax=1.0, etmax=0.1, ofmt="abacus", dstyle='uniform', write_d=False, diag=0):
    # removing the exists files
    flist = glob.glob('*.' + ofmt)
    for f in flist:
        os.remove(f)

    # read-in by ase
    #atoms = ase.io.read(fin)
    #natoms = atoms.get_number_of_atoms()
    #cell0 = atoms.get_cell()
    
    stru = get_abacus_STRU(fin)
    natoms = sum(stru["atom_numbs"])
    cell0 = stru['cells']

    # creat nfile ofmt files.
    for fid in range(1, nfile + 1):
        # Use copy(), otherwise it will modify the input atoms every time.
        stru_d = stru.copy()

        # random flux for atomic positions
        if write_d:
            fw = open('disp-' + str(fid) + '.dat', 'w')
        dpos = np.zeros((natoms, 3))
        for i in range(natoms):
            dr = gen_random_disturb(dmax, -0.5, 0.5, dstyle)
            dpos[i, :] = dr
            if write_d:
                dnorm = np.linalg.norm(dr)
                fw.write('%d\t%f\t%f\t%f\t%f\n' %
                         (i + 1, dr[0], dr[1], dr[2], dnorm))
                fw.flush()

        # random flux for volumes
        cell = np.dot(cell0, gen_random_emat(etmax, diag))
        stru_d['cells'] = cell
        if write_d:
            fout_c = 'cell-' + str(fid) + '.dat'
            np.savetxt(fout_c, cell, '%f')

        # determine new cell & atomic positions randomiziations
        stru_d['coords'] += dpos

        # pre-converting the Atoms to be in low tri-angular cell matrix
        cell_new = io_lammps.convert_cell(cell)
        #pos_new = io_lammps.convert_positions(pos, cell, cell_new)
        stru_d['cells'] = cell_new

        convert_mat = np.linalg.inv(cell).dot(cell_new)
        stru_d['coords'] = np.matmul(stru_d['coords'], convert_mat)


        # Writing it
        fout = fin + str(fid) + '.' + ofmt
        print("Creating %s ..." % fout)
        ret = make_abacus_scf_stru(stru_d, stru_d['pp_files'], stru_d['orb_files'], stru_d['dpks_descriptor'])
        with open(fout, "w") as fp:
            fp.write(ret)
        if write_d:
            fw.close()
    return


def create_random_alloys(fin, alloy_dist, ifmt='vasp', ofmt='vasp'):
    '''
    In fact, atomsk also gives us the convinient tool to do this
    '''
    # alloy_dist = {'Zr': 0.80, 'Nb': 0.20}
    atomic_symbols = alloy_dist.keys()
    atomic_ratios = alloy_dist.values()

    # renormalize the atomic_ratio
    atomic_ratios = atomic_ratios / atomic_ratios.sum()

    atoms = ase.io.read(fin, format=ifmt)

    # decide the exactly number of atoms for each types.
    natoms = atoms.num_of_numbers()
    num_for_each = []
    for r in atomic_ratios:
        num_for_each.append(int(r * natoms))

    # ========= decide which atom-ID to be substituted. ===========
    # a natoms vector with random sequence
    all_ids = np.random.permutation(natoms)
    ids_for_each = []
    for i in range(num_for_each):
        num = num_for_each[i]
        if i == 0:
            id_start = 0
        else:
            id_start = num_for_each[i - 1]
        ids_for_each.append(all_ids[id_start, id_start + num])

    # subsitute them, by replacing their chemical symbols of A by B
    new_chemical_symbols = atoms.get_chemical_symbols()
    for i in range(1, num_for_each):  # ignore the 0-th element (it is itself)
        symbol = atomic_symbols[i]
        ids = ids_for_each[i]
        for idx in ids:
            new_chemical_symbols[idx] = symbol
    # set the new chemical_symbols to the atoms
    atoms.set_chemical_symbols(new_chemical_symbols)

    # write it as ofmt
    fout = fin.split(',')[0] + '_random' + ofmt
    ase.io.write(fout, atoms, format=ofmt, vasp5 = True)
    return


def RandomDisturbParser():
    parser = argparse.ArgumentParser(
        description="Script to generate random disturb configurations")
    parser.add_argument('fin', type=str, help="input file name")
    parser.add_argument('nfile', type=int,
                        help='number of files to be created')
    parser.add_argument('dmax', type=float, help='dmax')
    parser.add_argument('-etmax', type=float, default=0,
                        help='etmax for random strain tensor generations')
    parser.add_argument('-diag', type=int, default=0,
                        help='only diagonal elements of strain tensors are randomized?')
    parser.add_argument('-ofmt', type=str, default='lmp',
                        help='output fileformat')
    parser.add_argument('-dstyle', type=str, default='uniform',
                        help='random distribution style [uniform?]')
    parser.add_argument('-wd', '--write_disp', type=int,
                        default=0, help='write displacement information?')
    return parser.parse_args()


#
if __name__ == "__main__":
    args = RandomDisturbParser()
    fin = args.fin
    nfile = args.nfile
    dmax = args.dmax
    ofmt = args.ofmt
    dstyle = args.dstyle
    write_d = args.write_disp
    etmax = args.etmax
    diag = args.diag

    if write_d == 0:
        write_d = False
    else:
        write_d = True

    # main program
    #create_disturbs_atomsk(fin, nfile, dmax, ofmt)
    #create_disturbs_ase(fin, nfile, dmax, ofmt, dstyle, write_d)
    if ofmt == "vasp":
        create_disturbs_ase_dev(fin, nfile, dmax, etmax,
                            ofmt, dstyle, write_d, diag)
    elif ofmt == "abacus":
        create_disturbs_abacus_dev(fin, nfile, dmax, etmax,
                            ofmt, dstyle, write_d, diag)
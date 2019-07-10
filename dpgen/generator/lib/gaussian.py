#!/usr/bin/python3


import uuid
import itertools
import numpy as np
import dpdata
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
try:
    import openbabel
except ImportError:
    pass
try:
    from ase import Atoms
except ImportError:
    pass


def _crd2frag(symbols, crds, pbc=False, cell=None):
    atomnumber = len(symbols)
    if pbc:
        all_atoms = Atoms(symbols = symbols, positions = crds, pbc=True, cell=cell)
        repeated_atoms = all_atoms.repeat(2)[atomnumber:]
        tree = cKDTree(crds)
        d = tree.query(repeated_atoms.get_positions(), k=1)[0]
        nearest = d < 5
        ghost_atoms = repeated_atoms[nearest]
        realnumber = np.where(nearest)[0] % atomnumber
        all_atoms += ghost_atoms
    xyzstring = ''.join((f"{atomnumber}\nDPGEN\n", "\n".join(
        ['{:2s} {:22.15f} {:22.15f} {:22.15f}'.format(s, x, y, z)
            for s, (x, y, z) in zip(symbols, crds)])))
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('xyz', 'mol2')
    mol = openbabel.OBMol()
    conv.ReadString(mol, xyzstring)
    mol2string = conv.WriteString(mol)
    bonds = []
    linecontent = -1
    for line in mol2string.split('\n'):
        if line.startswith("@<TRIPOS>BOND"):
            linecontent = 0
        else:
            if linecontent == 0:
                s = line.split()
                if len(s) > 3:
                    l = 12 if s[3] == 'ar' else int(s[3])
                    a, b = int(s[1])-1, int(s[2])-1
                    if a >= atomnumber and b >= atomnumber:
                        # duplicated
                        continue
                    elif a >= atomnumber:
                        a = realnumber[a-atomnumber]
                    elif b >= atomnumber:
                        b = realnumber[b-atomnumber]
                    bonds.append([a, b, l])
    bonds = np.array(bonds, ndmin=2).reshape((-1, 3))
    graph = csr_matrix(
        (bonds[:, 2], (bonds[:, 0], bonds[:, 1])), shape=(atomnumber, atomnumber))
    frag_numb, frag_index = connected_components(graph, 0)
    return frag_numb, frag_index


def detect_multiplicity(symbols):
    # only support C, H, O at present
    # oxygen -> 3
    if np.count_nonzero(symbols == ["O"]) == 2 and symbols.size == 2:
        return 3
    n_total = np.count_nonzero(symbols == ["H"])
    return n_total % 2 + 1


def make_gaussian_input(sys_data, fp_params):
    coordinates = sys_data['coords'][0]
    atom_names = sys_data['atom_names']
    #atom_numbs = sys_data['atom_numbs']
    atom_types = sys_data['atom_types']
    # get atom symbols list
    symbols = [atom_names[atom_type] for atom_type in atom_types]
    nproc = fp_params['nproc']
    keywords = fp_params['keywords']
    if type(keywords) == str:
        keywords = [keywords]
    # assume default charge is zero and default spin multiplicity is 1
    charge = fp_params['charge'] if 'charge' in fp_params else 0
    mult_auto = False
    frag = False
    if 'multiplicity' in fp_params:
        if type(fp_params['multiplicity']) == int:
            multiplicity = fp_params['multiplicity']
        elif fp_params['multiplicity'] == 'auto':
            mult_auto = True
        elif fp_params['multiplicity'] == 'frag':
            mult_auto = True
            frag = True
        else:
            raise RuntimeError('The keyword "multiplicity" is illegal.')
    else:
        multiplicity = 1

    if mult_auto:
        frag_numb, frag_index = _crd2frag(symbols, coordinates)
        if frag_numb == 1:
            frag = False
        mult_frags = []
        for i in range(frag_numb):
            idx = frag_index == i
            mult_frags.append(detect_multiplicity(np.array(symbols)[idx]))
        if frag:
            multiplicity = sum(mult_frags) - frag_numb + 1
            chargekeywords_frag = "%d %d" % (charge, multiplicity) + \
                ''.join([' %d %d' % (charge, mult_frag)
                         for mult_frag in mult_frags])
        else:
            multi_frags = np.array(mult_frags)
            multiplicity = 1 + \
                np.count_nonzero(multi_frags == 2) % 2 + \
                np.count_nonzero(multi_frags == 3) * 2
    buff = []
    # keywords, e.g., force b3lyp/6-31g**
    if frag:
        keywords[0] = '{} guess=fragment={}'.format(
            keywords[0], frag_numb)

    chkkeywords = []
    if len(keywords)>1:
        chkkeywords.append('%chk={}.chk'.format(str(uuid.uuid1())))

    nprockeywords = '%nproc={:d}'.format(nproc)
    titlekeywords = 'DPGEN'
    chargekeywords = '{} {}'.format(charge, multiplicity)

    buff = [*chkkeywords, nprockeywords, '#{}'.format(
        keywords[0]), '', titlekeywords, '', (chargekeywords_frag if frag else chargekeywords)]

    for ii, (symbol, coordinate) in enumerate(zip(symbols, coordinates)):
        if frag:
            buff.append("%s(Fragment=%d) %f %f %f" %
                        (symbol, frag_index[ii] + 1, *coordinate))
        else:
            buff.append("%s %f %f %f" % (symbol, *coordinate))
    if 'basis_set' in fp_params:
        # custom basis set
        buff.extend(['', fp_params['basis_set'], ''])
    for kw in itertools.islice(keywords, 1, None):
        buff.extend(['\n--link1--', *chkkeywords, nprockeywords,
                    '#{}'.format(kw), '', titlekeywords, '', chargekeywords, ''])
    buff.append('\n')
    return '\n'.join(buff)

def take_cluster(old_conf_name, type_map, idx, cutoff):
    sys = dpdata.System(old_conf_name, fmt = 'lammps/dump', type_map = type_map)
    atom_names = sys['atom_names']
    atom_types = sys['atom_types']
    cell = sys['cells'][0]
    coords = sys['coords'][0]
    symbols = [atom_names[atom_type] for atom_type in atom_types]
    # detect fragment 
    frag_numb, frag_index = _crd2frag(symbols, coords, True, cell)
    # get_distances
    all_atoms = Atoms(symbols = symbols, positions = coords, pbc=True, cell=cell)
    distances = all_atoms.get_distances(idx, range(len(all_atoms)), mic=True)
    cutoff_atoms_idx = np.where(distances < cutoff)[0]
    # make cutoff atoms in molecules
    taken_atoms_idx = []
    for ii in range(frag_numb):
        frag_atoms_idx = np.where(frag_index == ii)[0]
        if np.any(np.isin(frag_atoms_idx, cutoff_atoms_idx)):
            taken_atoms_idx.append(frag_atoms_idx)
    all_taken_atoms_idx = np.concatenate(taken_atoms_idx)
    # wrap
    cutoff_atoms = all_atoms[all_taken_atoms_idx]
    cutoff_atoms.wrap(
        center=coords[idx] /
        cutoff_atoms.get_cell_lengths_and_angles()[0: 3],
        pbc=True)
    coords = cutoff_atoms.get_positions()
    sys.data['coords'] = np.array([coords])
    sys.data['atom_types'] = atom_types[all_taken_atoms_idx]
    for ii, _ in enumerate(atom_names):
        sys.data['atom_numbs'][ii] = np.count_nonzero(sys.data['atom_types']==ii)
    return sys

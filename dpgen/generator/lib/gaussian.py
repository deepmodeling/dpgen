#!/usr/bin/python3


import uuid
import itertools
import warnings
import numpy as np
import dpdata
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
try:
    # expect openbabel >= 3.1.0
    from openbabel import openbabel
except ImportError:
    pass
try:
    from ase import Atoms, Atom
    from ase.data import atomic_numbers
except ImportError:
    pass


def _crd2frag(symbols, crds, pbc=False, cell=None, return_bonds=False):
    atomnumber = len(symbols)
    all_atoms = Atoms(symbols = symbols, positions = crds, pbc=pbc, cell=cell)
    # Use openbabel to connect atoms
    mol = openbabel.OBMol()
    mol.BeginModify()
    for idx, (num, position) in enumerate(zip(all_atoms.get_atomic_numbers(), all_atoms.positions)):
        atom = mol.NewAtom(idx)
        atom.SetAtomicNum(int(num))
        atom.SetVector(*position)
    # Apply period boundry conditions
    # openbabel#1853, supported in v3.1.0
    if pbc:
        uc = openbabel.OBUnitCell()
        uc.SetData(
            openbabel.vector3(cell[0][0], cell[0][1], cell[0][2]),
            openbabel.vector3(cell[1][0], cell[1][1], cell[1][2]),
            openbabel.vector3(cell[2][0], cell[2][1], cell[2][2]),
        )
        mol.CloneData(uc)
        mol.SetPeriodicMol()
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()
    mol.EndModify()
    bonds = []
    for ii in range(mol.NumBonds()):
        bond = mol.GetBond(ii)
        a = bond.GetBeginAtom().GetId()
        b = bond.GetEndAtom().GetId()
        bo = bond.GetBondOrder()
        bonds.extend([[a, b, bo], [b, a, bo]])
    bonds = np.array(bonds, ndmin=2).reshape((-1, 3))
    graph = csr_matrix(
        (bonds[:, 2], (bonds[:, 0], bonds[:, 1])), shape=(atomnumber, atomnumber))
    frag_numb, frag_index = connected_components(graph, 0)
    if return_bonds:
        return frag_numb, frag_index, graph
    return frag_numb, frag_index


def _crd2mul(symbols, crds):
    atomnumber = len(symbols)
    xyzstring = ''.join((f"{atomnumber}\nDPGEN\n", "\n".join(
        ['{:2s} {:22.15f} {:22.15f} {:22.15f}'.format(s, x, y, z)
            for s, (x, y, z) in zip(symbols, crds)])))
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats('xyz', 'gjf')
    mol = openbabel.OBMol()
    conv.ReadString(mol, xyzstring)
    gjfstring = conv.WriteString(mol)
    mul = int(gjfstring.split('\n')[5].split()[1])
    return mul  


def detect_multiplicity(symbols):
    # currently only support charge=0
    # oxygen -> 3
    if np.count_nonzero(symbols == ["O"]) == 2 and symbols.size == 2:
        return 3
    # calculates the total number of electrons, assumes they are paired as much as possible
    n_total = sum([atomic_numbers[s] for s in symbols])
    return n_total % 2 + 1


def make_gaussian_input(sys_data, fp_params):
    coordinates = sys_data['coords'][0]
    atom_names = sys_data['atom_names']
    #atom_numbs = sys_data['atom_numbs']
    atom_types = sys_data['atom_types']
    # get atom symbols list
    symbols = [atom_names[atom_type] for atom_type in atom_types]
    nproc = fp_params['nproc']

    if 'keywords_high_multiplicity' in fp_params and _crd2mul(symbols, coordinates)>=3:
        # multiplicity >= 3, meaning at least 2 radicals
        keywords = fp_params['keywords_high_multiplicity']
    else:
        keywords = fp_params['keywords']

    if type(keywords) == str:
        keywords = [keywords]
    else:
        keywords = keywords.copy()

    # assume default charge is zero and default spin multiplicity is 1
    if 'charge' in sys_data.keys():
        charge = sys_data['charge']
    else:
        charge = fp_params.get('charge', 0)
        
    use_fragment_guesses = False
    multiplicity = fp_params.get('multiplicity', 'auto')
    if type(multiplicity) == int:
        multiplicity = fp_params['multiplicity']
        mult_auto = False
    elif multiplicity == 'auto':
        mult_auto = True
    else:
        raise RuntimeError('The keyword "multiplicity" is illegal.')

    if fp_params.get("fragment_guesses", False):
        # Initial guess generated from fragment guesses
        # New feature of Gaussian 16
        use_fragment_guesses = True
        if not mult_auto:
            warnings.warn("Automatically set multiplicity to auto!")
            mult_auto = True

    if mult_auto:
        frag_numb, frag_index = _crd2frag(symbols, coordinates)
        if frag_numb == 1:
            use_fragment_guesses = False
        mult_frags = []
        for i in range(frag_numb):
            idx = frag_index == i
            mult_frags.append(detect_multiplicity(np.array(symbols)[idx]))
        if use_fragment_guesses:
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
    if use_fragment_guesses:
        keywords[0] = '{} guess=fragment={}'.format(
            keywords[0], frag_numb)

    chkkeywords = []
    if len(keywords)>1:
        chkkeywords.append('%chk={}.chk'.format(str(uuid.uuid1())))

    nprockeywords = '%nproc={:d}'.format(nproc)
    titlekeywords = 'DPGEN'
    chargekeywords = '{} {}'.format(charge, multiplicity)

    buff = [*chkkeywords, nprockeywords, '#{}'.format(
        keywords[0]), '', titlekeywords, '', (chargekeywords_frag if use_fragment_guesses else chargekeywords)]

    for ii, (symbol, coordinate) in enumerate(zip(symbols, coordinates)):
        if use_fragment_guesses:
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

def take_cluster(old_conf_name, type_map, idx, jdata):
    cutoff = jdata['cluster_cutoff']
    cutoff_hard = jdata.get('cluster_cutoff_hard', None)
    sys = dpdata.System(old_conf_name, fmt = 'lammps/dump', type_map = type_map)
    atom_names = sys['atom_names']
    atom_types = sys['atom_types']
    cell = sys['cells'][0]
    coords = sys['coords'][0]
    symbols = [atom_names[atom_type] for atom_type in atom_types]
    # detect fragment 
    frag_numb, frag_index, graph = _crd2frag(symbols, coords, True, cell, return_bonds=True)
    # get_distances
    all_atoms = Atoms(symbols = symbols, positions = coords, pbc=True, cell=cell)
    all_atoms[idx].tag = 1
    distances = all_atoms.get_distances(idx, range(len(all_atoms)), mic=True)
    distancescutoff = distances < cutoff
    cutoff_atoms_idx = np.where(distancescutoff)[0]
    if cutoff_hard is not None:
        distancescutoff_hard = distances < cutoff_hard
        cutoff_atoms_idx_hard = np.where(distancescutoff_hard)[0]
    # make cutoff atoms in molecules
    taken_atoms_idx = []
    added = []
    for ii in range(frag_numb):
        frag_atoms_idx = np.where(frag_index == ii)[0]
        if cutoff_hard is not None:
            # drop atoms out of the hard cutoff anyway
            frag_atoms_idx = np.intersect1d(frag_atoms_idx, cutoff_atoms_idx_hard)
        if np.any(np.isin(frag_atoms_idx, cutoff_atoms_idx)):
            if 'cluster_minify' in jdata and jdata['cluster_minify']:
                # support for organic species
                take_frag_idx=[]
                for aa in frag_atoms_idx:
                    if np.any(np.isin(aa, cutoff_atoms_idx)):
                        # atom is in the soft cutoff
                        # pick up anyway
                        take_frag_idx.append(aa)
                    elif np.count_nonzero(np.logical_and(distancescutoff, graph.toarray()[aa]==1)):
                        # atom is between the hard cutoff and the soft cutoff
                        # and has a single bond with the atom inside
                        if all_atoms[aa].symbol == 'H':
                            # for atom H: just add it
                            take_frag_idx.append(aa)
                        else:
                            # for other atoms (C, O, etc.): replace it with a ghost H atom
                            near_atom_idx = np.nonzero(np.logical_and(distancescutoff, graph.toarray()[aa]>0))[0][0]
                            vector = all_atoms[aa].position - all_atoms[near_atom_idx].position
                            new_position = all_atoms[near_atom_idx].position + vector / np.linalg.norm(vector) * 1.09
                            added.append(Atom('H', new_position))
                    elif np.count_nonzero(np.logical_and(distancescutoff, graph.toarray()[aa]>1)):
                        # if that atom has a double bond with the atom inside
                        # just pick up the whole fragment (within the hard cutoff)
                        # TODO: use a more fantastic method
                        take_frag_idx=frag_atoms_idx
                        break
            else:
                take_frag_idx = frag_atoms_idx
            taken_atoms_idx.append(take_frag_idx)
    all_taken_atoms_idx = np.concatenate(taken_atoms_idx)
    # wrap
    cutoff_atoms = sum(added, all_atoms[all_taken_atoms_idx])
    cutoff_atoms.wrap(
        center=coords[idx] /
        cutoff_atoms.get_cell_lengths_and_angles()[0: 3],
        pbc=True)
    coords = cutoff_atoms.get_positions()
    sys.data['coords'] = np.array([coords])
    sys.data['atom_types'] = np.array(list(atom_types[all_taken_atoms_idx]) + [atom_names.index('H')]*len(added))
    sys.data['atom_pref'] = np.array([cutoff_atoms.get_tags()])
    for ii, _ in enumerate(atom_names):
        sys.data['atom_numbs'][ii] = np.count_nonzero(sys.data['atom_types']==ii)
    return sys

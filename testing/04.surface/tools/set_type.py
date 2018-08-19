#!/usr/bin/env python3

import sys

def set_type (conf_file, type_id) :
    type_id = int(type_id)
    with open(conf_file, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines :
        if "atoms" in ii :
            natoms = int(ii.split()[0])
    new_lines = lines
    ntypes = type_id
    # revise ntypes
    idx_ntypes = -1
    for idx, ii in enumerate(lines) :
        if 'atom types' in ii :
            idx_ntypes = idx
    if idx_ntypes == -1 :
        raise RuntimeError("cannot find the entry 'atom types' in ", conf_file)
    words = lines[idx_ntypes].split()
    words[0] = str(ntypes)
    new_lines[idx_ntypes] = " ".join(words)
    # find number of atoms
    idx_atom_entry = -1
    for idx, ii in enumerate(lines) :
        if 'Atoms' in ii :
            idx_atom_entry = idx
    if idx_atom_entry == -1 :
        raise RuntimeError("cannot find the entry 'Atoms' in ", conf_file)
    # revise atom type
    for idx in range(idx_atom_entry+2, idx_atom_entry+2+natoms) :
        ii = lines[idx]
        words = ii.split()
        assert(len(words) >= 5)
        new_id = type_id
        words[1] = str(new_id)
        ii = " ".join(words)
        new_lines[idx] = ii
    # set masses
    new_lines.append("Masses")
    new_lines.append("")
    for ii in range(type_id) :
        new_lines.append("%d 1" % (ii+1))
    new_lines.append("")
    with open(conf_file, 'w') as fp:
        fp.write("\n".join(new_lines))

set_type(sys.argv[1], sys.argv[2])

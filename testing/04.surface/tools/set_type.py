#!/usr/bin/env python3

import sys

def set_type (conf_file, type_id) :
    with open(conf_file, 'r') as fp:
        lines = fp.read().split('\n')
    for ii in lines :
        if "atoms" in ii :
            natoms = int(ii.split()[0])
    # find number of atoms
    idx_atom_entry = -1
    for idx, ii in enumerate(lines) :
        if 'Atoms' in ii :
            idx_atom_entry = idx
    if idx_atom_entry == -1 :
        raise RuntimeError("cannot find the entry 'Atoms' in ", conf_file)
    # revise atom type
    new_lines = lines
    for idx in range(idx_atom_entry+2, idx_atom_entry+2+natoms) :
        ii = lines[idx]
        words = ii.split()
        assert(len(words) >= 5)
        new_id = type_id
        words[1] = str(new_id)
        ii = " ".join(words)
        new_lines[idx] = ii
    with open(conf_file, 'w') as fp:
        fp.write("\n".join(new_lines))

set_type(sys.argv[1], sys.argv[2])

#!/usr/bin/env python3

import numpy as np

ry2ev = 13.605693009
bohr2ang = 0.52917721067
kbar2evperang3 = 1./1602

def get_block (lines, keyword, skip = 0) :
    ret = []
    for idx,ii in enumerate(lines) :
        if keyword in ii :
            blk_idx = idx + 1 + skip
            while len(lines[blk_idx]) != 0 and blk_idx != len(lines):
                ret.append(lines[blk_idx])
                blk_idx += 1
            break
    return ret

def get_types (lines) :
    ret = []
    blk = get_block(lines, 'ATOMIC_SPECIES')
    for ii in blk:
        ret.append(ii.split()[0])
    return ret

def get_cell (lines) :
    ret = []
    blk = get_block(lines, 'CELL_PARAMETERS')
    for ii in blk:
        ret.append([float(jj) for jj in ii.split()[0:3]])
    ret = np.array(ret)
    return ret

def get_coords (lines) :
    ret = []
    blk = get_block(lines, 'ATOMIC_POSITIONS')
    for ii in blk:
        ret.append([float(jj) for jj in ii.split()[1:4]])
    ret = np.array(ret)
    return ret

def get_natoms (lines) :
    types = get_types (lines)
    names = []
    blk = get_block(lines, 'ATOMIC_POSITIONS')
    for ii in blk:
        names.append(ii.split()[0])
    natoms = []
    for ii in types :
        natoms.append(names.count(ii))
    return natoms

def get_energy (lines) :
    for ii in lines :
        if '!    total energy' in ii :
            return ry2ev * float(ii.split('=')[1].split()[0])
    return None

def get_force (lines) :
    blk = get_block(lines, 'Forces acting on atoms', skip = 1)
    ret = []
    for ii in blk:
        ret.append([float(jj) for jj in ii.split('=')[1].split()])
    ret = np.array(ret)
    ret *= (ry2ev / bohr2ang)
    return ret

def get_stress (lines) :
    blk = get_block(lines, 'total   stress')
    ret = []
    for ii in blk:
        ret.append([float(jj) for jj in ii.split()[3:6]])
    ret = np.array(ret)
    ret *= kbar2evperang3
    return ret
    
def write_config(types, natoms, coord, cell, energy, force, stress) :
    ret = ""
    ret += "#N %d 1\n" % sum(natoms)
    ret += "#C "
    for ii in types :
        ret += str(ii) + " "
    ret += "\n"
    ret += "##\n"
    ret += "#X %20.12e %20.12e %20.12e\n" % (cell[0][0], cell[0][1], cell[0][2])
    ret += "#Y %20.12e %20.12e %20.12e\n" % (cell[1][0], cell[1][1], cell[1][2])
    ret += "#Z %20.12e %20.12e %20.12e\n" % (cell[2][0], cell[2][1], cell[2][2])
    ret += "#W 1.0\n"
    ret += "#E %.12e\n" % (energy / (np.sum(natoms)))
    ret += "#S %.12e %.12e %.12e %.12e %.12e %.12e\n" %  \
           (stress[0][0], stress[1][1], stress[2][2], \
            stress[0][1], stress[1][2], stress[0][2])
    cc = 0
    ret += "#F\n"
    for idx,ii in enumerate(natoms):
        for jj in range(ii) :
            ret += "%d %.12e %.12e %.12e %.12e %.12e %.12e\n" % \
                   (idx, 
                    coord[cc][0], coord[cc][1], coord[cc][2],
                    force[cc][0], force[cc][1], force[cc][2]
                   )
            cc += 1
    return ret


def cvt_1frame (fin, fout, fconfig):
    outlines = open(fout, 'r').read().split('\n')
    inlines = open(fin, 'r').read().split('\n')
    types       = (get_types (inlines))
    natoms      = (get_natoms(inlines))
    coords      = (get_coords(inlines))
    cell        = (get_cell  (inlines))
    energy      = (get_energy(outlines))
    force       = (get_force (outlines))
    stress      = (get_stress(outlines))

    ret = write_config(types, natoms, coords, cell, energy, force, stress)
    open(fconfig, 'w').write (ret)

cvt_1frame('input', 'output', 'test.configs')

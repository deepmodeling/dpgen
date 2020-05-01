#!/usr/bin/env python3

import numpy as np
import dpdata
import ase,os


def make_one(out_dir) :
    # [0.5, 1)
    [aa,bb,cc] = np.random.random(3) * 0.5 + 0.5
    # [1, 179)
    [alpha,beta,gamma] = np.random.random(3) * (178 / 180) + 1
    # make cell
    cell = ase.geometry.cellpar_to_cell([aa,bb,cc,alpha,beta,gamma])
    sys = dpdata.System('POSCAR')
    sys['cells'][0] = cell
    os.makedirs(out_dir, exist_ok=True)
    sys.to_vasp_poscar(os.path.join(out_dir, 'POSCAR'))
    

ntest = 30
for ii in range(ntest) :
     out_dir = 'test.%03d' % ii
     make_one(out_dir)

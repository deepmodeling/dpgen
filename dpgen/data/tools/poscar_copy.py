#!/usr/bin/env ovitos

from ovito.io import *
from ovito.modifiers import *

import numpy as np
import argparse

def copy_system (ncopy, fin, fout) :
    nx = ncopy[0]
    ny = ncopy[1]
    nz = ncopy[2]
    node = import_file(fin)
    pbc = ShowPeriodicImagesModifier(adjust_box = True,
                                     num_x = nx, 
                                     num_y = ny,
                                     num_z = nz, 
                                     replicate_x = True,
                                     replicate_y = True,
                                     replicate_z = True
    )
    node.modifiers.append(pbc)
    node.compute()
    export_file(node, fout, 'vasp')

parser = argparse.ArgumentParser(
    description="Copy system")
parser.add_argument('-n', '--ncopy', type=int, nargs = 3,
                    help="the number of copies in each direction")
parser.add_argument('INPUT', type=str,
                    help="the input file")
parser.add_argument('OUTPUT', type=str,
                    help="the output file")
args = parser.parse_args()

copy_system(args.ncopy, args.INPUT, args.OUTPUT)

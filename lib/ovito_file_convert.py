#!/usr/bin/env ovitos
'''
This Script is adapted from Alexander Stukowski, the author of OVITO.
See: http://forum.ovito.org/index.php?topic=131.0 for details.
'''
import os
import sys
import argparse
import numpy as np

from ovito.io import *

supp_ofmt = ['lammps_dump', 'lammps_data', 'vasp']
supp_exts = ['dump', 'lmp', 'poscar/POSCAR']

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--ofmt", type=str, 
                    help="the output format, supported: " + str(supp_ofmt))
parser.add_argument("INPUT", type=str, 
                    help="the input file")
parser.add_argument("OUTPUT", type=str, 
                    help="the output file, supported ext: " + str(supp_exts))
args = parser.parse_args()

fin = args.INPUT
fout = args.OUTPUT
if args.ofmt is not None :
    ofmt = args.ofmt
else :
    ext = fout.split('.')[-1]
    if ext == 'dump' :
        ofmt = 'lammps_dump'
    elif ext == 'lmp' :
        ofmt = 'lammps_data'
    elif ext == 'poscar' or ext == 'POSCAR' :
        ofmt = 'vasp'
if not ofmt in supp_ofmt :
    raise RuntimeError ("output format " + ofmt + " is not supported. use one of " + str(supp_ofmt))

columns = None
if ofmt == "lammps_dump" :
    columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"]

node = import_file(fin)
if columns is not None :
    export_file(node, fout, ofmt, columns = columns)
else :
    export_file(node, fout, ofmt)
    

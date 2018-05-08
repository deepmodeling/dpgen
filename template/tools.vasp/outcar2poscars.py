#!/usr/bin/env python3

import os
import sys
import shutil
import ase.io

# read OUTCAR
foutcar = sys.argv[1]
print("Reading OUTCAR ...")
cfgs = ase.io.read(foutcar, index=':')

dirname = 'POSCARs'
if os.path.exists(dirname):
    shutil.rmtree(dirname)

os.mkdir(dirname)

for idx, cfg in enumerate(cfgs):
    print('Writing image %d ...' % idx)
    fout = 'POSCAR-' + str(idx)
    ase.io.write(fout, cfg, format='vasp', direct=True, vasp5=True)
    shutil.move(fout, dirname)

print('\n All Done.')

#!/usr/bin/env python2

import argparse
import re
import sys
import numpy as np

def get_line_c (fname) :
    myfile = open(fname, "r")
    mylist = []
    for line in myfile :
        if line[0] == '#' :
            mylist.append (line)
    return mylist

def get_key_float (lines, key, nwords = 3) :
    mylist = []
    for line in lines:
        if re.search(key, line):
            words = line.split()[1:1+nwords]
            box = [ float (ii) for ii in words ]
            mylist.append (box)
    return (np.array(mylist))

def _main () :
    parser = argparse.ArgumentParser(
        description="*** Train a model. ***")
    parser.add_argument('INPUT', default = "test.configs",
                        help='the input json database ')
    args = parser.parse_args()

    lc = get_line_c (args.INPUT)

    # get natom
    natom = int(lc[0].split()[1])

    # get box
    boxx = get_key_float (lc, '#X', 3)
    boxy = get_key_float (lc, '#Y', 3)
    boxz = get_key_float (lc, '#Z', 3)
    box = np.concatenate ((boxx, boxy, boxz), axis = 1)
    nframe = box.shape[0]

    # get vol
    vol = np.zeros (nframe)
    for ii in range (nframe) :
        vol[ii] = np.linalg.det (np.reshape(box[ii], [3,3]))

    # get energy
    ener = get_key_float (lc, '#E', 1)
    ener *= natom
    
    # get virial 
    virial = get_key_float (lc, '#S', 6)
    # pref = - 1602 * 1e3 / 1.602176621e6 * 0.5
    pref = 1602 * 1e3 / 1.602176621e6
    for ii in range (nframe) :
        virial[ii] *= pref * vol[ii]
    virial9 = np.zeros ([nframe, 9])
    virial9[:,0] = virial[:,0]
    virial9[:,4] = virial[:,1]
    virial9[:,8] = virial[:,2]
    virial9[:,1] = virial[:,3]
    virial9[:,2] = virial[:,5]
    virial9[:,3] = virial[:,3]
    virial9[:,5] = virial[:,4]
    virial9[:,6] = virial[:,5]
    virial9[:,7] = virial[:,4]

    # get coord, force, type
    data = np.loadtxt (args.INPUT)
    coord = np.reshape (data[:,1:4], [-1, natom * 3])
    force = np.reshape (data[:,4:7], [-1, natom * 3])
    mytype = np.reshape (data[:,0], [-1, natom])
    mytype = mytype[0]

    idx = np.arange (box.shape[0])
    # np.random.shuffle(idx)
    # print idx
    
    np.savetxt ('box.raw', box[idx])
    np.savetxt ('energy.raw', ener[idx])
    np.savetxt ('virial.raw', virial9[idx])
    np.savetxt ('coord.raw', coord[idx])
    np.savetxt ('force.raw', force[idx])
    np.savetxt ('type.raw', mytype, fmt = '%d')

if __name__ == "__main__" :
    _main()

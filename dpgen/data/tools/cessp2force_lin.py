#!/usr/bin/env python2
################################################################
#
# vasp2force:
#   convert vasp output data into potfit reference configurations
#
################################################################
#
#   Copyright 2014
#             Institute for Theoretical and Applied Physics
#             University of Stuttgart, D-70550 Stuttgart, Germany
#             http://potfit.sourceforge.net/
#
#################################################################
#
#   This file is part of potfit.
#
#   potfit is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   potfit is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with potfit; if not, see <http://www.gnu.org/licenses/>.
#
#################################################################

import argparse
import copy
import gzip
import os
import sys


def get_outcar_files(directory, recursive):
    # walk directory (recursively) and return all OUTCAR* files
    # return list of outcars' path
    sys.stderr.write(
        'Searching directory %s for OUTCAR* files ...\n' % directory)
    outcars = []
    if not recursive:
        for item in os.listdir(directory):
            if item.startswith('OUTCAR'):
                outcars.append(os.path.join(directory, item))
    else:
        for root, SubFolders, files in os.walk(directory):
            for item in files:
                if item.startswith('OUTCAR'):
                    outcars.append(os.path.join(root, item))
    if len(outcars) == 0:
        sys.stderr.write(
            'Could not find any OUTCAR files in this directory.\n')
    else:
        sys.stderr.write('Found the following files:\n')
        sys.stderr.write('  {}\n'.format(('\n  ').join(outcars)))
        return outcars
    return outcars


def uniq(seq):
    # return unique list without changing the order
    # from http://stackoverflow.com/questions/480214
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]


def scan_outcar_file(file_handle):
    # scans an outcar file and returns a list with
    # - number of configurations
    # - atom types

    # first try TOTAL-FORCE
    configs = 0
    atom_types = []
    title = []
    potcar = []
    ipt = []
    for line in file_handle:
        if line.startswith('|'):
            continue
        if 'TOTAL-FORCE' in line:
            configs += 1
        if 'VRHFIN' in line:
            atom_types.append(
                line.split()[1].replace('=', '').replace(':', ''))
        if 'title' in line:
            title.append(line.split()[3][0:2])
        if 'POTCAR' in line:
            potcar.append(line.split()[2][0:2])
        if 'ions per type' in line:
            ipt = [int(s) for s in line.split()[4:]]

    potcar = uniq(potcar)

    if atom_types:
        return [configs, atom_types, ipt]
    elif title:
        return [configs, title, ipt]
    elif potcar:
        return [configs, potcar, ipt]
    else:
        sys.stderr.write(
            'Could not determine atom types in file %s.\n' % filename)
        sys.exit()


def process_outcar_file_v5_dev(outcars, data, numbers, types, max_types, elements=None, windex=None, fout='potfit.configs'):
    fw = open(fout, 'w')
    for i in range(len(outcars)):
        if outcars[i].endswith('.gz'):
            f = gzip.open(outcars[i], 'rb')
        else:
            f = open(outcars[i], 'r')

        # Writing current OUTCAR's information into potfit format files.
        nconfs = data[i][0]
        natoms = sum(data[i][2])  # ipt
        if windex == None:
            windex = range(nconfs)
        if windex == 'final':
            windex = [nconfs - 1]

        # reading current OUTCAR
        print("Reading %s ..." % outcars[i])
        count = -1
        line = f.readline()
        while line != '':
            line = f.readline()
            if 'Iteration' in line:
                energy = 0
                box_x = []
                box_y = []
                box_z = []
                stress = []
                atom_data = []
            # if 'energy  without' in line:
            #     # appears in each electronic-iteration steps
            #     energy = float(line.split()[6]) / natoms
            if 'free  energy   TOTEN' in line:
                energy = float (line.split()[4]) / natoms
                if count in windex:
                    fw.write("#N %s 1\n" % natoms)
                    fw.write('#C ')
                    if elements:
                        fw.write("%s " % numbers[0])
                        for j in range(1, max_types):
                            fw.write('%s\t' % numbers[j])
                    else:
                        fw.write(" %s" % data[i][1][0])
                        for j in range(1, max_types):
                            fw.write(' %s' % data[i][1][j])
                    fw.write("\n")
                    fw.write("## force file generated from file %s config %d\n" % (
                        outcars[i], count))
                    fw.write("#X %13.8f %13.8f %13.8f\n" %
                             (box_x[0], box_x[1], box_x[2]))
                    fw.write("#Y %13.8f %13.8f %13.8f\n" %
                             (box_y[0], box_y[1], box_y[2]))
                    fw.write("#Z %13.8f %13.8f %13.8f\n" %
                             (box_z[0], box_z[1], box_z[2]))
                    fw.write("#W %f\n" % (args.weight))
                    fw.write("#E %.10f\n" % (energy))
                    if stress:
                        fw.write("#S ")
                        for num in range(6):
                            fw.write('%8.7g\t' % (stress[num]))
                        fw.write('\n')
                    fw.write("#F\n")
                    fw.flush()
                    for adata in atom_data:
                        fw.write("%d %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g\n" %
                                 (adata[0], adata[1], adata[2], adata[3], adata[4], adata[5], adata[6]))
            if 'VOLUME and BASIS' in line:
                for do in range(5):
                    # SKIP 5 lines
                    line = f.readline()
                box_x = [float(s)
                         for s in line.replace('-', ' -').split()[0:3]]
                line = f.readline()
                box_y = [float(s)
                         for s in line.replace('-', ' -').split()[0:3]]
                line = f.readline()
                box_z = [float(s)
                         for s in line.replace('-', ' -').split()[0:3]]
            if 'in kB' in line:
                stress = [float(s) / 1602 for s in line.split()[2:8]]
            if 'TOTAL-FORCE' in line:
                # only appears in Ion-iteration steps
                line = f.readline()  # skip 1 line
                adata = [0] * 7
                for num in range(len(data[i][2])):
                    for k in range(data[i][2][num]):
                        line = [float(s) for s in f.readline().split()[0:6]]
                        if args.c:
                            adata[0] = int(types[data[i][1][num]])
                        else:
                            adata[0] = int(num)
                        adata[1] = float(line[0])
                        adata[2] = float(line[1])
                        adata[3] = float(line[2])
                        adata[4] = float(line[3])
                        adata[5] = float(line[4])
                        adata[6] = float(line[5])
                        atom_data.append(copy.copy(adata))
                count += 1
    fw.close()
    return


def Parser():
    parser = argparse.ArgumentParser(
        description='''Converts vasp output data into potfit reference configurations.''')

    parser.add_argument('-c', type=str, required=False,
                        help='list of chemical species to use, e.g. -c Mg=0,Zn=1')
    parser.add_argument(
        '-e', type=str, required=False, help='file with single atom energies (NYI)')
    parser.add_argument('-r', '--recursive', action='store_true',
                        help='scan recursively for OUTCAR files')
    parser.add_argument('-f', '--final', action='store_true',
                        help='use only the final configuration from OUTCAR')
    parser.add_argument('-sr', '--configs_range', type=str,
                        help='range of the configurations to use')
    parser.add_argument('-w', '--weight', type=float, default=1.0,
                        help='set configuration weight for all configurations')
    parser.add_argument(
        'files', type=str, nargs='*', help='list of OUTCAR files (plain or gzipped)')

    args = parser.parse_args()
    return args

#
if __name__ == "__main__":
    # Check for sane arguments
    args = Parser()
    if args.weight < 0:
        sys.stderr.write("The weight needs to be positive!\n")
        sys.exit()

    # determine all OUTCAR files
    outcars = []
    if not args.files:
        outcars = get_outcar_files('.', args.recursive)

    for item in args.files:
        if os.path.isdir(item):
            outcars += get_outcar_files(item, args.recursive)
        else:
            outcars.append(os.path.abspath(item))

    # remove duplicate entries from the outcars table
    outcars = uniq(outcars)

    # read metadata from all outcar files
    data = []
    max_types = 1
    for item in outcars:
        if item.endswith('.gz'):
            f = gzip.open(item, 'rb')
        else:
            f = open(item, 'r')
        data.append(scan_outcar_file(f))
        f.close()
        max_types = max(max_types, len(data[-1][1]))

    # check the types string
    types = dict()
    numbers = dict()
    if args.c:
        if len(args.c.split(',')) > max_types:
            sys.stderr.write(
                "\nERROR: There are too many items in you -c string!\n")
            sys.exit()
        if len(args.c.split(',')) < max_types:
            sys.stderr.write(
                "\nERROR: There are not enough items in you -c string!\n")
            sys.exit()
        for item in args.c.split(','):
            if len(item.split('=')) != 2:
                sys.stderr.write("\nERROR: Could not read the -c string.\n")
                sys.stderr.write("Maybe a missing or extra '=' sign?\n")
                sys.exit()
            else:
                try:
                    name = str(item.split('=')[0])
                    number = int(item.split('=')[1])
                except Exception:
                    sys.stderr.write("\nERROR: Could not read the -c string\n")
                    sys.exit()
                if number >= max_types:
                    sys.stderr.write(
                        "\nERROR: The atom type for %s is invalid!\n" % name)
                    sys.exit()
                if name in types:
                    sys.stderr.write(
                        "\nERROR: Duplicate atom type found in -c string\n")
                    sys.exit()
                if number in numbers:
                    sys.stderr.write(
                        "\nERROR: Duplicate atom number found in -c string\n")
                    sys.exit()
                types[name] = number
                numbers[number] = name

    # process all the outcar files
    sr = args.configs_range
    if sr == None:
        windex = None
    else:
        sr = sr.split()
        sr0 = int(sr[0])
        sr1 = int(sr[1])
        windex = range(sr0, sr1 + 1)

    if args.final:
        windex = 'final'

    process_outcar_file_v5_dev(
        outcars, data, numbers, types, max_types, windex=windex, fout='test.configs')

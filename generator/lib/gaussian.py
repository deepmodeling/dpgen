#!/usr/bin/python3


def make_gaussian_input(sys_data, fp_params):
    # keywords, e.g., force b3lyp/6-31g**
    keywords = fp_params['keywords']
    nproc = fp_params['nproc']
    # assume charge is zero and spin multiplicity is 1
    buff = ['%nproc={:d}'.format(nproc), '#force {}'.format(
        keywords), '', 'dpgen', '', '0 1']
    coordinates = sys_data['coords'][0]
    atom_names = (sys_data['atom_names'])
    atom_numbs = (sys_data['atom_numbs'])
    cc = 0
    for ii, atom_name in enumerate(atom_names):
        for _ in range(atom_numbs[ii]):
            buff.append("%s %f %f %f" % (atom_name, *coordinates[cc]))
            cc += 1
    if 'basis_set' in fp_params:
        # custom basis set
        buff.extend(['', fp_params['basis_set'], ''])
    buff.append('\n')
    return '\n'.join(buff)

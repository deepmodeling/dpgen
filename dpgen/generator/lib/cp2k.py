import dpdata
import numpy as np

def make_section(section_name, section_value = None):
    if section_value == None :
        temp_section = '&' + section_name + '\n'
        temp_section += '&END ' + section_name + '\n'
    else :
        temp_section = '&' + section_name + ' ' + section_value + '\n'
        temp_section += '&END ' + section_name + '\n'
    return temp_section

def section_add_subsection(section_string, subsection_string):
    section_string, section_end = section_string.rsplit('\n', 2)[0:2]
    section_string += '\n' + subsection_string + section_end + '\n'
    return section_string

def section_add_keyword_and_value(section_string, keyword, keyword_value):
    section_string, section_end = section_string.rsplit('\n', 2)[0:2]
    section_string += '\n' + keyword + ' ' + keyword_value + '\n' + section_end + '\n'
    return section_string

def make_cp2k_xyz(sys_data):
    #get structral information
    atom_names = sys_data['atom_names']
    atom_types = sys_data['atom_types']

    #write coordinate to xyz file used by cp2k input
    coord_list = sys_data['coords'][0]
    u = np.array(atom_names)
    atom_list = u[atom_types]
    x = '\n'
    for kind, coord in zip(atom_list, coord_list) :
        x += str(kind) + ' ' + str(coord[:])[1:-1] + '\n'
    return x

def make_cp2k_input(sys_data, fp_params):

    #covert cell to cell string
    cell = sys_data['cells'][0]
    cell = np.reshape(cell, [3,3])
    cell_a = np.array2string(cell[0,:])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1,:])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2,:])
    cell_c = cell_c[1:-1]


    #made global section
    global_section = make_section('GLOBAL')
    global_section = section_add_keyword_and_value(global_section, 'PROJECT', 'DPGEN')

    #made force_eval section
    mgrid_section = make_section('MGRID')
    mgrid_section = section_add_keyword_and_value(mgrid_section, 'CUTOFF', fp_params['cutoff'])
    mgrid_section = section_add_keyword_and_value(mgrid_section, 'REL_CUTOFF', fp_params['rel_cutoff'])
    mgrid_section = section_add_keyword_and_value(mgrid_section, 'NGRIDS', '5')

    qs_section = make_section('QS')
    qs_section = section_add_keyword_and_value(qs_section, 'EPS_DEFAULT', '1.0E-12')

    ot_section = make_section('OT')
    ot_section = section_add_keyword_and_value(ot_section, 'MINIMIZER', 'DIIS')
    ot_section = section_add_keyword_and_value(ot_section, 'PRECONDITIONER', 'FULL_SINGLE_INVERSE')

    outer_scf_section = make_section('OUTER_SCF')
    outer_scf_section = section_add_keyword_and_value(outer_scf_section, 'EPS_SCF', '1.0E-6')
    outer_scf_section = section_add_keyword_and_value(outer_scf_section, 'MAX_SCF', '10')

    scf_section = make_section('SCF')
    scf_section = section_add_keyword_and_value(scf_section, 'SCF_GUESS', 'ATOMIC')
    scf_section = section_add_keyword_and_value(scf_section, 'EPS_SCF', '1.0E-6')
    scf_section = section_add_keyword_and_value(scf_section, 'MAX_SCF', '50')
    scf_section = section_add_subsection(scf_section, ot_section)
    scf_section = section_add_subsection(scf_section,outer_scf_section)


    xc_functional_section = make_section('XC_FUNCTIONAL', fp_params['functional'])

    pair_potential_section = make_section('PAIR_POTENTIAL')
    pair_potential_section = section_add_keyword_and_value(pair_potential_section, 'TYPE', 'DFTD3')
    pair_potential_section = section_add_keyword_and_value(pair_potential_section, 'PARAMETER_FILE_NAME', fp_params['pair_potential_path'])
    pair_potential_section = section_add_keyword_and_value(pair_potential_section, 'REFERENCE_FUNCTIONAL', fp_params['pair_ref_functional'])


    vdw_potential_section = make_section('VDW_POTENTIAL')
    vdw_potential_section = section_add_keyword_and_value(vdw_potential_section, 'DISPERSION_FUNCTIONAL', 'PAIR_POTENTIAL')
    vdw_potential_section = section_add_subsection(vdw_potential_section, pair_potential_section)


    xc_section = make_section('XC')
    xc_section = section_add_subsection(xc_section, xc_functional_section)
    if fp_params['pair_potential_path'] !="None":
        xc_section = section_add_subsection(xc_section, vdw_potential_section)


    dft_section = make_section('DFT')
    dft_section = section_add_keyword_and_value(dft_section, 'BASIS_SET_FILE_NAME', fp_params['basis_path'])
    dft_section = section_add_keyword_and_value(dft_section, 'POTENTIAL_FILE_NAME', fp_params['pp_path'])
    dft_section = section_add_keyword_and_value(dft_section, 'CHARGE', '0')
    dft_section = section_add_keyword_and_value(dft_section, 'UKS', 'F')
    dft_section = section_add_keyword_and_value(dft_section, 'MULTIPLICITY', '1')
    dft_section = section_add_subsection(dft_section, mgrid_section)
    dft_section = section_add_subsection(dft_section, qs_section)
    dft_section = section_add_subsection(dft_section, scf_section)
    dft_section = section_add_subsection(dft_section, xc_section)

    cell_section = make_section('CELL')
    cell_section = section_add_keyword_and_value(cell_section, 'A', cell_a)
    cell_section = section_add_keyword_and_value(cell_section, 'B', cell_b)
    cell_section = section_add_keyword_and_value(cell_section, 'C', cell_c)

    coord_section = make_section('COORD')
    coord_section = section_add_keyword_and_value(coord_section, '@include', 'coord.xyz')

    subsys_section = make_section('SUBSYS')
    subsys_section = section_add_subsection(subsys_section, cell_section)
    subsys_section = section_add_subsection(subsys_section, coord_section)

    for kind, basis, potential in zip(fp_params['element_list'], fp_params['basis_list'], fp_params['pp_list']) :
        kind_section = make_section('KIND', kind)
        kind_section = section_add_keyword_and_value(kind_section, 'BASIS_SET', basis)
        kind_section = section_add_keyword_and_value(kind_section, 'POTENTIAL', potential)
        subsys_section = section_add_subsection(subsys_section, kind_section)

    forces_section = make_section('FORCES', 'ON')

    print_section = make_section('PRINT')
    print_section = section_add_subsection(print_section, forces_section)

    force_eval_section = make_section('FORCE_EVAL')
    force_eval_section = section_add_keyword_and_value(force_eval_section, 'METHOD', 'QS')
    force_eval_section = section_add_keyword_and_value(force_eval_section, 'STRESS_TENSOR', 'ANALYTICAL')
    force_eval_section = section_add_subsection(force_eval_section, dft_section)
    force_eval_section = section_add_subsection(force_eval_section, subsys_section)
    force_eval_section = section_add_subsection(force_eval_section, print_section)
    return global_section + force_eval_section





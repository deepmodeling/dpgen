import numpy as np
from dpdata.abacus.scf import get_cell, get_coords
from dpgen.auto_test.lib.vasp import reciprocal_box
import os
bohr2ang = 0.52917721067
def make_abacus_scf_kpt(fp_params):
    # Make KPT file for abacus pw scf calculation.
    # KPT file is the file containing k points infomation in ABACUS scf calculation.
    k_points = [1, 1, 1, 0, 0, 0]
    if "k_points" in fp_params:
        k_points = fp_params["k_points"]
        if len(k_points) != 6:
            raise RuntimeError("k_points has to be a list containig 6 integers specifying MP k points generation.")
    ret = "K_POINTS\n0\nGamma\n"
    for i in range(6):
        ret += str(k_points[i]) + " "
    return ret

def make_abacus_scf_input(fp_params):
    # Make INPUT file for abacus pw scf calculation.
    ret = "INPUT_PARAMETERS\n"
    ret += "calculation scf\n"
    assert(fp_params['ntype'] >= 0 and type(fp_params["ntype"]) == int),  "'ntype' should be a positive integer."
    ret += "ntype %d\n" % fp_params['ntype']
    #ret += "pseudo_dir ./\n"
    if "ecutwfc" in fp_params:
        assert(fp_params["ecutwfc"] >= 0) ,  "'ntype' should be non-negative."
        ret += "ecutwfc %f\n" % fp_params["ecutwfc"]
    if "dr2" in fp_params:
        ret += "dr2 %e\n" % fp_params["dr2"]
    if "niter" in fp_params:
        assert(fp_params['niter'] >= 0 and type(fp_params["niter"])== int), "'niter' should be a positive integer."
        ret += "niter %d\n" % fp_params["niter"]    
    if "basis_type" in fp_params:
        assert(fp_params["basis_type"] in ["pw", "lcao", "lcao_in_pw"]) , "'basis_type' must in 'pw', 'lcao' or 'lcao_in_pw'."
        ret+= "basis_type %s\n" % fp_params["basis_type"]
    if "dft_functional" in fp_params:
        ret += "dft_functional %s\n" % fp_params["dft_functional"]
    if "gamma_only" in fp_params:
        assert(fp_params["gamma_only"] ==1 ) , "'gamma_only' should be 1. Multi-k algorithm will be supported after the KPT generator is completed."
        ret+= "gamma_only %d\n" % fp_params["gamma_only"]  
    if "mixing_type" in fp_params:
        assert(fp_params["mixing_type"] in ["plain", "kerker", "pulay", "pulay-kerker", "broyden"])
        ret += "mixing_type %s\n" % fp_params["mixing_type"]
    if "mixing_beta" in fp_params:
        assert(fp_params["mixing_beta"] >= 0 and fp_params["mixing_beta"] < 1), "'mixing_beta' should between 0 and 1."
        ret += "mixing_beta %f\n" % fp_params["mixing_beta"]
    if "symmetry" in fp_params:
        assert(fp_params["symmetry"] == 0 or fp_params["symmetry"] == 1), "'symmetry' should be either 0 or 1."
        ret += "symmetry %d\n" % fp_params["symmetry"]
    if "nbands" in fp_params:
        assert(fp_params["nbands"] > 0 and type(fp_params["nbands"]) == int), "'nbands' should be a positive integer."
        ret += "nbands %d\n" % fp_params["nbands"]
    if "nspin" in fp_params:
        assert(fp_params["nspin"] == 1 or fp_params["nspin"] == 2 or fp_params["nspin"] == 4), "'nspin' can anly take 1, 2 or 4"
        ret += "nspin %d\n" % fp_params["nspin"]
    if "ks_solver" in fp_params:
        assert(fp_params["ks_solver"] in ["cg", "dav", "lapack", "genelpa", "hpseps", "scalapack_gvx"]), "'ks_sover' should in 'cgx', 'dav', 'lapack', 'genelpa', 'hpseps', 'scalapack_gvx'."
        ret += "ks_solver %s\n" % fp_params["ks_solver"]
    if "smearing" in fp_params:
        assert(fp_params["smearing"] in ["gaussian", "fd", "fixed", "mp", "mp2", "mv"]), "'smearing' should in 'gaussian', 'fd', 'fixed', 'mp', 'mp2', 'mv'. "
        ret += "smearing %s\n" % fp_params["smearing"]
    if "sigma" in fp_params:
        assert(fp_params["sigma"] >= 0), "'sigma' should be non-negative."
        ret += "sigma %f\n" % fp_params["sigma"]
    if "force" in fp_params:
        assert(fp_params["force"] == 0  or fp_params["force"] == 1), "'force' should be either 0 or 1."
        ret += "force %d\n" % fp_params["force"]
    if "stress" in fp_params:
        assert(fp_params["stress"] == 0  or fp_params["stress"] == 1), "'stress' should be either 0 or 1."
        ret += "stress %d\n" % fp_params["stress"]    
    #paras for deepks
    if "out_descriptor" in fp_params:
        assert(fp_params["out_descriptor"] == 0 or fp_params["out_descriptor"] == 1), "'out_descriptor' should be either 0 or 1."
        ret += "out_descriptor %d\n" % fp_params["out_descriptor"]
    if "lmax_descriptor" in fp_params:
        assert(fp_params["lmax_descriptor"] >= 0),  "'lmax_descriptor' should be  a positive integer."
        ret += "lmax_descriptor %d\n" % fp_params["lmax_descriptor"]
    if "deepks_scf" in fp_params:
        assert(fp_params["deepks_scf"] == 0  or fp_params["deepks_scf"] == 1), "'deepks_scf' should be either 0 or 1."
        ret += "deepks_scf %d\n" % fp_params["deepks_scf"]
    if "model_file" in fp_params:
        ret += "model_file %s\n" % fp_params["model_file"]
    return ret

def make_abacus_scf_stru(sys_data, fp_pp_files, fp_params = None):
    atom_names = sys_data['atom_names']
    atom_numbs = sys_data['atom_numbs']
    assert(len(atom_names) == len(fp_pp_files)), "the number of pp_files must be equal to the number of atom types. "
    assert(len(atom_names) == len(atom_numbs)), "Please check the name of atoms. "
    cell = sys_data["cells"].reshape([3, 3])
    coord = sys_data['coords'].reshape([sum(atom_numbs), 3])
    #volume = np.linalg.det(cell)
    #lattice_const = np.power(volume, 1/3)

    ret = "ATOMIC_SPECIES\n"
    for iatom in range(len(atom_names)):
        if 'atom_masses' not in sys_data:
            ret += atom_names[iatom] + " 1.00 " + fp_pp_files[iatom] + "\n"
        else:
            ret += atom_names[iatom] + " %.3f "%sys_data['atom_masses'][iatom] + fp_pp_files[iatom] + "\n"

    if fp_params is not None and "lattice_constant" in fp_params:
        ret += "\nLATTICE_CONSTANT\n"
        ret += str(fp_params["lattice_constant"]) + "\n\n" # in Bohr, in this way coord and cell are in Angstrom 
    else:
        ret += "\nLATTICE_CONSTANT\n"
        ret += str(1/bohr2ang) + "\n\n"

    ret += "LATTICE_VECTORS\n"
    for ix in range(3):
        for iy in range(3):
            ret += str(cell[ix][iy]) + " "
        ret += "\n"
    ret += "\n"

    ret += "ATOMIC_POSITIONS\n"
    ret += "Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)\n"
    #ret += "\n"
    natom_tot = 0
    for iele in range(len(atom_names)):
        ret += atom_names[iele] + "\n"
        ret += "0.0\n"
        ret += str(atom_numbs[iele]) + "\n"
        for iatom in range(atom_numbs[iele]):
            ret += "%.12f %.12f %.12f %d %d %d\n" % (coord[natom_tot, 0], coord[natom_tot, 1], coord[natom_tot, 2], 1, 1, 1)
            natom_tot += 1
    assert(natom_tot == sum(atom_numbs))

    if fp_params is not None and "basis_type" in fp_params and fp_params["basis_type"]=="lcao":
        ret +="\nNUMERICAL_ORBITAL\n"
        assert(len(fp_params["orb_files"])==len(atom_names))
        for iatom in range(len(atom_names)):
            ret += fp_params["orb_files"][iatom] +"\n"

    if fp_params is not None and "deepks_scf" in fp_params and fp_params["out_descriptor"]==1:
        ret +="\nNUMERICAL_DESCRIPTOR\n"
        ret +=fp_params["proj_file"][0]+"\n"

    return ret

def get_abacus_input_parameters(INPUT):
    with open(INPUT) as fp:
        inlines = fp.read().split("\n")
    input_parameters = {}
    for line in inlines:
        if line.split() == [] or len(line.split()) < 2 :
            continue
        parameter_name = line.split()[0]
        parameter_value = line.split()[1]
        input_parameters[parameter_name] = parameter_value
    return input_parameters

def get_mass_from_STRU(geometry_inlines, inlines, atom_names):
    nele = None
    for line in inlines:
        if line.split() == []:
            continue
        if "ntype" in line and "ntype" == line.split()[0]:
            nele = int(line.split()[1])
    assert(nele is not None)
    mass_list = [0 for i in atom_names]
    pp_file_list = [i for i in atom_names]
    for iline, line in enumerate(geometry_inlines):
        if line.split() == []:
            continue
        if "ATOMIC_SPECIES" == line.split()[0]:
            for iele1 in range(1, 1+nele):
                for iele2 in range(nele):
                    if geometry_inlines[iline+iele1].split()[0] == atom_names[iele2]:
                        mass_list[iele2] = float(geometry_inlines[iline+iele1].split()[1])
                        pp_file_list[iele2] = geometry_inlines[iline+iele1].split()[2]
    for iele in range(len(mass_list)):
        assert(mass_list[iele] > 0)
    return mass_list, pp_file_list

def get_natoms_from_stru(geometry_inlines):
    key_words_list = ["ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS"]
    keyword_sequence = []
    keyword_line_index = []
    atom_names = []
    atom_numbs = []
    for iline, line in enumerate(geometry_inlines):
        if line.split() == []:
            continue
        have_key_word = False
        for keyword in key_words_list:
            if keyword in line and keyword == line.split()[0]:
                keyword_sequence.append(keyword)
                keyword_line_index.append(iline)
    assert(len(keyword_line_index) == len(keyword_sequence))
    assert(len(keyword_sequence) > 0)
    keyword_line_index.append(len(geometry_inlines))
    for idx, keyword in enumerate(keyword_sequence):
        if keyword == "ATOMIC_POSITIONS":
            iline = keyword_line_index[idx]+2
            while iline < keyword_line_index[idx+1]-1:
                atom_names.append(geometry_inlines[iline].split()[0])
                atom_numbs.append(int(geometry_inlines[iline+2].split()[0]))
                iline += 3+atom_numbs[-1]
    return atom_names, atom_numbs
def get_abacus_STRU(STRU, INPUT = None, n_ele = None):
    # read in geometry from STRU file. n_ele is the number of elements.
    # Either n_ele or INPUT should be provided.
    with open(STRU, 'r') as fp:
        geometry_inlines = fp.read().split('\n')
    for iline, line in enumerate(geometry_inlines):
        if line.split() == []:
            del geometry_inlines[iline]
    geometry_inlines.append("")
    celldm, cell = get_cell(geometry_inlines) 
    if n_ele is None and INPUT is not None:
        assert(os.path.isfile(INPUT)), "file %s should exists" % INPUT
        with open(INPUT, 'r') as fp:
            inlines = fp.read().split('\n')
        atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines, inlines) 
    elif n_ele is not None and INPUT is None:
        assert(n_ele > 0)
        inlines = ["ntype %d" %n_ele]
        atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines, inlines)
    else:
        atom_names, atom_numbs = get_natoms_from_stru(geometry_inlines)
        inlines = ["ntype %d" %len(atom_numbs)]
        atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines, inlines)
    masses, pp_files = get_mass_from_STRU(geometry_inlines, inlines, atom_names)
    data = {}
    data['atom_names'] = atom_names
    data['atom_numbs'] = natoms
    data['atom_types'] = types
    data['cells'] = cell
    data['coords'] = coords
    data['atom_masses'] = masses # Notice that this key is not defined in dpdata system. 
    data['pp_files'] = pp_files
    return data

def make_supercell_abacus(from_struct, super_cell):
    if "types" in from_struct:
        from_struct["types"] = from_struct["types"] * super_cell[0] * super_cell[1] * super_cell[2]
    for ix in range(super_cell[0]):
        for iy in range(super_cell[1]):
            for iz in range(super_cell[2]):
                if ix == 0 and iy == 0 and iz == 0:
                    continue
                for ia in range(sum(from_struct["atom_numbs"])):
                    coord = from_struct['coords'][ia] + from_struct['cells'][0]*ix + from_struct['cells'][1]*iy + from_struct['cells'][2]*iz
                    from_struct['coords'] = np.vstack([from_struct['coords'], coord])
    from_struct["atom_numbs"] = [i * super_cell[0] * super_cell[1] * super_cell[2] for i in from_struct["atom_numbs"]]
    from_struct['cells'][0] *= super_cell[0]
    from_struct['cells'][1] *= super_cell[1]
    from_struct['cells'][2] *= super_cell[2]
    return from_struct

def make_kspacing_kpoints_stru(stru, kspacing) :
    # adapted from dpgen.autotest.lib.vasp.make_kspacing_kpoints
    if type(kspacing) is not list:
        kspacing = [kspacing, kspacing, kspacing]
    box = stru['cells']
    rbox = reciprocal_box(box)
    kpoints = [max(1,(np.ceil(2 * np.pi * np.linalg.norm(ii) / ks).astype(int))) for ii,ks in zip(rbox,kspacing)]
    kpoints += [0, 0, 0]
    return kpoints

if __name__ == "__main__":
    fp_params = {"k_points": [1, 1, 1, 0, 0, 0]}
    ret = make_abacus_scf_kpt(fp_params)
    print(ret)
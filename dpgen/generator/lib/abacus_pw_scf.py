import numpy as np
bohr2ang = 0.52917721067
def make_abacus_pw_scf_kpt(fp_params):
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

def make_abacus_pw_scf_input(fp_params):
    # Make INPUT file for abacus pw scf calculation.
    ret = "INPUT_PARAMETERS\n"
    assert(fp_params['ntype'] >= 0 and type(fp_params["ntype"]) == int)
    ret += "ntype %d\n" % fp_params['ntype']
    ret += "pseudo_dir ./\n"
    if "ecutwfc" in fp_params:
        assert(fp_params["ecutwfc"] >= 0)
        ret += "ecutwfc %f\n" % fp_params["ecutwfc"]
    if "mixing_type" in fp_params:
        assert(fp_params["mixing_type"] in ["plain", "kerker", "pulay", "pulay-kerker", "broyden"])
        ret += "mixing_type %s\n" % fp_params["mixing_type"]
    if "mixing_beta" in fp_params:
        assert(fp_params["mixing_beta"] >= 0 and fp_params["mixing_beta"] < 1)
        ret += "mixing_beta %f\n" % fp_params["mixing_beta"]
    if "symmetry" in fp_params:
        assert(fp_params["symmetry"] == 0 or fp_params["symmetry"] == 1)
        ret += "symmetry %d\n" % fp_params["symmetry"]
    if "nbands" in fp_params:
        assert(fp_params["nbands"] > 0 and type(fp_params["nbands"]) == int)
        ret += "nbands %f\n" % fp_params["nbands"]
    if "nspin" in fp_params:
        assert(fp_params["nspin"] == 1 or fp_params["nspin"] == 2 or fp_params["nspin"] == 4)
        ret += "nspin %d\n" % fp_params["nspin"]
    if "ks_solver" in fp_params:
        assert(fp_params["ks_solver"] in ["cg", "dav", "lapack", "genelpa", "hpseps", "scalapack_gvx"])
        ret += "ks_solver %s\n" % fp_params["ks_solver"]
    if "smearing" in fp_params:
        assert(fp_params["smearing"] in ["gauss", "fd", "fixed", "mp", "mp2", "mv"])
        ret += "smearing %s\n" % fp_params["smearing"]
    if "sigma" in fp_params:
        assert(fp_params["sigma"] >= 0)
        ret += "sigma %f\n" % fp_params["sigma"]
    ret += "force 1\nstress 1\n"
    return ret

def make_abacus_pw_scf_stru(sys_data, fp_pp_files):
    atom_names = sys_data['atom_names']
    atom_numbs = sys_data['atom_numbs']
    assert(len(atom_names) == len(fp_pp_files))
    assert(len(atom_names) == len(atom_numbs))
    cell = sys_data["cells"][0].reshape([3, 3])
    coord = sys_data['coords'][0]
    #volume = np.linalg.det(cell)
    #lattice_const = np.power(volume, 1/3)
    lattice_const = 1/bohr2ang # in Bohr, in this way coord and cell are in Angstrom 

    ret = "ATOMIC_SPECIES\n"
    for iatom in range(len(atom_names)):
        ret += atom_names[iatom] + " 1.00 " + fp_pp_files[iatom] + "\n"
    
    ret += "\nLATTICE_CONSTANT\n"
    ret += str(lattice_const) + "\n\n"

    ret += "LATTICE_VECTORS\n"
    for ix in range(3):
        for iy in range(3):
            ret += str(cell[ix][iy]) + " "
        ret += "\n"
    ret += "\n"

    ret += "ATOMIC_POSITIONS\n"
    ret += "Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)\n"
    natom_tot = 0
    for iele in range(len(atom_names)):
        ret += atom_names[iele] + "\n"
        ret += "0.0\n"
        ret += str(atom_numbs[iele]) + "\n"
        for iatom in range(atom_numbs[iele]):
            ret += "%.12f %.12f %.12f %d %d %d\n" % (coord[natom_tot, 0], coord[natom_tot, 1], coord[natom_tot, 2], 0, 0, 0)
            natom_tot += 1
    assert(natom_tot == sum(atom_numbs))

    return ret


if __name__ == "__main__":
    fp_params = {"k_points": [1, 1, 1, 0, 0, 0]}
    ret = make_abacus_pw_scf_kpt(fp_params)
    print(ret)
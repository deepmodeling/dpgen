import copy
import os
import re

import numpy as np
from dpdata.abacus.scf import get_cell, get_coords, get_nele_from_stru

from dpgen.auto_test.lib import vasp

bohr2ang = 0.52917721067


def make_abacus_scf_kpt(fp_params):
    # Make KPT file for abacus pw scf calculation.
    # KPT file is the file containing k points infomation in ABACUS scf calculation.
    k_points = [1, 1, 1, 0, 0, 0]
    if "k_points" in fp_params:
        k_points = fp_params["k_points"]
        if len(k_points) != 6:
            raise RuntimeError(
                "k_points has to be a list containig 6 integers specifying MP k points generation."
            )
    ret = "K_POINTS\n0\nGamma\n"
    for i in range(6):
        ret += str(k_points[i]) + " "
    return ret


def make_abacus_scf_input(fp_params):
    # Make INPUT file for abacus pw scf calculation.
    ret = "INPUT_PARAMETERS\n"
    ret += "calculation scf\n"
    for key in fp_params:
        if key == "ecutwfc":
            fp_params["ecutwfc"] = float(fp_params["ecutwfc"])
            assert fp_params["ecutwfc"] >= 0, "'ecutwfc' should be non-negative."
            ret += "ecutwfc %f\n" % fp_params["ecutwfc"]
        elif key == "kspacing":
            if isinstance(fp_params["kspacing"], (int, float)):
                fp_params["kspacing"] = [float(fp_params["kspacing"])]
            elif isinstance(fp_params["kspacing"], (list, tuple)):
                fp_params["kspacing"] = list(fp_params["kspacing"])
            elif isinstance(fp_params["kspacing"], str):
                fp_params["kspacing"] = [
                    float(i) for i in fp_params["kspacing"].split()
                ]
            assert len(fp_params["kspacing"]) in [
                1,
                3,
            ], "'kspacing' only accept a float, or a list of one or three float, or a string of one or three float"
            ret += "kspacing "
            for ikspacing in fp_params["kspacing"]:
                assert ikspacing >= 0, "'kspacing' should be non-negative."
                ret += "%f " % ikspacing
            ret += "\n"
        elif key == "scf_thr":
            fp_params["scf_thr"] = float(fp_params["scf_thr"])
            ret += "scf_thr %e\n" % fp_params["scf_thr"]
        elif key == "scf_nmax":
            fp_params["scf_nmax"] = int(fp_params["scf_nmax"])
            assert (
                fp_params["scf_nmax"] >= 0 and type(fp_params["scf_nmax"]) == int
            ), "'scf_nmax' should be a positive integer."
            ret += "scf_nmax %d\n" % fp_params["scf_nmax"]
        elif key == "basis_type":
            assert fp_params["basis_type"] in [
                "pw",
                "lcao",
                "lcao_in_pw",
            ], "'basis_type' must in 'pw', 'lcao' or 'lcao_in_pw'."
            ret += "basis_type %s\n" % fp_params["basis_type"]
        elif key == "dft_functional":
            ret += "dft_functional %s\n" % fp_params["dft_functional"]
        elif key == "gamma_only":
            if type(fp_params["gamma_only"]) == str:
                fp_params["gamma_only"] = int(eval(fp_params["gamma_only"]))
            assert (
                fp_params["gamma_only"] == 0 or fp_params["gamma_only"] == 1
            ), "'gamma_only' should be either 0 or 1."
            ret += "gamma_only %d\n" % fp_params["gamma_only"]
        elif key == "mixing_type":
            assert fp_params["mixing_type"] in [
                "plain",
                "kerker",
                "pulay",
                "pulay-kerker",
                "broyden",
            ]
            ret += "mixing_type %s\n" % fp_params["mixing_type"]
        elif key == "mixing_beta":
            fp_params["mixing_beta"] = float(fp_params["mixing_beta"])
            assert (
                fp_params["mixing_beta"] >= 0 and fp_params["mixing_beta"] < 1
            ), "'mixing_beta' should between 0 and 1."
            ret += "mixing_beta %f\n" % fp_params["mixing_beta"]
        elif key == "symmetry":
            if type(fp_params["symmetry"]) == str:
                fp_params["symmetry"] = int(eval(fp_params["symmetry"]))
            assert (
                fp_params["symmetry"] == 0 or fp_params["symmetry"] == 1
            ), "'symmetry' should be either 0 or 1."
            ret += "symmetry %d\n" % fp_params["symmetry"]
        elif key == "nbands":
            fp_params["nbands"] = int(fp_params["nbands"])
            assert (
                fp_params["nbands"] > 0 and type(fp_params["nbands"]) == int
            ), "'nbands' should be a positive integer."
            ret += "nbands %d\n" % fp_params["nbands"]
        elif key == "nspin":
            fp_params["nspin"] = int(fp_params["nspin"])
            assert (
                fp_params["nspin"] == 1
                or fp_params["nspin"] == 2
                or fp_params["nspin"] == 4
            ), "'nspin' can anly take 1, 2 or 4"
            ret += "nspin %d\n" % fp_params["nspin"]
        elif key == "ks_solver":
            assert fp_params["ks_solver"] in [
                "cg",
                "dav",
                "lapack",
                "genelpa",
                "hpseps",
                "scalapack_gvx",
            ], "'ks_sover' should in 'cgx', 'dav', 'lapack', 'genelpa', 'hpseps', 'scalapack_gvx'."
            ret += "ks_solver %s\n" % fp_params["ks_solver"]
        elif key == "smearing_method":
            assert fp_params["smearing_method"] in [
                "gauss",
                "gaussian",
                "fd",
                "fixed",
                "mp",
                "mp2",
                "mv",
            ], "'smearing_method' should in 'gauss', 'gaussian', 'fd', 'fixed', 'mp', 'mp2', 'mv'. "
            ret += "smearing_method %s\n" % fp_params["smearing_method"]
        elif key == "smearing_sigma":
            fp_params["smearing_sigma"] = float(fp_params["smearing_sigma"])
            assert (
                fp_params["smearing_sigma"] >= 0
            ), "'smearing_sigma' should be non-negative."
            ret += "smearing_sigma %f\n" % fp_params["smearing_sigma"]
        elif key == "cal_force":
            if type(fp_params["cal_force"]) == str:
                fp_params["cal_force"] = int(eval(fp_params["cal_force"]))
            assert (
                fp_params["cal_force"] == 0 or fp_params["cal_force"] == 1
            ), "'cal_force' should be either 0 or 1."
            ret += "cal_force %d\n" % fp_params["cal_force"]
        elif key == "cal_stress":
            if type(fp_params["cal_stress"]) == str:
                fp_params["cal_stress"] = int(eval(fp_params["cal_stress"]))
            assert (
                fp_params["cal_stress"] == 0 or fp_params["cal_stress"] == 1
            ), "'cal_stress' should be either 0 or 1."
            ret += "cal_stress %d\n" % fp_params["cal_stress"]
        # paras for deepks
        elif key == "deepks_out_labels":
            if type(fp_params["deepks_out_labels"]) == str:
                fp_params["deepks_out_labels"] = int(
                    eval(fp_params["deepks_out_labels"])
                )
            assert (
                fp_params["deepks_out_labels"] == 0
                or fp_params["deepks_out_labels"] == 1
            ), "'deepks_out_labels' should be either 0 or 1."
            ret += "deepks_out_labels %d\n" % fp_params["deepks_out_labels"]
        elif key == "deepks_descriptor_lmax":
            fp_params["deepks_descriptor_lmax"] = int(
                fp_params["deepks_descriptor_lmax"]
            )
            assert (
                fp_params["deepks_descriptor_lmax"] >= 0
            ), "'deepks_descriptor_lmax' should be  a positive integer."
            ret += "deepks_descriptor_lmax %d\n" % fp_params["deepks_descriptor_lmax"]
        elif key == "deepks_scf":
            if type(fp_params["deepks_scf"]) == str:
                fp_params["deepks_scf"] = int(eval(fp_params["deepks_scf"]))
            assert (
                fp_params["deepks_scf"] == 0 or fp_params["deepks_scf"] == 1
            ), "'deepks_scf' should be either 0 or 1."
            ret += "deepks_scf %d\n" % fp_params["deepks_scf"]
        elif key == "deepks_model":
            ret += "deepks_model %s\n" % fp_params["deepks_model"]
        elif key[0] == "_":
            pass
        elif key == "calculation":
            pass
        else:
            ret += "%s %s\n" % (key, str(fp_params[key]))
    return ret


def make_abacus_scf_stru(
    sys_data,
    fp_pp_files,
    fp_orb_files=None,
    fp_dpks_descriptor=None,
    fp_params=None,
    type_map=None,
):
    atom_names = sys_data["atom_names"]
    atom_numbs = sys_data["atom_numbs"]
    if type_map == None:
        type_map = atom_names

    assert len(atom_names) == len(atom_numbs), "Please check the name of atoms. "
    cell = sys_data["cells"].reshape([3, 3])
    coord = sys_data["coords"].reshape([sum(atom_numbs), 3])
    # volume = np.linalg.det(cell)
    # lattice_const = np.power(volume, 1/3)

    ret = "ATOMIC_SPECIES\n"
    for iatom in range(len(atom_names)):
        assert atom_names[iatom] in type_map, (
            "element %s is not defined in type_map" % atom_names[iatom]
        )
        idx = type_map.index(atom_names[iatom])
        if "atom_masses" not in sys_data:
            ret += atom_names[iatom] + " 1.00 " + fp_pp_files[idx] + "\n"
        else:
            ret += (
                atom_names[iatom]
                + " %.3f " % sys_data["atom_masses"][iatom]
                + fp_pp_files[idx]
                + "\n"
            )

    if fp_params is not None and "lattice_constant" in fp_params:
        ret += "\nLATTICE_CONSTANT\n"
        ret += (
            str(fp_params["lattice_constant"]) + "\n\n"
        )  # in Bohr, in this way coord and cell are in Angstrom
    else:
        ret += "\nLATTICE_CONSTANT\n"
        ret += str(1 / bohr2ang) + "\n\n"

    ret += "LATTICE_VECTORS\n"
    for ix in range(3):
        for iy in range(3):
            ret += str(cell[ix][iy]) + " "
        ret += "\n"
    ret += "\n"

    ret += "ATOMIC_POSITIONS\n"
    ret += "Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)\n"
    # ret += "\n"
    natom_tot = 0
    for iele in range(len(atom_names)):
        ret += atom_names[iele] + "\n"
        ret += "0.0\n"
        ret += str(atom_numbs[iele]) + "\n"
        for iatom in range(atom_numbs[iele]):
            ret += "%.12f %.12f %.12f %d %d %d\n" % (
                coord[natom_tot, 0],
                coord[natom_tot, 1],
                coord[natom_tot, 2],
                1,
                1,
                1,
            )
            natom_tot += 1
    assert natom_tot == sum(atom_numbs)

    if fp_orb_files is not None:
        ret += "\nNUMERICAL_ORBITAL\n"
        assert len(fp_orb_files) == len(type_map)
        for iatom in range(len(atom_names)):
            idx = type_map.index(atom_names[iatom])
            ret += fp_orb_files[idx] + "\n"

    if fp_dpks_descriptor is not None:
        ret += "\nNUMERICAL_DESCRIPTOR\n"
        ret += "%s\n" % fp_dpks_descriptor

    return ret


def get_abacus_input_parameters(INPUT):
    with open(INPUT) as fp:
        inlines = fp.read().split("\n")
    input_parameters = {}
    for line in inlines:
        sline = re.split("[ \t]", line.split("#")[0].strip(), maxsplit=1)
        if len(sline) == 2:
            input_parameters[sline[0].strip()] = sline[1].strip()
    fp.close()
    return input_parameters


def get_mass_from_STRU(geometry_inlines, atom_names):
    nele = get_nele_from_stru(geometry_inlines)
    assert nele > 0
    mass_list = [0 for i in atom_names]
    pp_file_list = [i for i in atom_names]
    for iline, line in enumerate(geometry_inlines):
        if line.split() == []:
            continue
        if "ATOMIC_SPECIES" == line.split()[0]:
            for iele1 in range(1, 1 + nele):
                for iele2 in range(nele):
                    if geometry_inlines[iline + iele1].split()[0] == atom_names[iele2]:
                        mass_list[iele2] = float(
                            geometry_inlines[iline + iele1].split()[1]
                        )
                        pp_file_list[iele2] = geometry_inlines[iline + iele1].split()[2]
    for iele in range(len(mass_list)):
        assert mass_list[iele] > 0
    return mass_list, pp_file_list


def get_natoms_from_stru(geometry_inlines):
    key_words_list = [
        "ATOMIC_SPECIES",
        "NUMERICAL_ORBITAL",
        "LATTICE_CONSTANT",
        "LATTICE_VECTORS",
        "ATOMIC_POSITIONS",
        "NUMERICAL_DESCRIPTOR",
    ]
    keyword_sequence = []
    keyword_line_index = []
    atom_names = []
    atom_numbs = []
    tmp_line = []
    for i in geometry_inlines:
        if i.strip() != "":
            tmp_line.append(i)
    for iline, line in enumerate(tmp_line):
        if line.split() == []:
            continue
        have_key_word = False
        for keyword in key_words_list:
            if keyword in line and keyword == line.split()[0]:
                keyword_sequence.append(keyword)
                keyword_line_index.append(iline)
    assert len(keyword_line_index) == len(keyword_sequence)
    assert len(keyword_sequence) > 0
    keyword_line_index.append(len(tmp_line))
    for idx, keyword in enumerate(keyword_sequence):
        if keyword == "ATOMIC_POSITIONS":
            iline = keyword_line_index[idx] + 2
            while iline < keyword_line_index[idx + 1] - 1:
                atom_names.append(tmp_line[iline].split()[0])
                atom_numbs.append(int(tmp_line[iline + 2].split()[0]))
                iline += 3 + atom_numbs[-1]
    return atom_names, atom_numbs


def get_additional_from_STRU(geometry_inlines, nele):
    dpks_descriptor_kw = "NUMERICAL_DESCRIPTOR"
    orb_file_kw = "NUMERICAL_ORBITAL"
    dpks_descriptor = None
    orb_file = None
    for iline in range(len(geometry_inlines)):
        if len(geometry_inlines[iline]) > 0:
            if orb_file_kw == geometry_inlines[iline].split()[0]:
                orb_file = []
                for iele in range(nele):
                    orb_file.append(geometry_inlines[iline + iele + 1].strip())
            if dpks_descriptor_kw == geometry_inlines[iline].split()[0]:
                dpks_descriptor = geometry_inlines[iline + 1].strip()
    return orb_file, dpks_descriptor


def get_abacus_STRU(STRU, INPUT=None, n_ele=None):
    # read in geometry from STRU file. n_ele is the number of elements.
    # Either n_ele or INPUT should be provided.
    with open(STRU, "r") as fp:
        geometry_inlines = fp.read().split("\n")
    for iline, line in enumerate(geometry_inlines):
        if line.split() == [] or len(line) == 0:
            del geometry_inlines[iline]
    geometry_inlines.append("")
    celldm, cell = get_cell(geometry_inlines)
    atom_names, natoms, types, coords = get_coords(celldm, cell, geometry_inlines)
    masses, pp_files = get_mass_from_STRU(geometry_inlines, atom_names)
    orb_files, dpks_descriptor = get_additional_from_STRU(geometry_inlines, len(masses))
    data = {}
    data["atom_names"] = atom_names
    data["atom_numbs"] = natoms
    data["atom_types"] = types
    data["cells"] = cell
    data["coords"] = coords
    data[
        "atom_masses"
    ] = masses  # Notice that this key is not defined in dpdata system.
    data["pp_files"] = pp_files
    data["orb_files"] = orb_files
    data["dpks_descriptor"] = dpks_descriptor
    return data


def make_supercell_abacus(from_struct, super_cell):
    to_struct = copy.deepcopy(from_struct)

    if "atom_types" in from_struct:
        new_types = []
        # to_struct["atom_types"] = (
        #    from_struct["atom_types"] * super_cell[0] * super_cell[1] * super_cell[2]
        # )
        for idx_atm, ina in enumerate(from_struct["atom_numbs"]):
            new_types += (
                [idx_atm for i in range(ina)]
                * super_cell[0]
                * super_cell[1]
                * super_cell[2]
            )
        to_struct["atom_types"] = new_types
    new_coord = None
    for ia in range(sum(from_struct["atom_numbs"])):
        for ix in range(super_cell[0]):
            for iy in range(super_cell[1]):
                for iz in range(super_cell[2]):
                    # if ix == 0 and iy == 0 and iz == 0:
                    #    continue

                    coord = (
                        from_struct["coords"][ia]
                        + from_struct["cells"][0] * ix
                        + from_struct["cells"][1] * iy
                        + from_struct["cells"][2] * iz
                    )
                    if ia == 0 and ix == 0 and iy == 0 and iz == 0:
                        new_coord = np.array(coord)
                    else:
                        new_coord = np.vstack([new_coord, coord])
    to_struct["coords"] = new_coord
    new_numbs = [
        i * super_cell[0] * super_cell[1] * super_cell[2]
        for i in from_struct["atom_numbs"]
    ]
    to_struct["atom_numbs"] = new_numbs
    to_struct["cells"][0] *= super_cell[0]
    to_struct["cells"][1] *= super_cell[1]
    to_struct["cells"][2] *= super_cell[2]
    return to_struct


def make_kspacing_kpoints_stru(stru, kspacing):
    # adapted from dpgen.autotest.lib.vasp.make_kspacing_kpoints
    if type(kspacing) is not list:
        kspacing = [kspacing, kspacing, kspacing]
    box = stru["cells"]
    rbox = vasp.reciprocal_box(box)
    kpoints = [
        max(1, (np.ceil(2 * np.pi * np.linalg.norm(ii) / ks).astype(int)))
        for ii, ks in zip(rbox, kspacing)
    ]
    kpoints += [0, 0, 0]
    return kpoints


if __name__ == "__main__":
    fp_params = {"k_points": [1, 1, 1, 0, 0, 0]}
    ret = make_abacus_scf_kpt(fp_params)
    print(ret)

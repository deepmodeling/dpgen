import copy
import os
import re

import numpy as np

from dpdata.abacus.stru import get_frame_from_stru, make_unlabeled_stru

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


def make_abacus_scf_input(fp_params, extra_file_path=""):
    # Make INPUT file for abacus pw scf calculation.
    # put extra files (such as: deepks_model) to extra_file_path folder
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
            assert (
                len(fp_params["kspacing"])
                in [
                    1,
                    3,
                ]
            ), "'kspacing' only accept a float, or a list of one or three float, or a string of one or three float"
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
            assert fp_params["scf_nmax"] >= 0 and isinstance(
                fp_params["scf_nmax"], int
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
            if isinstance(fp_params["gamma_only"], str):
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
            if isinstance(fp_params["symmetry"], str):
                fp_params["symmetry"] = int(eval(fp_params["symmetry"]))
            assert (
                fp_params["symmetry"] == 0 or fp_params["symmetry"] == 1
            ), "'symmetry' should be either 0 or 1."
            ret += "symmetry %d\n" % fp_params["symmetry"]
        elif key == "nbands":
            fp_params["nbands"] = int(fp_params["nbands"])
            assert fp_params["nbands"] > 0 and isinstance(
                fp_params["nbands"], int
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
            assert (
                fp_params["ks_solver"]
                in [
                    "cg",
                    "dav",
                    "lapack",
                    "genelpa",
                    "hpseps",
                    "scalapack_gvx",
                ]
            ), "'ks_sover' should in 'cgx', 'dav', 'lapack', 'genelpa', 'hpseps', 'scalapack_gvx'."
            ret += "ks_solver %s\n" % fp_params["ks_solver"]
        elif key == "smearing_method":
            assert (
                fp_params["smearing_method"]
                in [
                    "gauss",
                    "gaussian",
                    "fd",
                    "fixed",
                    "mp",
                    "mp2",
                    "mv",
                ]
            ), "'smearing_method' should in 'gauss', 'gaussian', 'fd', 'fixed', 'mp', 'mp2', 'mv'. "
            ret += "smearing_method %s\n" % fp_params["smearing_method"]
        elif key == "smearing_sigma":
            fp_params["smearing_sigma"] = float(fp_params["smearing_sigma"])
            assert (
                fp_params["smearing_sigma"] >= 0
            ), "'smearing_sigma' should be non-negative."
            ret += "smearing_sigma %f\n" % fp_params["smearing_sigma"]
        elif key == "cal_force":
            if isinstance(fp_params["cal_force"], str):
                fp_params["cal_force"] = int(eval(fp_params["cal_force"]))
            assert (
                fp_params["cal_force"] == 0 or fp_params["cal_force"] == 1
            ), "'cal_force' should be either 0 or 1."
            ret += "cal_force %d\n" % fp_params["cal_force"]
        elif key == "cal_stress":
            if isinstance(fp_params["cal_stress"], str):
                fp_params["cal_stress"] = int(eval(fp_params["cal_stress"]))
            assert (
                fp_params["cal_stress"] == 0 or fp_params["cal_stress"] == 1
            ), "'cal_stress' should be either 0 or 1."
            ret += "cal_stress %d\n" % fp_params["cal_stress"]
        # paras for deepks
        elif key == "deepks_out_labels":
            if isinstance(fp_params["deepks_out_labels"], str):
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
            if isinstance(fp_params["deepks_scf"], str):
                fp_params["deepks_scf"] = int(eval(fp_params["deepks_scf"]))
            assert (
                fp_params["deepks_scf"] == 0 or fp_params["deepks_scf"] == 1
            ), "'deepks_scf' should be either 0 or 1."
            ret += "deepks_scf %d\n" % fp_params["deepks_scf"]
        elif key == "deepks_model":
            ret += "deepks_model %s\n" % os.path.join(
                extra_file_path, os.path.split(fp_params["deepks_model"])[1]
            )
        elif key[0] == "_":
            pass
        elif key == "calculation":
            pass
        else:
            ret += f"{key} {str(fp_params[key])}\n"
    return ret


def make_abacus_scf_stru(
    sys_data,
    fp_pp_files,
    fp_orb_files=None,
    fp_dpks_descriptor=None,
    type_map=None,
    pporb="",  # pull all pp orb dpks files to pporb folder
):
    # re-construct the path of files by pporb + file name
    fp_pp_files = [os.path.join(pporb, i) for i in fp_pp_files]
    if fp_orb_files is not None:
        fp_orb_files = [os.path.join(pporb, i) for i in fp_orb_files]
    if fp_dpks_descriptor is not None:
        fp_dpks_descriptor = os.path.join(pporb, fp_dpks_descriptor)
    
    c = make_unlabeled_stru(sys_data, 0, pp_files=fp_pp_files, numerical_orbital=fp_orb_files, numerical_descriptor=fp_dpks_descriptor)
    
    return c


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


def get_abacus_STRU(STRU, INPUT=None, n_ele=None):
    # read in geometry from STRU file. n_ele is the number of elements.
    # Either n_ele or INPUT should be provided.
    data = get_frame_from_stru(STRU)
    data["atom_masses"] = data.pop("masses")
    data["cells"] = data.pop("cells")[0]
    data["coords"] = data.pop("coords")[0]
    if "orb_files" not in data:
        data["orb_files"] = None
    if "dpks_descriptor" not in data:
        data["dpks_descriptor"] = None
    return data


def make_supercell_abacus(from_struct, super_cell):
    to_struct = copy.deepcopy(from_struct)

    if "atom_types" in from_struct:
        new_types = []
        # to_struct["atom_types"] = (
        #    from_struct["atom_types"] * super_cell[0] * super_cell[1] * super_cell[2]
        # )
        for idx_atm in from_struct["atom_types"]:
            new_types += [idx_atm] * super_cell[0] * super_cell[1] * super_cell[2]
        to_struct["atom_types"] = new_types
    to_atom_num = (
        sum(from_struct["atom_numbs"]) * super_cell[0] * super_cell[1] * super_cell[2]
    )
    new_coord = np.zeros((to_atom_num, 3))
    idx_atm = 0
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
                    new_coord[idx_atm] = coord
                    idx_atm += 1

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
    if not isinstance(kspacing, list):
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

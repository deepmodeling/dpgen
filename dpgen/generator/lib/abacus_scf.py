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
            ret += "ecutwfc {:f}\n".format(fp_params["ecutwfc"])
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
            ], (
                "'kspacing' only accept a float, or a list of one or three float, or a string of one or three float"
            )
            ret += "kspacing "
            for ikspacing in fp_params["kspacing"]:
                assert ikspacing >= 0, "'kspacing' should be non-negative."
                ret += f"{ikspacing:f} "
            ret += "\n"
        elif key == "scf_thr":
            fp_params["scf_thr"] = float(fp_params["scf_thr"])
            ret += "scf_thr {:e}\n".format(fp_params["scf_thr"])
        elif key == "scf_nmax":
            fp_params["scf_nmax"] = int(fp_params["scf_nmax"])
            assert fp_params["scf_nmax"] >= 0 and isinstance(
                fp_params["scf_nmax"], int
            ), "'scf_nmax' should be a positive integer."
            ret += "scf_nmax %d\n" % fp_params["scf_nmax"]  # noqa: UP031
        elif key == "basis_type":
            assert fp_params["basis_type"] in [
                "pw",
                "lcao",
                "lcao_in_pw",
            ], "'basis_type' must in 'pw', 'lcao' or 'lcao_in_pw'."
            ret += "basis_type {}\n".format(fp_params["basis_type"])
        elif key == "dft_functional":
            ret += "dft_functional {}\n".format(fp_params["dft_functional"])
        elif key == "gamma_only":
            if isinstance(fp_params["gamma_only"], str):
                fp_params["gamma_only"] = int(eval(fp_params["gamma_only"]))
            assert fp_params["gamma_only"] == 0 or fp_params["gamma_only"] == 1, (
                "'gamma_only' should be either 0 or 1."
            )
            ret += "gamma_only %d\n" % fp_params["gamma_only"]  # noqa: UP031
        elif key == "mixing_type":
            assert fp_params["mixing_type"] in [
                "plain",
                "kerker",
                "pulay",
                "pulay-kerker",
                "broyden",
            ]
            ret += "mixing_type {}\n".format(fp_params["mixing_type"])
        elif key == "mixing_beta":
            fp_params["mixing_beta"] = float(fp_params["mixing_beta"])
            assert fp_params["mixing_beta"] >= 0 and fp_params["mixing_beta"] < 1, (
                "'mixing_beta' should between 0 and 1."
            )
            ret += "mixing_beta {:f}\n".format(fp_params["mixing_beta"])
        elif key == "symmetry":
            if isinstance(fp_params["symmetry"], str):
                fp_params["symmetry"] = int(eval(fp_params["symmetry"]))
            assert fp_params["symmetry"] == 0 or fp_params["symmetry"] == 1, (
                "'symmetry' should be either 0 or 1."
            )
            ret += "symmetry %d\n" % fp_params["symmetry"]  # noqa: UP031
        elif key == "nbands":
            fp_params["nbands"] = int(fp_params["nbands"])
            assert fp_params["nbands"] > 0 and isinstance(fp_params["nbands"], int), (
                "'nbands' should be a positive integer."
            )
            ret += "nbands %d\n" % fp_params["nbands"]  # noqa: UP031
        elif key == "nspin":
            fp_params["nspin"] = int(fp_params["nspin"])
            assert (
                fp_params["nspin"] == 1
                or fp_params["nspin"] == 2
                or fp_params["nspin"] == 4
            ), "'nspin' can anly take 1, 2 or 4"
            ret += "nspin %d\n" % fp_params["nspin"]  # noqa: UP031
        elif key == "ks_solver":
            assert fp_params["ks_solver"] in [
                "cg",
                "dav",
                "lapack",
                "genelpa",
                "hpseps",
                "scalapack_gvx",
            ], (
                "'ks_sover' should in 'cgx', 'dav', 'lapack', 'genelpa', 'hpseps', 'scalapack_gvx'."
            )
            ret += "ks_solver {}\n".format(fp_params["ks_solver"])
        elif key == "smearing_method":
            assert fp_params["smearing_method"] in [
                "gauss",
                "gaussian",
                "fd",
                "fixed",
                "mp",
                "mp2",
                "mv",
            ], (
                "'smearing_method' should in 'gauss', 'gaussian', 'fd', 'fixed', 'mp', 'mp2', 'mv'. "
            )
            ret += "smearing_method {}\n".format(fp_params["smearing_method"])
        elif key == "smearing_sigma":
            fp_params["smearing_sigma"] = float(fp_params["smearing_sigma"])
            assert fp_params["smearing_sigma"] >= 0, (
                "'smearing_sigma' should be non-negative."
            )
            ret += "smearing_sigma {:f}\n".format(fp_params["smearing_sigma"])
        elif key == "cal_force":
            if isinstance(fp_params["cal_force"], str):
                fp_params["cal_force"] = int(eval(fp_params["cal_force"]))
            assert fp_params["cal_force"] == 0 or fp_params["cal_force"] == 1, (
                "'cal_force' should be either 0 or 1."
            )
            ret += "cal_force %d\n" % fp_params["cal_force"]  # noqa: UP031
        elif key == "cal_stress":
            if isinstance(fp_params["cal_stress"], str):
                fp_params["cal_stress"] = int(eval(fp_params["cal_stress"]))
            assert fp_params["cal_stress"] == 0 or fp_params["cal_stress"] == 1, (
                "'cal_stress' should be either 0 or 1."
            )
            ret += "cal_stress %d\n" % fp_params["cal_stress"]  # noqa: UP031
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
            ret += "deepks_out_labels %d\n" % fp_params["deepks_out_labels"]  # noqa: UP031
        elif key == "deepks_descriptor_lmax":
            fp_params["deepks_descriptor_lmax"] = int(
                fp_params["deepks_descriptor_lmax"]
            )
            assert fp_params["deepks_descriptor_lmax"] >= 0, (
                "'deepks_descriptor_lmax' should be  a positive integer."
            )
            ret += "deepks_descriptor_lmax %d\n" % fp_params["deepks_descriptor_lmax"]  # noqa: UP031
        elif key == "deepks_scf":
            if isinstance(fp_params["deepks_scf"], str):
                fp_params["deepks_scf"] = int(eval(fp_params["deepks_scf"]))
            assert fp_params["deepks_scf"] == 0 or fp_params["deepks_scf"] == 1, (
                "'deepks_scf' should be either 0 or 1."
            )
            ret += "deepks_scf %d\n" % fp_params["deepks_scf"]  # noqa: UP031
        elif key == "deepks_model":
            ret += "deepks_model {}\n".format(
                os.path.join(
                    extra_file_path, os.path.split(fp_params["deepks_model"])[1]
                )
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
    sys_data_copy = copy.deepcopy(sys_data)
    # re-construct the path of files by pporb + file name
    # when element in sys_data is part of type_map/fp_pp_files
    # we need to only pass the pp_file in sys_data, but not all pp_files
    if type_map is None:
        type_map = sys_data_copy["atom_names"]

    missing_atoms = set(sys_data_copy["atom_names"]) - set(type_map)
    if len(missing_atoms) > 0:
        raise ValueError(
            f"Some atoms in sys_data are not in type_map: {missing_atoms}. "
            "Please provide a valid type_map."
        )

    if len(fp_pp_files) != len(type_map):
        raise ValueError(
            "The length of fp_pp_files should be equal to the length of type_map."
        )
    if fp_orb_files is not None and len(fp_orb_files) != len(type_map):
        raise ValueError(
            "The length of fp_orb_files should be equal to the length of type_map."
        )

    fp_pp_files = [
        os.path.join(pporb, fp_pp_files[type_map.index(atom_name)])
        for atom_name in sys_data_copy["atom_names"]
    ]

    if fp_orb_files is not None:
        fp_orb_files = [
            os.path.join(pporb, fp_orb_files[type_map.index(atom_name)])
            for atom_name in sys_data_copy["atom_names"]
        ]

    if fp_dpks_descriptor is not None:
        fp_dpks_descriptor = os.path.join(pporb, fp_dpks_descriptor)

    # we need to make sure that the shape of cells and coords are the same
    # and if they are 2D, we need to convert them to 3D
    cells = np.array(sys_data["cells"])
    coords = np.array(sys_data["coords"])
    assert len(cells.shape) == len(coords.shape), (
        "cells and coords should have the same shape."
    )

    if len(cells.shape) == 2:
        sys_data_copy["cells"] = np.array([cells])
        sys_data_copy["coords"] = np.array([coords])
    c = make_unlabeled_stru(
        sys_data_copy,
        0,
        pp_file=fp_pp_files,
        numerical_orbital=fp_orb_files,
        numerical_descriptor=fp_dpks_descriptor,
    )

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


def get_abacus_STRU(STRU):
    """Read STRU file and return a dictionary containing the structure information.

    Args:
        STRU (str): The path of STRU file.

    Returns
    -------
    dict: A dictionary containing the structure information.
    {
        "atom_names": list of str,
        "atom_numbs": list of int,
        "atom_masses": list of float,
        "coords": np.ndarray,
        "cells": np.ndarray,
        "pp_files": list of str,
        "orb_files": list of str,
        "dpks_descriptor": str,
    }
    """
    data = get_frame_from_stru(STRU)
    data["atom_masses"] = data.pop("masses")
    data["cells"] = data.pop("cells")[0]
    data["coords"] = data.pop("coords")[0]
    assert "pp_files" in data, "pp_files should be provided in STRU file."
    if None in data["pp_files"]:
        data["pp_files"] = None
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
        to_struct["atom_types"] = np.array(new_types)

    # expand move, spins
    for ikey in ["move", "spins"]:
        if ikey in from_struct:
            new_list = []
            for ia in range(sum(from_struct["atom_numbs"])):
                new_list += (
                    [from_struct[ikey][0][ia]]
                    * super_cell[0]
                    * super_cell[1]
                    * super_cell[2]
                )
            to_struct[ikey] = np.array([new_list])

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

import numpy as np

atomic_numbers = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Nh": 113,
    "Fl": 114,
    "Mc": 115,
    "Lv": 116,
    "Ts": 117,
    "Og": 118,
}


default_config = {
    "GLOBAL": {"PROJECT": "DPGEN"},
    "FORCE_EVAL": {
        "METHOD": "QS",
        "STRESS_TENSOR": "ANALYTICAL",
        "DFT": {
            "BASIS_SET_FILE_NAME": "./cp2k_basis_pp_file/BASIS_MOLOPT",
            "POTENTIAL_FILE_NAME": "./cp2k_basis_pp_file/GTH_POTENTIALS",
            "CHARGE": 0,
            "UKS": "F",
            "MULTIPLICITY": 1,
            "MGRID": {"CUTOFF": 400, "REL_CUTOFF": 50, "NGRIDS": 4},
            "QS": {"EPS_DEFAULT": "1.0E-12"},
            "SCF": {"SCF_GUESS": "ATOMIC", "EPS_SCF": "1.0E-6", "MAX_SCF": 50},
            "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
        },
        "SUBSYS": {
            "CELL": {"A": "10 .0 .0", "B": ".0 10 .0", "C": ".0 .0 10"},
            "COORD": {"@include": "coord.xyz"},
            "KIND": {
                "_": ["H", "C", "N"],
                "POTENTIAL": ["GTH-PBE-q1", "GTH-PBE-q4", "GTH-PBE-q5"],
                "BASIS_SET": ["DZVP-MOLOPT-GTH", "DZVP-MOLOPT-GTH", "DZVP-MOLOPT-GTH"],
            },
        },
        "PRINT": {"FORCES": {"_": "ON"}, "STRESS_TENSOR": {"_": "ON"}},
    },
}


def update_dict(old_d, update_d):
    """A method to recursive update dict
    :old_d: old dictionary
    :update_d: some update value written in dictionary form.
    """
    import collections.abc

    for k, v in update_d.items():
        if (
            k in old_d
            and isinstance(old_d[k], dict)
            and isinstance(update_d[k], collections.abc.Mapping)
        ):
            update_dict(old_d[k], update_d[k])
        else:
            old_d[k] = update_d[k]


def iterdict(d, out_list, flag=None, indent=0):
    """
    :doc: a recursive expansion of dictionary into cp2k input
    :k: current key
    :v: current value
    :d: current dictionary under expansion
    :flag: used to record dictionary state. if flag is None,
    it means we are in top level dict. flag is a string.
    :indent: intent for current section.
    """
    for k, v in d.items():
        k = str(k)  # cast key into string
        # if value is dictionary
        if isinstance(v, dict):
            # flag == None, it is now in top level section of cp2k
            if flag is None:
                out_list.append("&" + k)
                out_list.append("&END " + k)
                iterdict(v, out_list, k, indent + 2)
            # flag is not None, now it has name of section
            else:
                index = out_list.index(" " * (indent - 2) + "&END " + flag)
                out_list.insert(index, " " * indent + "&" + k + " #" + flag)
                out_list.insert(index + 1, " " * indent + "&END " + k + " #" + flag)
                # the flag now contains its parent section name, separed by "#".
                iterdict(v, out_list, k + " #" + flag, indent + 2)
        elif isinstance(v, list):
            #            print("we have encountered the repeat section!")
            index = out_list.index(" " * (indent - 2) + "&" + flag)
            # delete the current constructed repeat section
            del out_list[index : index + 2]
            # do a loop over key and corresponding list
            k_tmp_list = []
            v_list_tmp_list = []
            for k_tmp, v_tmp in d.items():
                k_tmp_list.append(str(k_tmp))
                v_list_tmp_list.append(v_tmp)
            for repeat_keyword in zip(*v_list_tmp_list):
                out_list.insert(index, " " * (indent - 2) + "&" + flag)
                out_list.insert(index + 1, " " * (indent - 2) + "&END " + flag)
                for idx, k_tmp in enumerate(k_tmp_list):
                    if k_tmp == "_":
                        out_list[index] = (
                            " " * (indent - 2)
                            + "&"
                            + flag.split(" #")[0]
                            + " "
                            + repeat_keyword[idx]
                        )
                    else:
                        out_list.insert(
                            index + 1,
                            " " * (indent) + k_tmp + " " + repeat_keyword[idx],
                        )
            break

        else:
            v = str(v)
            if flag is None:
                out_list.append(k + " " + v)
                print(k, ":", v)
            else:
                if k == "_":
                    index = out_list.index(" " * (indent - 2) + "&" + flag)
                    out_list[index] = (
                        " " * (indent - 2) + "&" + flag.split(" #")[0] + " " + v
                    )

                else:
                    index = out_list.index(" " * (indent - 2) + "&END " + flag)
                    out_list.insert(index, " " * indent + k + " " + v)


def calculate_multiplicity(atom_names, atom_types, charge=0):
    """
    Calculate the multiplicity based on atom species, quantities, and system charge.

    This function provides a basic heuristic for determining multiplicity:
    - Even number of electrons -> singlet (multiplicity = 1)
    - Odd number of electrons -> doublet (multiplicity = 2)

    Note: This approach assumes that an odd electron count always results in a doublet state.
    It does not account for systems with multiple unpaired electrons, which can have higher
    multiplicities (e.g., triplet, quartet, etc.). Users should be aware of this limitation
    and use the function accordingly.

    :param atom_names: List of element symbols.
    :param atom_types: List of atom type indices.
    :param charge: System charge (default: 0).
    :return: Multiplicity.
    """
    # Calculate the total number of electrons
    total_electrons = 0
    for idx in atom_types:
        element = atom_names[idx]
        try:
            total_electrons += atomic_numbers[element]
        except KeyError:
            raise ValueError(f"Unknown element '{element}' encountered in atom_names.")

    # Subtract/add electrons based on system charge
    # Positive charge means we remove electrons, negative charge means we add electrons
    total_electrons -= charge

    # Determine multiplicity based on the total number of electrons
    # Even number of electrons -> singlet (multiplicity = 1)
    # Odd number of electrons -> doublet (multiplicity = 2)
    if total_electrons % 2 == 0:
        return 1
    else:
        return 2


def make_cp2k_input(sys_data, fp_params):
    # covert cell to cell string
    cell = sys_data["cells"][0]
    cell = np.reshape(cell, [3, 3])
    cell_a = np.array2string(cell[0, :])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1, :])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2, :])
    cell_c = cell_c[1:-1]

    atom_names = sys_data["atom_names"]
    atom_types = sys_data["atom_types"]
    # Get system charge if provided, default to 0
    charge = sys_data.get("charge", 0)
    dft_params = fp_params.get("FORCE_EVAL", {}).get("DFT", {})
    if "MULTIPLICITY" in dft_params:
        multiplicity = dft_params["MULTIPLICITY"]
    else:
        multiplicity = calculate_multiplicity(atom_names, atom_types, charge)

    # get update from user
    user_config = fp_params
    # get update from cell
    cell_config = {
        "FORCE_EVAL": {"SUBSYS": {"CELL": {"A": cell_a, "B": cell_b, "C": cell_c}}}
    }
    # get update for multiplicity
    multiplicity_config = {"FORCE_EVAL": {"DFT": {"MULTIPLICITY": multiplicity}}}
    update_dict(default_config, user_config)
    update_dict(default_config, cell_config)
    update_dict(default_config, multiplicity_config)
    # output list
    input_str = []
    iterdict(default_config, input_str)
    string = "\n".join(input_str)
    return string


def make_cp2k_xyz(sys_data):
    # get structral information
    atom_names = sys_data["atom_names"]
    atom_types = sys_data["atom_types"]

    # write coordinate to xyz file used by cp2k input
    coord_list = sys_data["coords"][0]
    u = np.array(atom_names)
    atom_list = u[atom_types]
    x = "\n"
    for kind, coord in zip(atom_list, coord_list):
        x += str(kind) + " " + str(coord[:])[1:-1] + "\n"
    return x


def make_cp2k_input_from_external(sys_data, exinput_path):
    # read the input content as string
    with open(exinput_path) as f:
        exinput = f.readlines()

    # find the ABC cell string
    for line_idx, line in enumerate(exinput):
        if "ABC" in line:
            delete_cell_idx = line_idx
            delete_cell_line = line

    # remove the useless CELL line
    exinput.remove(delete_cell_line)

    # insert the cell information
    # covert cell to cell string
    cell = sys_data["cells"][0]
    cell = np.reshape(cell, [3, 3])
    cell_a = np.array2string(cell[0, :])
    cell_a = cell_a[1:-1]
    cell_b = np.array2string(cell[1, :])
    cell_b = cell_b[1:-1]
    cell_c = np.array2string(cell[2, :])
    cell_c = cell_c[1:-1]

    exinput.insert(delete_cell_idx, "A  " + cell_a + "\n")
    exinput.insert(delete_cell_idx + 1, "B  " + cell_b + "\n")
    exinput.insert(delete_cell_idx + 2, "C  " + cell_c + "\n")

    return "".join(exinput)

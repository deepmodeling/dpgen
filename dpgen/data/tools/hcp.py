import numpy as np


def numb_atoms():
    return 2


def gen_box():
    box = np.array(
        [[1, 0, 0], [0.5, 0.5 * np.sqrt(3), 0], [0, 0, 2.0 * np.sqrt(2.0 / 3.0)]]
    )
    return box


def poscar_unit(latt):
    box = gen_box()
    ret = ""
    ret += f"HCP : a = {latt:f} / sqrt(2)\n"
    ret += "%.16f\n" % (latt / np.sqrt(2))
    ret += f"{box[0][0]:.16f} {box[0][1]:.16f} {box[0][2]:.16f}\n"
    ret += f"{box[1][0]:.16f} {box[1][1]:.16f} {box[1][2]:.16f}\n"
    ret += f"{box[2][0]:.16f} {box[2][1]:.16f} {box[2][2]:.16f}\n"
    ret += "X\n"
    ret += "%d\n" % numb_atoms()  # noqa: UP031
    ret += "Direct\n"
    ret += f"{0:.16f} {0:.16f} {0:.16f}\n"
    ret += f"{1.0 / 3:.16f} {1.0 / 3:.16f} {1.0 / 2:.16f}\n"
    return ret

import numpy as np


def numb_atoms():
    return 2


def gen_box():
    box = [
        [0.000000, 1.000000, 1.000000],
        [1.000000, 0.000000, 1.000000],
        [1.000000, 1.000000, 0.000000],
    ]
    return np.array(box)


def poscar_unit(latt):
    box = gen_box()
    ret = ""
    ret += "DIAMOND\n"
    ret += "%.16f\n" % (latt)
    ret += f"{box[0][0]:.16f} {box[0][1]:.16f} {box[0][2]:.16f}\n"
    ret += f"{box[1][0]:.16f} {box[1][1]:.16f} {box[1][2]:.16f}\n"
    ret += f"{box[2][0]:.16f} {box[2][1]:.16f} {box[2][2]:.16f}\n"
    ret += "Type\n"
    ret += "%d\n" % numb_atoms()
    ret += "Direct\n"
    ret += "{:.16f} {:.16f} {:.16f}\n".format(
        0.12500000000000,
        0.12500000000000,
        0.12500000000000,
    )
    ret += "{:.16f} {:.16f} {:.16f}\n".format(
        0.87500000000000,
        0.87500000000000,
        0.87500000000000,
    )
    return ret

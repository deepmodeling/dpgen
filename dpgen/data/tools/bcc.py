import numpy as np


def numb_atoms():
    return 2


def gen_box():
    return np.eye(3)


def poscar_unit(latt):
    box = gen_box()
    ret = ""
    ret += f"BCC : a = {latt:f} \n"
    ret += f"{latt:.16f}\n"
    ret += f"{box[0][0]:.16f} {box[0][1]:.16f} {box[0][2]:.16f}\n"
    ret += f"{box[1][0]:.16f} {box[1][1]:.16f} {box[1][2]:.16f}\n"
    ret += f"{box[2][0]:.16f} {box[2][1]:.16f} {box[2][2]:.16f}\n"
    ret += "X\n"
    ret += "%d\n" % numb_atoms()
    ret += "Direct\n"
    ret += f"{0.0:.16f} {0.0:.16f} {0.0:.16f}\n"
    ret += f"{0.5:.16f} {0.5:.16f} {0.5:.16f}\n"
    return ret

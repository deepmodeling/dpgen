import numpy as np

def voigt_to_stress(inpt) :
    ret = np.zeros((3,3))
    ret[0][0] = inpt[0]
    ret[1][1] = inpt[1]
    ret[2][2] = inpt[2]
    ret[0][1] = inpt[3]
    ret[0][2] = inpt[4]
    ret[1][2] = inpt[5]
    ret[2][0] = ret[0][2]
    ret[1][0] = ret[0][1]
    ret[2][1] = ret[1][2]
    return ret

import numpy as np

def numb_atoms () :
    return 2

def gen_box () :    
    return np.eye(3)

def poscar_unit (latt) :
    box = gen_box()
    ret  = ""
    ret += "BCC : a = %f \n" % latt
    ret += "%.16f\n" % (latt)
    ret += "%.16f %.16f %.16f\n" % (box[0][0], box[0][1], box[0][2])
    ret += "%.16f %.16f %.16f\n" % (box[1][0], box[1][1], box[1][2])
    ret += "%.16f %.16f %.16f\n" % (box[2][0], box[2][1], box[2][2])
    ret += "Type\n"
    ret += "%d\n" % numb_atoms()
    ret += "Direct\n"
    ret += "%.16f %.16f %.16f\n" % (0.0, 0.0, 0.0)
    ret += "%.16f %.16f %.16f\n" % (0.5, 0.5, 0.5)
    return ret

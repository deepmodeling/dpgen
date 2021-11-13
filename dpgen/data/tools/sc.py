import numpy as np

def numb_atoms () :
    return 1

def gen_box () :    
    return np.eye(3)

def poscar_unit (latt) :
    box = gen_box()
    ret  = ""
    ret += "SC : a = %f \n" % latt
    ret += "%.16f\n" % (latt)
    ret += "%.16f %.16f %.16f\n" % (box[0][0], box[0][1], box[0][2])
    ret += "%.16f %.16f %.16f\n" % (box[1][0], box[1][1], box[1][2])
    ret += "%.16f %.16f %.16f\n" % (box[2][0], box[2][1], box[2][2])
    ret += "Type\n"
    ret += "%d\n" % numb_atoms()
    ret += "Direct\n"
    ret += "%.16f %.16f %.16f\n" % (0.0, 0.0, 0.0)
    return ret

def STRU_unit(latt):

    box = gen_box() * latt
    stru_data = {}
    stru_data['cells'] = box
    stru_data['atom_numbs'] = [1]
    stru_data['coords'] = np.array([0.0, 0.0, 0.0])
    stru_data['types'] = [0]
    return stru_data
import numpy as np

def numb_atoms () :
    return 2

def gen_box () :    
    box = [[0.000000, 1.000000, 1.000000],
           [1.000000, 0.000000, 1.000000],
           [1.000000, 1.000000, 0.000000]
    ]
    return np.array(box)

def poscar_unit (latt) :
    box = gen_box()
    ret  = ""
    ret += "DIAMOND\n"
    ret += "%.16f\n" % (latt)
    ret += "%.16f %.16f %.16f\n" % (box[0][0], box[0][1], box[0][2])
    ret += "%.16f %.16f %.16f\n" % (box[1][0], box[1][1], box[1][2])
    ret += "%.16f %.16f %.16f\n" % (box[2][0], box[2][1], box[2][2])
    ret += "Type\n"
    ret += "%d\n" % numb_atoms()
    ret += "Direct\n"
    ret += "%.16f %.16f %.16f\n" % (0.12500000000000,   0.12500000000000,   0.12500000000000)
    ret += "%.16f %.16f %.16f\n" % (0.87500000000000,   0.87500000000000,   0.87500000000000)
    return ret

def STRU_unit(latt):

    box = gen_box() * latt
    stru_data = {}
    stru_data['cells'] = box
    stru_data['atom_numbs'] = [2]
    stru_data['coords'] = np.matmul(np.array([[0.125, 0.125, 0.125], [0.875, 0.875, 0.875]]), box)
    stru_data['types'] = [0, 0]
    return stru_data